from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .backends.crisprware import build_crisprware_index
from .genome_sources import default_cache_dir, download_ncbi_genome


@dataclass(frozen=True)
class ReferenceIndexResult:
    accession: str
    chromosomes: tuple[str, ...]
    fasta_path: Path
    index_prefix: Path
    pam: str
    spacer_length: int
    pam_side: str
    cached: bool

    def to_dict(self) -> dict:
        return {
            "accession": self.accession,
            "chromosomes": list(self.chromosomes),
            "fasta_path": str(self.fasta_path),
            "index_prefix": str(self.index_prefix),
            "index_file": str(self.index_prefix) + ".crot",
            "pam": self.pam,
            "spacer_length": self.spacer_length,
            "pam_side": self.pam_side,
            "cached": self.cached,
        }


def _safe_component(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "-", value.strip())
    return value.strip("-") or "value"


def ncbi_index_paths(
    accession: str,
    *,
    chromosomes: Iterable[str] | None = None,
    pam: str = "NGG",
    spacer_length: int = 20,
    pam_side: str = "3prime",
    cache_dir: str | Path | None = None,
    fasta_name: str | None = None,
) -> tuple[Path, Path | None]:
    """Return the index output directory and expected crispr-ots prefix."""
    chromosome_list = tuple(
        chromosome.strip() for chromosome in (chromosomes or ()) if chromosome.strip()
    )
    cache_root = Path(cache_dir) if cache_dir is not None else default_cache_dir()
    chr_key = "-".join(_safe_component(item) for item in chromosome_list) or "all"
    key = "__".join(
        [
            _safe_component(accession.upper()),
            chr_key,
            _safe_component(pam.upper()),
            str(spacer_length),
            _safe_component(pam_side),
        ]
    )
    output_dir = cache_root.expanduser() / "indexes" / key
    if fasta_name is None:
        return output_dir, None
    stem = Path(fasta_name).stem
    prefix = output_dir / f"{stem}_crisprots" / f"{stem}_crisprots"
    return output_dir, prefix


def build_or_reuse_ncbi_index(
    accession: str,
    *,
    chromosomes: Iterable[str] | None = None,
    pam: str = "NGG",
    spacer_length: int = 20,
    pam_side: str = "3prime",
    cache_dir: str | Path | None = None,
    crisprware_executable: str = "crisprware",
    bin_width: int | None = None,
    force: bool = False,
) -> ReferenceIndexResult:
    """Download/cache an NCBI FASTA and build or reuse its crispr-ots index."""
    if pam_side not in {"3prime", "5prime"}:
        raise ValueError("pam_side must be 3prime or 5prime")
    chromosome_list = tuple(
        chromosome.strip() for chromosome in (chromosomes or ()) if chromosome.strip()
    )
    download = download_ncbi_genome(
        accession,
        chromosomes=chromosome_list,
        cache_dir=cache_dir,
        force=False,
    )
    output_dir, prefix = ncbi_index_paths(
        accession,
        chromosomes=chromosome_list,
        pam=pam,
        spacer_length=spacer_length,
        pam_side=pam_side,
        cache_dir=cache_dir,
        fasta_name=download.fasta_path.name,
    )
    assert prefix is not None
    index_file = Path(str(prefix) + ".crot")
    manifest = output_dir / "grnalib-index.json"

    if index_file.exists() and manifest.exists() and not force:
        return ReferenceIndexResult(
            accession=accession.upper(),
            chromosomes=chromosome_list,
            fasta_path=download.fasta_path,
            index_prefix=prefix,
            pam=pam.upper(),
            spacer_length=spacer_length,
            pam_side=pam_side,
            cached=True,
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    build_crisprware_index(
        download.fasta_path,
        pam=pam,
        spacer_length=spacer_length,
        pam_5_prime=pam_side == "5prime",
        bin_width=bin_width,
        output_directory=output_dir,
        executable=crisprware_executable,
    )
    if not index_file.exists():
        raise RuntimeError(
            f"CRISPRware returned success but the expected index file was not found: {index_file}"
        )
    manifest.write_text(
        json.dumps(
            {
                "accession": accession.upper(),
                "chromosomes": list(chromosome_list),
                "fasta_path": str(download.fasta_path),
                "index_prefix": str(prefix),
                "pam": pam.upper(),
                "spacer_length": spacer_length,
                "pam_side": pam_side,
                "bin_width": bin_width,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return ReferenceIndexResult(
        accession=accession.upper(),
        chromosomes=chromosome_list,
        fasta_path=download.fasta_path,
        index_prefix=prefix,
        pam=pam.upper(),
        spacer_length=spacer_length,
        pam_side=pam_side,
        cached=False,
    )
