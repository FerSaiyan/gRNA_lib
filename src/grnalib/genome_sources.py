from __future__ import annotations

import os
import shutil
import tempfile
import urllib.parse
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .sequence import normalize_dna

NCBI_DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
ENSEMBL_REST_BASE = "https://rest.ensembl.org"
USER_AGENT = "gRNA-Library/0.2 (+https://github.com/FerSaiyan/gRNA_lib)"


@dataclass(frozen=True)
class FastaRecordInfo:
    name: str
    description: str
    length: int

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "description": self.description,
            "length": self.length,
        }


@dataclass(frozen=True)
class GenomeDownload:
    source: str
    accession: str
    fasta_path: Path
    records: tuple[FastaRecordInfo, ...]
    requested_chromosomes: tuple[str, ...] = ()

    def to_dict(self) -> dict:
        return {
            "source": self.source,
            "accession": self.accession,
            "fasta_path": str(self.fasta_path),
            "records": [record.to_dict() for record in self.records],
            "requested_chromosomes": list(self.requested_chromosomes),
        }


def default_cache_dir() -> Path:
    override = os.environ.get("GRNALIB_CACHE_DIR")
    if override:
        return Path(override).expanduser()
    return Path.home() / ".cache" / "grna-lib"


def fasta_records(path: str | Path) -> tuple[FastaRecordInfo, ...]:
    """Return FASTA record names and lengths without loading sequences into memory."""
    path = Path(path)
    records: list[FastaRecordInfo] = []
    name: str | None = None
    description = ""
    length = 0

    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append(FastaRecordInfo(name, description, length))
                description = line[1:].strip()
                if not description:
                    raise ValueError(f"FASTA header without a record name in {path}")
                name = description.split()[0]
                length = 0
            else:
                if name is None:
                    raise ValueError(f"FASTA sequence appears before the first header in {path}")
                length += len("".join(line.split()))

    if name is not None:
        records.append(FastaRecordInfo(name, description, length))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return tuple(records)


def read_fasta_record(path: str | Path, record: str | None = None) -> tuple[str, str]:
    """Read one FASTA record and return its name and normalized sequence.

    If the file contains more than one record, callers must select one by its
    first-header token. This avoids silently joining chromosome boundaries.
    """
    path = Path(path)
    infos = fasta_records(path)
    if record is None:
        if len(infos) != 1:
            names = ", ".join(info.name for info in infos[:8])
            suffix = "..." if len(infos) > 8 else ""
            raise ValueError(
                f"{path} contains {len(infos)} FASTA records; choose one record "
                f"({names}{suffix}) rather than concatenating them"
            )
        record = infos[0].name

    selected: list[str] = []
    active = False
    found = False
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].strip().split()[0]
                if active:
                    break
                active = current == record
                if active:
                    found = True
                continue
            if active:
                selected.append(line)

    if not found:
        raise ValueError(f"FASTA record {record!r} was not found in {path}")
    return record, normalize_dna("".join(selected))


def parse_fasta_text(text: str) -> list[tuple[str, str]]:
    """Parse FASTA text into normalized name/sequence records."""
    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, normalize_dna("".join(parts))))
            description = line[1:].strip()
            if not description:
                raise ValueError("FASTA header without a record name")
            name = description.split()[0]
            parts = []
        else:
            if name is None:
                raise ValueError("FASTA sequence appears before the first header")
            parts.append(line)
    if name is not None:
        records.append((name, normalize_dna("".join(parts))))
    if not records:
        raise ValueError("No FASTA records found")
    return records


def _download_to_file(url: str, destination: Path, *, timeout: int) -> None:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "application/zip,application/octet-stream,*/*",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response, destination.open("wb") as out:
            shutil.copyfileobj(response, out)
    except Exception as exc:
        raise RuntimeError(f"download failed from {url}: {exc}") from exc


def _find_genomic_fasta(dataset_zip: Path, destination: Path) -> Path:
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(dataset_zip) as archive:
        candidates = [
            member for member in archive.infolist()
            if not member.is_dir()
            and member.filename.startswith("ncbi_dataset/data/")
            and member.filename.endswith(("_genomic.fna", "genomic.fna"))
        ]
        if not candidates:
            raise RuntimeError("NCBI data package did not contain a genomic FASTA")
        if len(candidates) > 1:
            raise RuntimeError(
                "NCBI data package contained multiple genomic FASTA files; "
                "request one assembly accession at a time"
            )
        member = candidates[0]
        output = destination / Path(member.filename).name
        with archive.open(member) as src, output.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    return output


def download_ncbi_genome(
    accession: str,
    *,
    chromosomes: Iterable[str] | None = None,
    cache_dir: str | Path | None = None,
    timeout: int = 300,
    force: bool = False,
) -> GenomeDownload:
    """Download an NCBI Datasets genome FASTA by GCF/GCA assembly accession.

    Chromosome labels such as 1, X or MT are sent to NCBI so filtering occurs
    server-side instead of requiring a whole-assembly download first.
    """
    accession = accession.strip().upper()
    if not (accession.startswith("GCF_") or accession.startswith("GCA_")):
        raise ValueError("NCBI assembly accession must start with GCF_ or GCA_")

    chromosome_list = tuple(
        chromosome.strip() for chromosome in (chromosomes or ()) if chromosome.strip()
    )
    cache_root = Path(cache_dir) if cache_dir is not None else default_cache_dir()
    key = accession + (
        "__" + "-".join(chromosome_list).replace("/", "_")
        if chromosome_list else "__all"
    )
    target_dir = cache_root.expanduser() / "ncbi" / key

    existing = sorted(target_dir.glob("*_genomic.fna")) + sorted(target_dir.glob("genomic.fna"))
    if existing and not force:
        fasta_path = existing[0]
        return GenomeDownload(
            source="NCBI Datasets",
            accession=accession,
            fasta_path=fasta_path,
            records=fasta_records(fasta_path),
            requested_chromosomes=chromosome_list,
        )

    params: list[tuple[str, str]] = [
        ("include_annotation_type", "GENOME_FASTA"),
        ("hydrated", "FULLY_HYDRATED"),
    ]
    params.extend(("chromosomes", chromosome) for chromosome in chromosome_list)
    query = urllib.parse.urlencode(params)
    url = (
        f"{NCBI_DATASETS_BASE}/genome/accession/"
        f"{urllib.parse.quote(accession, safe='')}/download?{query}"
    )

    target_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="grnalib-ncbi-") as tmp:
        archive_path = Path(tmp) / "dataset.zip"
        _download_to_file(url, archive_path, timeout=timeout)
        try:
            fasta_path = _find_genomic_fasta(archive_path, target_dir)
        except zipfile.BadZipFile as exc:
            raise RuntimeError(
                "NCBI returned a response that was not a valid genome data package"
            ) from exc

    return GenomeDownload(
        source="NCBI Datasets",
        accession=accession,
        fasta_path=fasta_path,
        records=fasta_records(fasta_path),
        requested_chromosomes=chromosome_list,
    )


def fetch_ensembl_region(
    species: str,
    region: str,
    *,
    timeout: int = 60,
) -> str:
    """Fetch a genomic region from Ensembl REST (currently limited to 10 Mb/request)."""
    species = species.strip()
    region = region.strip()
    if not species or not region:
        raise ValueError("species and region are required")

    url = (
        f"{ENSEMBL_REST_BASE}/sequence/region/"
        f"{urllib.parse.quote(species, safe='')}/"
        f"{urllib.parse.quote(region, safe=':.-')}"
        "?content-type=text/plain"
    )
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "text/plain",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            text = response.read().decode("utf-8").strip()
    except Exception as exc:
        raise RuntimeError(f"Ensembl sequence request failed: {exc}") from exc

    if not text:
        raise RuntimeError("Ensembl returned an empty sequence")
    if text.startswith(">"):
        parsed = parse_fasta_text(text)
        if len(parsed) != 1:
            raise RuntimeError("Ensembl unexpectedly returned multiple FASTA records")
        return parsed[0][1]
    return normalize_dna(text)
