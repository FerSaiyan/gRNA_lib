from __future__ import annotations

import json
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
class NcbiAssembly:
    accession: str
    organism_name: str
    tax_id: int | None
    assembly_name: str
    assembly_level: str
    source_database: str
    refseq_category: str
    release_date: str
    submitter: str
    description: str
    synonym: str

    @property
    def is_reference(self) -> bool:
        return "reference genome" in self.refseq_category.lower()

    @property
    def is_representative(self) -> bool:
        return "representative genome" in self.refseq_category.lower()

    def to_dict(self) -> dict:
        return {
            "accession": self.accession,
            "organism_name": self.organism_name,
            "tax_id": self.tax_id,
            "assembly_name": self.assembly_name,
            "assembly_level": self.assembly_level,
            "source_database": self.source_database,
            "refseq_category": self.refseq_category,
            "release_date": self.release_date,
            "submitter": self.submitter,
            "description": self.description,
            "synonym": self.synonym,
            "is_reference": self.is_reference,
            "is_representative": self.is_representative,
        }


@dataclass(frozen=True)
class NcbiSequenceRecord:
    assembly_accession: str
    chromosome: str
    sequence_name: str
    refseq_accession: str
    genbank_accession: str
    ucsc_style_name: str
    length: int
    role: str
    assembly_unit: str
    location_type: str

    @property
    def display_name(self) -> str:
        return self.chromosome or self.sequence_name or self.refseq_accession or self.genbank_accession

    def to_dict(self) -> dict:
        return {
            "assembly_accession": self.assembly_accession,
            "chromosome": self.chromosome,
            "sequence_name": self.sequence_name,
            "refseq_accession": self.refseq_accession,
            "genbank_accession": self.genbank_accession,
            "ucsc_style_name": self.ucsc_style_name,
            "length": self.length,
            "role": self.role,
            "assembly_unit": self.assembly_unit,
            "location_type": self.location_type,
            "display_name": self.display_name,
        }

@dataclass(frozen=True)
class GenomeDownload:
    source: str
    accession: str
    fasta_path: Path
    records: tuple[FastaRecordInfo, ...]
    requested_chromosomes: tuple[str, ...] = ()
    cached: bool = False

    def to_dict(self) -> dict:
        return {
            "source": self.source,
            "accession": self.accession,
            "fasta_path": str(self.fasta_path),
            "records": [record.to_dict() for record in self.records],
            "requested_chromosomes": list(self.requested_chromosomes),
            "cached": self.cached,
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



def _get_json(url: str, *, timeout: int = 60) -> dict:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "application/json",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except Exception as exc:
        raise RuntimeError(f"NCBI metadata request failed: {exc}") from exc


def _value(mapping: dict | None, *keys, default=None):
    if not mapping:
        return default
    for key in keys:
        if key in mapping and mapping[key] is not None:
            return mapping[key]
    return default


def _assembly_from_report(report: dict) -> NcbiAssembly:
    organism = _value(report, "organism", default={})
    info = _value(report, "assemblyInfo", "assembly_info", default={})
    accession = str(_value(report, "accession", default=""))
    return NcbiAssembly(
        accession=accession,
        organism_name=str(_value(organism, "organismName", "organism_name", "sciName", "sci_name", default="")),
        tax_id=_value(organism, "taxId", "tax_id"),
        assembly_name=str(_value(info, "assemblyName", "assembly_name", default="")),
        assembly_level=str(_value(info, "assemblyLevel", "assembly_level", default="")),
        source_database=str(_value(report, "sourceDatabase", "source_database", default="")),
        refseq_category=str(_value(info, "refseqCategory", "refseq_category", default="")),
        release_date=str(_value(info, "releaseDate", "release_date", default="")),
        submitter=str(_value(info, "submitter", default="")),
        description=str(_value(info, "description", default="")),
        synonym=str(_value(info, "synonym", default="")),
    )


def search_ncbi_assemblies(
    taxon: str,
    *,
    limit: int = 20,
    source: str = "refseq",
    exact_match: bool = True,
    reference_only: bool = False,
    timeout: int = 60,
) -> tuple[NcbiAssembly, ...]:
    """Find current NCBI genome assemblies by scientific/common taxon name or TaxID."""
    taxon = taxon.strip()
    if not taxon:
        raise ValueError("taxon is required")
    if limit < 1 or limit > 1000:
        raise ValueError("limit must be between 1 and 1000")
    if source not in {"all", "refseq", "genbank"}:
        raise ValueError("source must be all, refseq, or genbank")

    params: list[tuple[str, str]] = [
        ("page_size", str(limit)),
        ("tax_exact_match", str(exact_match).lower()),
        ("filters.exclude_paired_reports", "true"),
        ("filters.exclude_atypical", "true"),
        ("filters.assembly_version", "current"),
        ("filters.assembly_source", source),
    ]
    if reference_only:
        params.append(("filters.reference_only", "true"))
    query = urllib.parse.urlencode(params)
    url = (
        f"{NCBI_DATASETS_BASE}/genome/taxon/"
        f"{urllib.parse.quote(taxon, safe='')}/dataset_report?{query}"
    )
    payload = _get_json(url, timeout=timeout)
    assemblies = [
        _assembly_from_report(report)
        for report in payload.get("reports", [])
        if _value(report, "accession")
    ]

    def order(item: NcbiAssembly):
        category = 0 if item.is_reference else 1 if item.is_representative else 2
        refseq = 0 if item.accession.startswith("GCF_") else 1
        level = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}.get(
            item.assembly_level, 4
        )
        return (category, refseq, level, item.assembly_name, item.accession)

    return tuple(sorted(assemblies, key=order))


def list_ncbi_sequences(
    accession: str,
    *,
    chromosomes_only: bool = False,
    timeout: int = 60,
) -> tuple[NcbiSequenceRecord, ...]:
    """List sequence/chromosome metadata for one NCBI assembly accession."""
    accession = accession.strip().upper()
    if not (accession.startswith("GCF_") or accession.startswith("GCA_")):
        raise ValueError("NCBI assembly accession must start with GCF_ or GCA_")

    params = urllib.parse.urlencode({"page_size": 1000})
    url = (
        f"{NCBI_DATASETS_BASE}/genome/accession/"
        f"{urllib.parse.quote(accession, safe='')}/sequence_reports?{params}"
    )
    payload = _get_json(url, timeout=timeout)
    records: list[NcbiSequenceRecord] = []
    for report in payload.get("reports", []):
        chromosome = str(_value(report, "chrName", "chr_name", default=""))
        role = str(_value(report, "role", default=""))
        if chromosomes_only and not chromosome:
            continue
        records.append(
            NcbiSequenceRecord(
                assembly_accession=str(_value(report, "assemblyAccession", "assembly_accession", default=accession)),
                chromosome=chromosome,
                sequence_name=str(_value(report, "sequenceName", "sequence_name", default="")),
                refseq_accession=str(_value(report, "refseqAccession", "refseq_accession", default="")),
                genbank_accession=str(_value(report, "genbankAccession", "genbank_accession", default="")),
                ucsc_style_name=str(_value(report, "ucscStyleName", "ucsc_style_name", default="")),
                length=int(_value(report, "length", default=0) or 0),
                role=role,
                assembly_unit=str(_value(report, "assemblyUnit", "assembly_unit", default="")),
                location_type=str(
                    _value(
                        report,
                        "assignedMoleculeLocationType",
                        "assigned_molecule_location_type",
                        default="",
                    )
                ),
            )
        )

    return tuple(
        sorted(
            records,
            key=lambda item: (
                0 if item.chromosome else 1,
                int(item.chromosome) if item.chromosome.isdigit() else 10_000,
                item.chromosome,
                item.sequence_name,
            ),
        )
    )


def ncbi_cache_status(
    accession: str,
    *,
    chromosomes: Iterable[str] | None = None,
    cache_dir: str | Path | None = None,
) -> dict:
    accession = accession.strip().upper()
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
    return {
        "accession": accession,
        "chromosomes": list(chromosome_list),
        "cached": bool(existing),
        "fasta_path": str(existing[0]) if existing else None,
    }


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
            cached=True,
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
        cached=False,
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
