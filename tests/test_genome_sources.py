from __future__ import annotations

import io
import zipfile

import pytest

from grnalib import genome_sources


def test_fasta_record_index_and_selection(tmp_path):
    path = tmp_path / "multi.fa"
    path.write_text(
        ">chr1 first chromosome\nACGTNN\nACGT\n>chr2 second chromosome\nTTTT\n",
        encoding="utf-8",
    )

    records = genome_sources.fasta_records(path)
    assert [(record.name, record.length) for record in records] == [
        ("chr1", 10),
        ("chr2", 4),
    ]

    with pytest.raises(ValueError, match="contains 2 FASTA records"):
        genome_sources.read_fasta_record(path)

    name, sequence = genome_sources.read_fasta_record(path, "chr2")
    assert name == "chr2"
    assert sequence == "TTTT"


def test_parse_fasta_text_normalizes_case_and_whitespace():
    records = genome_sources.parse_fasta_text(">target description\nacgt n\n")
    assert records == [("target", "ACGTN")]


def test_ncbi_download_builds_chromosome_filtered_request_and_caches(tmp_path, monkeypatch):
    seen = {}

    def fake_download(url, destination, *, timeout):
        seen["url"] = url
        with zipfile.ZipFile(destination, "w") as archive:
            archive.writestr(
                "ncbi_dataset/data/GCF_000001405.40/test_genomic.fna",
                ">NC_000023.11 chromosome X\nACGTACGT\n",
            )

    monkeypatch.setattr(genome_sources, "_download_to_file", fake_download)

    result = genome_sources.download_ncbi_genome(
        "GCF_000001405.40",
        chromosomes=["X"],
        cache_dir=tmp_path,
    )
    assert "include_annotation_type=GENOME_FASTA" in seen["url"]
    assert "chromosomes=X" in seen["url"]
    assert result.records[0].name == "NC_000023.11"
    assert result.records[0].length == 8
    assert result.fasta_path.exists()
    assert result.cached is False

    def should_not_download(*args, **kwargs):
        raise AssertionError("cache should have been used")

    monkeypatch.setattr(genome_sources, "_download_to_file", should_not_download)
    cached = genome_sources.download_ncbi_genome(
        "GCF_000001405.40",
        chromosomes=["X"],
        cache_dir=tmp_path,
    )
    assert cached.fasta_path == result.fasta_path
    assert cached.cached is True


def test_ensembl_region_request_returns_normalized_sequence(monkeypatch):
    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self):
            return b"acgtacgt\n"

    seen = {}

    def fake_urlopen(request, timeout):
        seen["url"] = request.full_url
        return Response()

    monkeypatch.setattr(genome_sources.urllib.request, "urlopen", fake_urlopen)
    sequence = genome_sources.fetch_ensembl_region(
        "homo_sapiens",
        "17:7668402..7668500:1",
    )
    assert sequence == "ACGTACGT"
    assert "/sequence/region/homo_sapiens/17:7668402..7668500:1" in seen["url"]


def test_ncbi_assembly_search_parses_and_prioritizes_reference(monkeypatch):
    payload = {
        "reports": [
            {
                "accession": "GCF_999999999.1",
                "sourceDatabase": "SOURCE_DATABASE_REFSEQ",
                "organism": {"taxId": 9606, "organismName": "Homo sapiens"},
                "assemblyInfo": {
                    "assemblyName": "Other",
                    "assemblyLevel": "Chromosome",
                    "refseqCategory": "representative genome",
                    "releaseDate": "2026-01-01",
                },
            },
            {
                "accession": "GCF_000001405.40",
                "sourceDatabase": "SOURCE_DATABASE_REFSEQ",
                "organism": {"taxId": 9606, "organismName": "Homo sapiens"},
                "assemblyInfo": {
                    "assemblyName": "GRCh38.p14",
                    "assemblyLevel": "Chromosome",
                    "refseqCategory": "reference genome",
                    "releaseDate": "2022-02-03",
                    "synonym": "hg38",
                },
            },
        ]
    }
    seen = {}

    def fake_get_json(url, *, timeout):
        seen["url"] = url
        return payload

    monkeypatch.setattr(genome_sources, "_get_json", fake_get_json)
    assemblies = genome_sources.search_ncbi_assemblies("Homo sapiens", limit=10)
    assert assemblies[0].accession == "GCF_000001405.40"
    assert assemblies[0].is_reference is True
    assert assemblies[0].synonym == "hg38"
    assert "tax_exact_match=true" in seen["url"]
    assert "filters.assembly_source=refseq" in seen["url"]


def test_ncbi_sequence_report_lists_chromosome_metadata(monkeypatch):
    payload = {
        "reports": [
            {
                "assemblyAccession": "GCF_000001405.40",
                "chrName": "17",
                "sequenceName": "17",
                "refseqAccession": "NC_000017.11",
                "genbankAccession": "CM000679.2",
                "ucscStyleName": "chr17",
                "length": 83257441,
                "role": "assembled-molecule",
                "assemblyUnit": "Primary Assembly",
                "assignedMoleculeLocationType": "Chromosome",
            },
            {
                "assemblyAccession": "GCF_000001405.40",
                "sequenceName": "KI270728.1",
                "refseqAccession": "NT_187361.1",
                "length": 1872759,
                "role": "unplaced-scaffold",
            },
        ]
    }

    monkeypatch.setattr(genome_sources, "_get_json", lambda url, timeout: payload)
    records = genome_sources.list_ncbi_sequences("GCF_000001405.40")
    assert records[0].chromosome == "17"
    assert records[0].ucsc_style_name == "chr17"
    assert records[0].refseq_accession == "NC_000017.11"

    chromosomes = genome_sources.list_ncbi_sequences(
        "GCF_000001405.40", chromosomes_only=True
    )
    assert [record.chromosome for record in chromosomes] == ["17"]


def test_ncbi_cache_status_reports_existing_download(tmp_path):
    target = tmp_path / "ncbi" / "GCF_000001405.40__17"
    target.mkdir(parents=True)
    fasta = target / "human_genomic.fna"
    fasta.write_text(">chr17\nACGT\n", encoding="utf-8")
    status = genome_sources.ncbi_cache_status(
        "GCF_000001405.40",
        chromosomes=["17"],
        cache_dir=tmp_path,
    )
    assert status["cached"] is True
    assert status["fasta_path"] == str(fasta)
