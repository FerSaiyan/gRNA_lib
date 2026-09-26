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

    def should_not_download(*args, **kwargs):
        raise AssertionError("cache should have been used")

    monkeypatch.setattr(genome_sources, "_download_to_file", should_not_download)
    cached = genome_sources.download_ncbi_genome(
        "GCF_000001405.40",
        chromosomes=["X"],
        cache_dir=tmp_path,
    )
    assert cached.fasta_path == result.fasta_path


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
