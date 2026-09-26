from __future__ import annotations

from pathlib import Path

from grnalib.genome_sources import FastaRecordInfo, GenomeDownload
from grnalib import reference_workflows


def test_ncbi_index_paths_include_reference_and_enzyme(tmp_path):
    output_dir, prefix = reference_workflows.ncbi_index_paths(
        "GCF_000001405.40",
        chromosomes=["17"],
        pam="NGG",
        spacer_length=20,
        pam_side="3prime",
        cache_dir=tmp_path,
        fasta_name="human_genomic.fna",
    )
    assert "GCF_000001405.40__17__NGG__20__3prime" in str(output_dir)
    assert prefix == output_dir / "human_genomic_crisprots" / "human_genomic_crisprots"


def test_build_or_reuse_ncbi_index_creates_manifest_and_reuses(tmp_path, monkeypatch):
    fasta = tmp_path / "human_genomic.fna"
    fasta.write_text(">chr17\nACGTACGT\n", encoding="utf-8")
    download = GenomeDownload(
        source="NCBI Datasets",
        accession="GCF_000001405.40",
        fasta_path=fasta,
        records=(FastaRecordInfo("chr17", "chr17", 8),),
        requested_chromosomes=("17",),
        cached=True,
    )
    monkeypatch.setattr(
        reference_workflows,
        "download_ncbi_genome",
        lambda *args, **kwargs: download,
    )

    calls = []

    def fake_build(path, **kwargs):
        calls.append((path, kwargs))
        outdir = Path(kwargs["output_directory"])
        prefix = outdir / "human_genomic_crisprots" / "human_genomic_crisprots"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        Path(str(prefix) + ".crot").write_bytes(b"index")

    monkeypatch.setattr(reference_workflows, "build_crisprware_index", fake_build)

    result = reference_workflows.build_or_reuse_ncbi_index(
        "GCF_000001405.40",
        chromosomes=["17"],
        pam="TTTV",
        spacer_length=23,
        pam_side="5prime",
        cache_dir=tmp_path,
    )
    assert result.cached is False
    assert Path(str(result.index_prefix) + ".crot").exists()
    assert calls[0][1]["pam_5_prime"] is True
    assert (result.index_prefix.parents[1] / "grnalib-index.json").exists()

    calls.clear()
    reused = reference_workflows.build_or_reuse_ncbi_index(
        "GCF_000001405.40",
        chromosomes=["17"],
        pam="TTTV",
        spacer_length=23,
        pam_side="5prime",
        cache_dir=tmp_path,
    )
    assert reused.cached is True
    assert calls == []
