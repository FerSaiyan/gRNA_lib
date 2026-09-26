from __future__ import annotations

import json
from pathlib import Path

import typer

from . import (
    PrimeEdit,
    design_guides,
    design_prime_edit,
    rank_guides,
    rank_prime_candidates,
    resolve_nuclease,
    score_prime_candidates_deepprime,
)
from .backends import (
    backend_status,
    build_crisprware_index,
    profile_with_crispr_ots,
    score_with_crisprware,
)
from .genome_sources import (
    download_ncbi_genome,
    fetch_ensembl_region,
    list_ncbi_sequences,
    search_ncbi_assemblies,
)
from .reference_workflows import build_or_reuse_ncbi_index

app = typer.Typer(help='Local-first CRISPR guide and prime-edit design toolkit.')
guide_app = typer.Typer(help='Guide-RNA design, ranking and off-target profiling.')
prime_app = typer.Typer(help='Prime-edit candidate design.')
genome_app = typer.Typer(help='Genome/off-target index operations.')
backend_app = typer.Typer(help='Optional scoring backend information.')
app.add_typer(guide_app, name='guide')
app.add_typer(prime_app, name='prime')
app.add_typer(genome_app, name='genome')
app.add_typer(backend_app, name='backend')


def _sequence(value: str | None, fasta: Path | None) -> str:
    if value:
        return value
    if not fasta:
        raise typer.BadParameter('Provide --sequence or --fasta')
    text = fasta.read_text()
    return ''.join(line.strip() for line in text.splitlines() if not line.startswith('>'))


def _dump(items) -> None:
    typer.echo(json.dumps(items, indent=2))


@backend_app.command('status')
def show_backend_status():
    _dump(backend_status())


@genome_app.command('index')
def genome_index(
    fasta: Path = typer.Option(..., exists=True, dir_okay=False),
    pam: str = 'NGG',
    spacer_length: int = 20,
    output_directory: Path = Path('.'),
    pam_side: str = typer.Option('3prime', help='3prime or 5prime'),
    bin_width: int | None = typer.Option(None, help='Optional crispr-ots bin width (1-15)'),
    crisprware_executable: str = 'crisprware',
):
    output_directory.mkdir(parents=True, exist_ok=True)
    build_crisprware_index(
        fasta,
        pam=pam,
        spacer_length=spacer_length,
        pam_5_prime=pam_side == '5prime',
        bin_width=bin_width,
        output_directory=output_directory,
        executable=crisprware_executable,
    )
    _dump({'ok': True, 'backend': 'crisprware', 'output_directory': str(output_directory)})


@genome_app.command('search')
def genome_search(
    taxon: str = typer.Option(..., help='Scientific/common name or NCBI TaxID'),
    limit: int = typer.Option(20, min=1, max=1000),
    source: str = typer.Option('refseq', help='refseq, genbank, or all'),
    exact_match: bool = typer.Option(True, help='Require an exact taxon match'),
    reference_only: bool = typer.Option(False, help='Return only reference assemblies'),
):
    """Find current NCBI assemblies without requiring an accession."""
    assemblies = search_ncbi_assemblies(
        taxon,
        limit=limit,
        source=source,
        exact_match=exact_match,
        reference_only=reference_only,
    )
    _dump([item.to_dict() for item in assemblies])


@genome_app.command('sequences')
def genome_sequences(
    accession: str = typer.Option(..., help='NCBI assembly accession'),
    chromosomes_only: bool = typer.Option(True, help='Hide unplaced/unlocalized sequences'),
):
    """List chromosome/sequence metadata for an NCBI assembly."""
    records = list_ncbi_sequences(accession, chromosomes_only=chromosomes_only)
    _dump([item.to_dict() for item in records])


@genome_app.command('index-ncbi')
def genome_index_ncbi(
    accession: str = typer.Option(..., help='NCBI assembly accession'),
    chromosome: list[str] | None = typer.Option(None, '--chromosome', help='Repeat for multiple chromosomes'),
    pam: str = 'NGG',
    spacer_length: int = 20,
    pam_side: str = typer.Option('3prime', help='3prime or 5prime'),
    bin_width: int | None = typer.Option(None, help='Optional crispr-ots bin width (1-15)'),
    cache_dir: Path | None = typer.Option(None),
    crisprware_executable: str = 'crisprware',
    force: bool = typer.Option(False, help='Rebuild even if a matching cached index exists'),
):
    """Download/cache an NCBI reference and build or reuse a matching off-target index."""
    result = build_or_reuse_ncbi_index(
        accession,
        chromosomes=chromosome,
        pam=pam,
        spacer_length=spacer_length,
        pam_side=pam_side,
        bin_width=bin_width,
        cache_dir=cache_dir,
        crisprware_executable=crisprware_executable,
        force=force,
    )
    _dump(result.to_dict())


@genome_app.command('fetch-ncbi')
def genome_fetch_ncbi(
    accession: str = typer.Option(..., help='NCBI assembly accession (GCF_ or GCA_)'),
    chromosome: list[str] | None = typer.Option(None, '--chromosome', help='Chromosome label; repeat for multiple'),
    cache_dir: Path | None = typer.Option(None, help='Optional cache directory'),
    force: bool = typer.Option(False, help='Re-download even if cached'),
):
    """Download and cache an assembly or selected chromosomes from NCBI Datasets."""
    result = download_ncbi_genome(
        accession,
        chromosomes=chromosome,
        cache_dir=cache_dir,
        force=force,
    )
    _dump(result.to_dict())


@genome_app.command('fetch-ensembl')
def genome_fetch_ensembl(
    species: str = typer.Option(..., help='Ensembl species name or alias'),
    region: str = typer.Option(..., help='Region such as 17:7668402..7687550:1'),
):
    """Fetch a genomic region from Ensembl REST and print the sequence as JSON."""
    sequence = fetch_ensembl_region(species, region)
    _dump({'source': 'Ensembl REST', 'species': species, 'region': region, 'sequence': sequence})


@guide_app.command('design')
def guide_design(
    sequence: str | None = typer.Option(None),
    fasta: Path | None = typer.Option(None, exists=True),
    nuclease: str = 'SpCas9',
    pam: str | None = None,
    spacer_length: int | None = None,
    pam_side: str | None = None,
):
    seq = _sequence(sequence, fasta)
    spec = resolve_nuclease(nuclease, pam=pam, spacer_length=spacer_length, pam_side=pam_side)
    _dump([g.to_dict() for g in design_guides(seq, spec)])


@guide_app.command('rank')
def guide_rank(
    sequence: str = typer.Option(..., help='Target sequence'),
    genome: Path | None = typer.Option(None, exists=True, help='Optional FASTA/text sequence for local off-target scan'),
    nuclease: str = 'SpCas9',
    pam: str | None = None,
    no_rs3: bool = False,
    crisprware_index: Path | None = typer.Option(None, help='Optional installed crispr-ots/GuideScan2 index'),
    chromosome: str | None = typer.Option(None, help='Chromosome/contig for the target sequence'),
    reference_start: int = typer.Option(0, help='0-based genomic coordinate of sequence base 0'),
    crisprware_executable: str = 'crisprware',
    threads: int = 4,
    mismatches: int = 3,
    rna_bulges: int = 0,
    dna_bulges: int = 0,
):
    spec = resolve_nuclease(nuclease, pam=pam)
    guides = design_guides(sequence, spec)
    genome_seq = _sequence(None, genome) if genome else None

    indexed = None
    source = None
    if crisprware_index is not None:
        if not chromosome:
            raise typer.BadParameter('--chromosome is required with --crisprware-index')
        result = score_with_crisprware(
            guides,
            spec,
            index=crisprware_index,
            chromosome=chromosome,
            reference_start=reference_start,
            executable=crisprware_executable,
            threads=threads,
            mismatches=mismatches,
            rna_bulges=rna_bulges,
            dna_bulges=dna_bulges,
        )
        indexed = result.scores
        source = f'{result.backend}:{result.index}:{result.specificity_column}'

    ranked = rank_guides(
        guides,
        spec=spec,
        genome_sequence=genome_seq,
        use_rs3=not no_rs3,
        indexed_specificity=indexed,
        indexed_specificity_source=source,
    )
    _dump([g.to_dict() for g in ranked])


@guide_app.command('profile')
def guide_profile(
    sequence: str = typer.Option(..., help='Target sequence containing candidate guides'),
    crispr_ots_index: Path = typer.Option(..., help='crispr-ots index prefix or .crot file'),
    chromosome: str = typer.Option(..., help='Chromosome/contig for the target sequence'),
    reference_start: int = typer.Option(0, help='0-based genomic coordinate of sequence base 0'),
    nuclease: str = 'SpCas9',
    pam: str | None = None,
    threads: int = 4,
    mismatches: int = 4,
    hit_limit: int = typer.Option(100, help='Representative loci retained per guide; counts are not truncated'),
    crispr_ots_executable: str = 'crispr-ots',
):
    """Profile whole-genome off-target mismatch classes and representative loci."""
    spec = resolve_nuclease(nuclease, pam=pam)
    guides = design_guides(sequence, spec)
    result = profile_with_crispr_ots(
        guides,
        spec,
        index=crispr_ots_index,
        chromosome=chromosome,
        reference_start=reference_start,
        executable=crispr_ots_executable,
        threads=threads,
        mismatches=mismatches,
        hit_limit=hit_limit,
    )
    _dump(result.to_dict())


@prime_app.command('design')
def prime_design(
    sequence: str = typer.Option(..., help='Reference sequence'),
    position: int = typer.Option(..., help='0-based edit start in reference'),
    ref: str = typer.Option('', help='Reference allele; empty for insertion'),
    alt: str = typer.Option('', help='Desired allele; empty for deletion'),
    nuclease: str = 'SpCas9',
    pam: str | None = None,
    limit: int = 100,
    scorer: str = typer.Option('structural', help='structural or deepprime'),
    pe_system: str = 'PE2max',
    cell_type: str = 'HEK293T',
):
    spec = resolve_nuclease(nuclease, pam=pam)
    edit = PrimeEdit(position, ref, alt)
    candidates = design_prime_edit(sequence, edit, spec)
    if scorer.lower() == 'deepprime':
        candidates = score_prime_candidates_deepprime(
            sequence,
            edit,
            candidates,
            pe_system=pe_system,
            cell_type=cell_type,
            pam=spec.pam,
        )
    elif scorer.lower() != 'structural':
        raise typer.BadParameter('--scorer must be structural or deepprime')
    candidates = rank_prime_candidates(candidates)
    _dump([x.to_dict() for x in candidates[:limit]])


def main() -> None:
    app()


if __name__ == '__main__':
    main()
