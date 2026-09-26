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
from .backends import backend_status, build_crisprware_index, score_with_crisprware

app = typer.Typer(help='Local-first CRISPR guide and prime-edit design toolkit.')
guide_app = typer.Typer(help='Guide-RNA design and ranking.')
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
    crisprware_executable: str = 'crisprware',
):
    output_directory.mkdir(parents=True, exist_ok=True)
    build_crisprware_index(
        fasta,
        pam=pam,
        spacer_length=spacer_length,
        output_directory=output_directory,
        executable=crisprware_executable,
    )
    _dump({'ok': True, 'backend': 'crisprware', 'output_directory': str(output_directory)})


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
