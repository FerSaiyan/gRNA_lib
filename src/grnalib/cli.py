from __future__ import annotations

import json
from pathlib import Path

import typer

from . import PrimeEdit, design_guides, design_prime_edit, rank_guides, resolve_nuclease

app = typer.Typer(help='Local-first CRISPR guide and prime-edit design toolkit.')
guide_app = typer.Typer(help='Guide-RNA design and ranking.')
prime_app = typer.Typer(help='Prime-edit candidate design.')
app.add_typer(guide_app, name='guide')
app.add_typer(prime_app, name='prime')


def _sequence(value: str | None, fasta: Path | None) -> str:
    if value:
        return value
    if not fasta:
        raise typer.BadParameter('Provide --sequence or --fasta')
    text = fasta.read_text()
    return ''.join(line.strip() for line in text.splitlines() if not line.startswith('>'))


def _dump(items) -> None:
    typer.echo(json.dumps(items, indent=2))


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
):
    spec = resolve_nuclease(nuclease, pam=pam)
    guides = design_guides(sequence, spec)
    genome_seq = _sequence(None, genome) if genome else None
    ranked = rank_guides(guides, spec=spec, genome_sequence=genome_seq, use_rs3=not no_rs3)
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
):
    spec = resolve_nuclease(nuclease, pam=pam)
    candidates = design_prime_edit(sequence, PrimeEdit(position, ref, alt), spec)
    _dump([x.to_dict() for x in candidates[:limit]])


def main() -> None:
    app()


if __name__ == '__main__':
    main()
