from __future__ import annotations

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
from .genome_sources import download_ncbi_genome, fetch_ensembl_region

try:
    from mcp.server.fastmcp import FastMCP
except ImportError:
    try:
        from fastmcp import FastMCP  # type: ignore
    except ImportError as exc:
        raise RuntimeError('Install gRNA-lib with the mcp extra: pip install .[mcp]') from exc

mcp = FastMCP('gRNA Library')


@mcp.tool()
def get_backend_status() -> dict:
    """Report whether optional RS3, DeepPrime/GenET, CRISPRware and crispr-ots backends are installed."""
    return backend_status()


@mcp.tool()
def design_grnas(sequence: str, nuclease: str = 'SpCas9', pam: str | None = None) -> list[dict]:
    """Enumerate guide candidates on both strands with normalized reference coordinates."""
    spec = resolve_nuclease(nuclease, pam=pam)
    return [g.to_dict() for g in design_guides(sequence, spec)]


@mcp.tool()
def rank_grnas(
    sequence: str,
    genome_sequence: str | None = None,
    nuclease: str = 'SpCas9',
    pam: str | None = None,
    crisprware_index: str | None = None,
    chromosome: str | None = None,
    reference_start: int = 0,
    threads: int = 4,
    mismatches: int = 3,
    rna_bulges: int = 0,
    dna_bulges: int = 0,
) -> list[dict]:
    """Rank guides with optional indexed whole-genome specificity."""
    spec = resolve_nuclease(nuclease, pam=pam)
    guides = design_guides(sequence, spec)
    indexed = None
    source = None
    if crisprware_index:
        if not chromosome:
            raise ValueError('chromosome is required when crisprware_index is supplied')
        result = score_with_crisprware(
            guides,
            spec,
            index=crisprware_index,
            chromosome=chromosome,
            reference_start=reference_start,
            threads=threads,
            mismatches=mismatches,
            rna_bulges=rna_bulges,
            dna_bulges=dna_bulges,
        )
        indexed = result.scores
        source = f'{result.backend}:{result.index}:{result.specificity_column}'
    return [
        g.to_dict()
        for g in rank_guides(
            guides,
            spec=spec,
            genome_sequence=genome_sequence,
            indexed_specificity=indexed,
            indexed_specificity_source=source,
        )
    ]


@mcp.tool()
def profile_genome_offtargets(
    sequence: str,
    crispr_ots_index: str,
    chromosome: str,
    reference_start: int = 0,
    nuclease: str = 'SpCas9',
    pam: str | None = None,
    mismatches: int = 4,
    threads: int = 4,
    hit_limit: int = 100,
    crispr_ots_executable: str = 'crispr-ots',
) -> dict:
    """Enumerate whole-genome off-target mismatch classes and representative loci."""
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
    return result.to_dict()


@mcp.tool()
def fetch_ncbi_genome(
    accession: str,
    chromosomes: list[str] | None = None,
) -> dict:
    """Download/cache an NCBI assembly or selected chromosomes and return local FASTA metadata."""
    result = download_ncbi_genome(accession, chromosomes=chromosomes)
    return result.to_dict()


@mcp.tool()
def fetch_ensembl_region_sequence(species: str, region: str) -> dict:
    """Fetch a genomic region from Ensembl REST."""
    sequence = fetch_ensembl_region(species, region)
    return {
        'source': 'Ensembl REST',
        'species': species,
        'region': region,
        'length': len(sequence),
        'sequence': sequence,
    }


@mcp.tool()
def build_offtarget_index(
    fasta: str,
    pam: str = 'NGG',
    spacer_length: int = 20,
    output_directory: str = '.',
) -> dict:
    """Build a CRISPRware crispr-ots off-target index from a local FASTA."""
    build_crisprware_index(
        fasta,
        pam=pam,
        spacer_length=spacer_length,
        output_directory=output_directory,
    )
    return {'ok': True, 'backend': 'crisprware', 'output_directory': output_directory}


@mcp.tool()
def design_prime_candidates(
    sequence: str,
    position: int,
    ref: str,
    alt: str,
    nuclease: str = 'SpCas9',
    pam: str | None = None,
    limit: int = 100,
    scorer: str = 'structural',
    pe_system: str = 'PE2max',
    cell_type: str = 'HEK293T',
) -> list[dict]:
    """Enumerate and rank pegRNAs, optionally using learned DeepPrime efficiency."""
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
        raise ValueError('scorer must be structural or deepprime')
    candidates = rank_prime_candidates(candidates)
    return [c.to_dict() for c in candidates[:limit]]


@mcp.tool()
def explain_candidate(candidate: dict) -> dict:
    """Return auditable ranking fields from a previously returned candidate."""
    return {
        'id': candidate.get('id'),
        'rank': candidate.get('rank'),
        'scores': candidate.get('scores', {}),
        'warnings': candidate.get('warnings', []),
        'interpretation': (
            'Guide ranking is hierarchical: hard filters, indexed specificity when available, '
            'then learned on-target activity. Prime ranking prefers a learned DeepPrime score '
            'when present and otherwise uses the structural fallback.'
        ),
    }


def main() -> None:
    mcp.run()


if __name__ == '__main__':
    main()
