from __future__ import annotations

from . import PrimeEdit, design_guides, design_prime_edit, rank_guides, resolve_nuclease

try:
    from mcp.server.fastmcp import FastMCP
except ImportError:  # pragma: no cover
    try:
        from fastmcp import FastMCP  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError('Install gRNA-lib with the mcp extra: pip install .[mcp]') from exc

mcp = FastMCP('gRNA Library')


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
) -> list[dict]:
    """Rank guides hierarchically and return component scores plus warnings."""
    spec = resolve_nuclease(nuclease, pam=pam)
    guides = design_guides(sequence, spec)
    return [g.to_dict() for g in rank_guides(guides, spec=spec, genome_sequence=genome_sequence)]


@mcp.tool()
def design_prime_candidates(
    sequence: str,
    position: int,
    ref: str,
    alt: str,
    nuclease: str = 'SpCas9',
    pam: str | None = None,
    limit: int = 100,
) -> list[dict]:
    """Enumerate pegRNA PBS/RTT combinations and PE3/PE3b nicking-guide candidates."""
    spec = resolve_nuclease(nuclease, pam=pam)
    candidates = design_prime_edit(sequence, PrimeEdit(position, ref, alt), spec)
    return [c.to_dict() for c in candidates[:limit]]


@mcp.tool()
def explain_candidate(candidate: dict) -> dict:
    """Return the auditable ranking fields from a previously returned candidate."""
    return {
        'id': candidate.get('id'),
        'rank': candidate.get('rank'),
        'scores': candidate.get('scores', {}),
        'warnings': candidate.get('warnings', []),
        'interpretation': (
            'Ranking is hierarchical: hard filters, specificity when available, then learned RS3 or the transparent fallback.'
        ),
    }


def main() -> None:
    mcp.run()


if __name__ == '__main__':
    main()
