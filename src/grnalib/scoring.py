from __future__ import annotations

from collections import Counter
from copy import deepcopy
from typing import Iterable, Mapping

from .design import design_guides
from .models import GuideCandidate, NucleaseSpec, OffTargetSummary
from .sequence import hamming_distance, normalize_dna


def heuristic_on_target(guide: GuideCandidate) -> float:
    """Transparent fallback score in [0, 1]; not a learned efficacy model."""
    gc = guide.gc_fraction
    gc_score = max(0.0, 1.0 - abs(gc - 0.50) / 0.50)
    score = 0.35 + 0.55 * gc_score
    if 'TTTT' in guide.spacer:
        score -= 0.60
    if guide.spacer.startswith('G'):
        score += 0.05
    return round(max(0.0, min(1.0, score)), 4)


def rs3_sequence_score(guide: GuideCandidate, tracr: str = 'Hsu2013') -> tuple[float | None, str | None]:
    if guide.context is None or len(guide.context) != 30:
        return None, '30-nt SpCas9 context unavailable'
    try:
        from rs3.seq import predict_seq  # type: ignore
    except Exception:
        return None, 'optional rs3 package is not installed'
    try:
        value = float(predict_seq([guide.context], sequence_tracr=tracr)[0])
        return value, None
    except Exception as exc:
        return None, f'RS3 scoring failed: {exc}'


def local_offtarget_summary(
    guide: GuideCandidate,
    genome_sequence: str,
    spec: NucleaseSpec,
    *,
    max_mismatches: int = 4,
) -> OffTargetSummary:
    genome = normalize_dna(genome_sequence)
    sites = design_guides(genome, spec)
    counts: Counter[int] = Counter()
    for site in sites:
        if len(site.spacer) != len(guide.spacer):
            continue
        d = hamming_distance(guide.spacer, site.spacer)
        if d > max_mismatches:
            continue
        counts[d] += 1

    if counts[0] > 0:
        counts[0] -= 1

    burden = (
        counts[0] * 2.0 + counts[1] * 1.0 + counts[2] * 0.25 +
        counts[3] * 0.10 + counts[4] * 0.05
    )
    return OffTargetSummary(
        exact_other_sites=counts[0],
        mismatches_1=counts[1],
        mismatches_2=counts[2],
        mismatches_3=counts[3],
        mismatches_4=counts[4],
        specificity_proxy=round(1.0 / (1.0 + burden), 4),
    )


def rank_guides(
    guides: Iterable[GuideCandidate],
    *,
    spec: NucleaseSpec,
    genome_sequence: str | None = None,
    use_rs3: bool = True,
    tracr: str = 'Hsu2013',
    indexed_specificity: Mapping[str, float] | None = None,
    indexed_specificity_source: str | None = None,
) -> list[GuideCandidate]:
    ranked: list[GuideCandidate] = []
    for original in guides:
        guide = deepcopy(original)
        guide.scores['sequence_quality'] = heuristic_on_target(guide)
        guide.scores['hard_filter_pass'] = 'TTTT' not in guide.spacer

        rs3_score = None
        guide.scores['on_target_model'] = 'heuristic-v2'
        if use_rs3 and spec.name == 'SpCas9':
            rs3_score, reason = rs3_sequence_score(guide, tracr=tracr)
            guide.scores['rs3_sequence'] = rs3_score
            if rs3_score is not None:
                guide.scores['on_target_model'] = f'RS3-sequence/{tracr}'
            if reason:
                guide.scores['rs3_note'] = reason

        local_specificity = None
        if genome_sequence:
            local = local_offtarget_summary(guide, genome_sequence, spec)
            guide.scores['local_off_target'] = local.to_dict()
            local_specificity = local.specificity_proxy

        external_score = None
        if indexed_specificity is not None:
            value = indexed_specificity.get(guide.id)
            if value is not None:
                external_score = float(value)
                guide.scores['indexed_specificity'] = {
                    'score': external_score,
                    'source': indexed_specificity_source or 'external-index',
                }

        if external_score is not None:
            specificity = external_score
            guide.scores['specificity_used_for_ranking'] = 'indexed'
        elif indexed_specificity is not None:
            # An indexed run was requested but this guide received no result.
            # Never interpret missing indexed evidence as perfect specificity.
            specificity = local_specificity if local_specificity is not None else -1.0
            guide.scores['indexed_specificity'] = {
                'score': None,
                'source': indexed_specificity_source or 'external-index',
                'note': 'indexed backend returned no specificity score for this guide',
            }
            guide.scores['specificity_used_for_ranking'] = 'indexed-missing'
        elif local_specificity is not None:
            specificity = local_specificity
            guide.scores['specificity_used_for_ranking'] = 'local-proxy'
        else:
            specificity = 1.0
            guide.scores['specificity_used_for_ranking'] = 'not-available'

        guide.scores['_sort'] = (
            1 if guide.scores['hard_filter_pass'] else 0,
            specificity,
            rs3_score if rs3_score is not None else guide.scores['sequence_quality'],
        )
        ranked.append(guide)

    ranked.sort(key=lambda g: g.scores['_sort'], reverse=True)
    for i, guide in enumerate(ranked, 1):
        guide.rank = i
        guide.scores.pop('_sort', None)
    return ranked
