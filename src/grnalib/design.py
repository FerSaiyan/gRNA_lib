from __future__ import annotations

from .models import GuideCandidate, NucleaseSpec
from .sequence import iupac_match, normalize_dna, reverse_complement


def _u6_spacer(spacer: str) -> str:
    # U6 promoters generally prefer a 5' G. Keep the genomic spacer immutable.
    return spacer if spacer.startswith('G') else 'G' + spacer


def _warnings(spacer: str, expressed_spacer: str) -> list[str]:
    warnings: list[str] = []
    if 'TTTT' in spacer:
        warnings.append('poly-T termination risk (TTTT)')
    gc = (spacer.count('G') + spacer.count('C')) / len(spacer)
    if gc < 0.30:
        warnings.append('low GC content (<30%)')
    elif gc >= 0.80:
        warnings.append('high GC content (>=80%)')
    if expressed_spacer != spacer:
        warnings.append("5' G added for U6 expression; genomic spacer unchanged")
    return warnings


def _oriented_candidates(oriented: str, spec: NucleaseSpec) -> list[dict]:
    out: list[dict] = []
    pam_len = len(spec.pam)
    L = spec.spacer_length
    for pam_start in range(0, len(oriented) - pam_len + 1):
        pam_seq = oriented[pam_start:pam_start + pam_len]
        if not iupac_match(pam_seq, spec.pam):
            continue
        if spec.pam_side == '3prime':
            start, end = pam_start - L, pam_start
            if start < 0:
                continue
            cut = end - spec.cut_offset
            context_start = start - 4
            context_end = pam_start + pam_len + 3
        else:
            start, end = pam_start + pam_len, pam_start + pam_len + L
            if end > len(oriented):
                continue
            cut = start + spec.cut_offset
            context_start = pam_start - 3
            context_end = end + 4
        spacer = oriented[start:end]
        context = oriented[context_start:context_end] if context_start >= 0 and context_end <= len(oriented) else None
        out.append({
            'spacer': spacer,
            'pam': pam_seq,
            'start': start,
            'end': end,
            'pam_start': pam_start,
            'pam_end': pam_start + pam_len,
            'cut': cut,
            'context': context,
        })
    return out


def design_guides(sequence: str, spec: NucleaseSpec) -> list[GuideCandidate]:
    seq = normalize_dna(sequence)
    rc = reverse_complement(seq)
    n = len(seq)
    results: list[GuideCandidate] = []

    for raw in _oriented_candidates(seq, spec):
        spacer = raw['spacer']
        expressed = _u6_spacer(spacer)
        results.append(GuideCandidate(
            spacer=spacer,
            expressed_spacer=expressed,
            pam=raw['pam'],
            strand='+',
            start=raw['start'],
            end=raw['end'],
            pam_start=raw['pam_start'],
            pam_end=raw['pam_end'],
            cut_site=raw['cut'],
            context=raw['context'],
            gc_fraction=(spacer.count('G') + spacer.count('C')) / len(spacer),
            warnings=_warnings(spacer, expressed),
        ))

    for raw in _oriented_candidates(rc, spec):
        spacer = raw['spacer']
        expressed = _u6_spacer(spacer)
        # Map reverse-complement half-open coordinates back to the reference.
        start, end = n - raw['end'], n - raw['start']
        pam_start, pam_end = n - raw['pam_end'], n - raw['pam_start']
        cut = n - raw['cut']
        results.append(GuideCandidate(
            spacer=spacer,
            expressed_spacer=expressed,
            pam=raw['pam'],
            strand='-',
            start=start,
            end=end,
            pam_start=pam_start,
            pam_end=pam_end,
            cut_site=cut,
            context=raw['context'],
            gc_fraction=(spacer.count('G') + spacer.count('C')) / len(spacer),
            warnings=_warnings(spacer, expressed),
        ))

    results.sort(key=lambda g: (g.start, g.strand, g.pam_start))
    return results
