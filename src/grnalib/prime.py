from __future__ import annotations

from copy import deepcopy
from typing import Iterable

from .design import design_guides
from .models import GuideCandidate, NickingGuide, NucleaseSpec, PrimeCandidate, PrimeEdit
from .scoring import heuristic_on_target
from .sequence import normalize_dna, reverse_complement


def apply_edit(reference: str, edit: PrimeEdit) -> str:
    seq = normalize_dna(reference)
    ref = normalize_dna(edit.ref) if edit.ref else ''
    alt = normalize_dna(edit.alt) if edit.alt else ''
    if edit.position < 0 or edit.position + len(ref) > len(seq):
        raise ValueError('Edit is outside the reference sequence')
    observed = seq[edit.position:edit.position + len(ref)]
    if observed != ref:
        raise ValueError(f"Edit REF mismatch: expected {ref!r} at {edit.position}, found {observed!r}")
    return seq[:edit.position] + alt + seq[edit.position + len(ref):]


def _oriented_edit_start(reference_len: int, edit: PrimeEdit, strand: str) -> int:
    if strand == '+':
        return edit.position
    return reference_len - (edit.position + len(edit.ref))


def _oriented_sequences(reference: str, edited: str, strand: str) -> tuple[str, str]:
    if strand == '+':
        return reference, edited
    return reverse_complement(reference), reverse_complement(edited)


def _oriented_nick(reference_len: int, guide: GuideCandidate) -> int:
    assert guide.cut_site is not None
    return guide.cut_site if guide.strand == '+' else reference_len - guide.cut_site


def _pam_disrupted(reference: str, edited: str, guide: GuideCandidate) -> bool:
    # Reference coordinates only remain directly comparable after equal-length edits.
    # Indel-aware PAM remapping is intentionally deferred rather than guessed.
    if len(reference) != len(edited):
        return False
    before = reference[guide.pam_start:guide.pam_end]
    after = edited[guide.pam_start:guide.pam_end]
    if guide.strand == '-':
        before, after = reverse_complement(before), reverse_complement(after)
    return before != after


def _classify_nick_guide(reference: str, edited: str, guide: GuideCandidate) -> str:
    # PE3b if the edited allele changes the candidate protospacer/PAM sequence at this locus.
    if len(reference) != len(edited):
        return 'PE3'
    if guide.strand == '+':
        before = reference[guide.start:guide.pam_end]
        after = edited[guide.start:guide.pam_end]
    else:
        lo, hi = min(guide.pam_start, guide.start), max(guide.pam_end, guide.end)
        before = reference[lo:hi]
        after = edited[lo:hi]
    return 'PE3b' if before != after else 'PE3'


def design_prime_edit(
    reference: str,
    edit: PrimeEdit,
    spec: NucleaseSpec,
    *,
    pbs_lengths: Iterable[int] = range(10, 16),
    rtt_lengths: Iterable[int] = range(10, 31),
    max_nick_to_edit: int = 40,
    pe3_min_distance: int = 40,
    pe3_max_distance: int = 150,
) -> list[PrimeCandidate]:
    if spec.name not in {'SpCas9', 'SpCas9-NG', 'Custom'} or spec.pam_side != '3prime':
        raise ValueError('Prime-edit extension generation currently supports 3-prime-PAM Cas9-family specs')

    ref = normalize_dna(reference)
    edited = apply_edit(ref, edit)
    guides = design_guides(ref, spec)
    edited_guides = design_guides(edited, spec)
    out: list[PrimeCandidate] = []
    n = len(ref)

    for guide in guides:
        if guide.cut_site is None:
            continue
        oriented_ref, oriented_edited = _oriented_sequences(ref, edited, guide.strand)
        nick = _oriented_nick(n, guide)
        edit_start = _oriented_edit_start(n, edit, guide.strand)
        if edit_start < nick:
            continue
        nick_to_edit = edit_start - nick
        if nick_to_edit > max_nick_to_edit:
            continue

        for pbs_len in pbs_lengths:
            if pbs_len <= 0 or nick - pbs_len < 0:
                continue
            pbs = reverse_complement(oriented_ref[nick - pbs_len:nick])
            for rtt_len in rtt_lengths:
                if rtt_len <= nick_to_edit:
                    continue
                if nick + rtt_len > len(oriented_edited):
                    continue
                rtt = reverse_complement(oriented_edited[nick:nick + rtt_len])
                extension = rtt + pbs
                candidate = PrimeCandidate(
                    spacer=deepcopy(guide),
                    pbs_sequence=pbs,
                    pbs_length=pbs_len,
                    rtt_sequence=rtt,
                    rtt_length=rtt_len,
                    extension_sequence=extension,
                    nick_to_edit=nick_to_edit,
                    pam_disrupted=_pam_disrupted(ref, edited, guide),
                )

                nicking: list[NickingGuide] = []
                for ng in edited_guides:
                    if ng.cut_site is None or ng.strand == guide.strand:
                        continue
                    distance = abs(ng.cut_site - guide.cut_site)
                    if pe3_min_distance <= distance <= pe3_max_distance:
                        nicking.append(NickingGuide(
                            guide=deepcopy(ng),
                            nick_distance=distance,
                            mode=_classify_nick_guide(ref, edited, ng),  # type: ignore[arg-type]
                        ))
                nicking.sort(key=lambda x: (0 if x.mode == 'PE3b' else 1, abs(x.nick_distance - 75)))
                candidate.nicking_guides = nicking[:10]

                # Transparent structural ranking until an optional learned prime-edit scorer is configured.
                structural = heuristic_on_target(guide)
                structural += 0.15 if candidate.pam_disrupted else 0.0
                structural += 0.10 if any(x.mode == 'PE3b' for x in candidate.nicking_guides) else 0.0
                structural -= min(0.25, nick_to_edit / 200)
                candidate.scores['structural_prime'] = round(max(0.0, min(1.0, structural)), 4)
                candidate.scores['model_note'] = (
                    'Structural fallback only. PRIDICT2/DeepPrime adapters can replace this score without changing the API.'
                )
                out.append(candidate)

    out.sort(key=lambda c: (
        c.pam_disrupted,
        any(x.mode == 'PE3b' for x in c.nicking_guides),
        c.scores['structural_prime'],
        -c.nick_to_edit,
        -c.rtt_length,
    ), reverse=True)
    for i, candidate in enumerate(out, 1):
        candidate.rank = i
    return out
