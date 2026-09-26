from __future__ import annotations

from copy import deepcopy
from importlib.metadata import PackageNotFoundError, version
from typing import Iterable

from .models import PrimeCandidate, PrimeEdit
from .prime import apply_edit
from .sequence import normalize_dna


def deepprime_input_sequence(reference: str, edit: PrimeEdit) -> str:
    """Convert a reference/edit pair to GenET DeepPrime 60-(REF/ALT)-61 notation."""
    refseq = normalize_dna(reference)
    ref = normalize_dna(edit.ref) if edit.ref else ''
    alt = normalize_dna(edit.alt) if edit.alt else ''
    apply_edit(refseq, PrimeEdit(edit.position, ref, alt))

    left = refseq[:edit.position]
    right = refseq[edit.position + len(ref):]
    if len(left) < 60 or len(right) < 61:
        raise ValueError(
            'DeepPrime requires at least 60 nt upstream and 61 nt downstream of the edit. '
            'Provide a longer reference context.'
        )
    if ref and alt and len(ref) != len(alt):
        raise ValueError('GenET DeepPrime supports substitutions, insertions, or deletions, not unequal complex replacements')
    edit_len = max(len(ref), len(alt))
    if not 1 <= edit_len <= 3:
        raise ValueError('GenET DeepPrime supports edit lengths of 1-3 nt')

    return left[-60:] + f'({ref}/{alt})' + right[:61]


def _row_key(row: dict) -> tuple[str, int, int, str] | None:
    try:
        spacer = str(row['Spacer']).upper().replace('U', 'T')
        pbs_len = int(row['PBS_len'])
        rtt_len = int(row['RTT_len'])
        extension = str(row['RT-PBS']).upper().replace('U', 'T')
    except (KeyError, TypeError, ValueError):
        return None
    return spacer, pbs_len, rtt_len, extension


def apply_deepprime_rows(
    candidates: Iterable[PrimeCandidate],
    rows: Iterable[dict],
    *,
    score_column: str,
    pe_system: str,
    cell_type: str,
    backend_version: str | None = None,
) -> list[PrimeCandidate]:
    """Pure matching layer used by the real GenET adapter and unit tests."""
    score_map: dict[tuple[str, int, int, str], float] = {}
    for row in rows:
        key = _row_key(row)
        if key is None:
            continue
        try:
            score_map[key] = float(row[score_column])
        except (KeyError, TypeError, ValueError):
            continue

    out: list[PrimeCandidate] = []
    for original in candidates:
        candidate = deepcopy(original)
        key = (
            candidate.spacer.spacer.upper(),
            candidate.pbs_length,
            candidate.rtt_length,
            candidate.extension_sequence.upper(),
        )
        score = score_map.get(key)
        payload = {
            'score': score,
            'pe_system': pe_system,
            'cell_type': cell_type,
            'backend': 'genet.DeepPrime',
            'version': backend_version,
        }
        if score is None:
            payload['note'] = 'candidate was not present in GenET DeepPrime output for this edit/configuration'
        candidate.scores['deepprime'] = payload
        out.append(candidate)
    return out


def score_prime_candidates_deepprime(
    reference: str,
    edit: PrimeEdit,
    candidates: Iterable[PrimeCandidate],
    *,
    pe_system: str = 'PE2max',
    cell_type: str = 'HEK293T',
    pam: str = 'NGG',
) -> list[PrimeCandidate]:
    """Apply the learned DeepPrime model through the optional MIT-licensed GenET package."""
    candidate_list = list(candidates)
    if not candidate_list:
        return []
    notation = deepprime_input_sequence(reference, edit)
    try:
        from genet.predict import DeepPrime  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            'DeepPrime scoring requires GenET. Use a Python 3.10 environment and install '
            'gRNA-lib[deepprime], or keep the structural scorer.'
        ) from exc

    pbs_lengths = [c.pbs_length for c in candidate_list]
    rtt_lengths = [c.rtt_length for c in candidate_list]
    if max(rtt_lengths) > 40:
        raise ValueError('DeepPrime supports RTT lengths up to 40 nt')

    model = DeepPrime(
        notation,
        pam=pam,
        pbs_min=max(1, min(pbs_lengths)),
        pbs_max=min(17, max(pbs_lengths)),
        rtt_min=max(0, min(rtt_lengths)),
        rtt_max=min(40, max(rtt_lengths)),
    )
    frame = model.predict(pe_system=pe_system, cell_type=cell_type)
    score_column = f'{pe_system}_score'
    if score_column not in frame.columns:
        raise RuntimeError(f'GenET DeepPrime output is missing expected column {score_column!r}')
    try:
        backend_version = version('genet')
    except PackageNotFoundError:
        backend_version = None
    return apply_deepprime_rows(
        candidate_list,
        frame.to_dict('records'),
        score_column=score_column,
        pe_system=pe_system,
        cell_type=cell_type,
        backend_version=backend_version,
    )


def rank_prime_candidates(candidates: Iterable[PrimeCandidate]) -> list[PrimeCandidate]:
    """Prefer learned DeepPrime scores when available, then the transparent structural score."""
    ranked = [deepcopy(c) for c in candidates]

    def key(candidate: PrimeCandidate):
        learned = candidate.scores.get('deepprime')
        learned_score = learned.get('score') if isinstance(learned, dict) else None
        has_learned = learned_score is not None
        structural = float(candidate.scores.get('structural_prime', 0.0))
        return (
            1 if has_learned else 0,
            float(learned_score) if has_learned else float('-inf'),
            structural,
            1 if candidate.pam_disrupted else 0,
            1 if any(x.mode == 'PE3b' for x in candidate.nicking_guides) else 0,
            -candidate.nick_to_edit,
        )

    ranked.sort(key=key, reverse=True)
    for i, candidate in enumerate(ranked, 1):
        candidate.rank = i
    return ranked
