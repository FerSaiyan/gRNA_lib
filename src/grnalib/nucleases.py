from __future__ import annotations

from .models import NucleaseSpec


NUCLEASES: dict[str, NucleaseSpec] = {
    'SpCas9': NucleaseSpec('SpCas9', 'NGG', 20, '3prime', 3),
    'SpCas9-NG': NucleaseSpec('SpCas9-NG', 'NG', 20, '3prime', 3),
    'SaCas9': NucleaseSpec('SaCas9', 'NNGRRT', 21, '3prime', 3),
    'AsCas12a': NucleaseSpec('AsCas12a', 'TTTV', 23, '5prime', 18),
}


def resolve_nuclease(
    name: str = 'SpCas9',
    *,
    pam: str | None = None,
    spacer_length: int | None = None,
    pam_side: str | None = None,
) -> NucleaseSpec:
    base = NUCLEASES.get(name, NUCLEASES['SpCas9'])
    if pam is None and spacer_length is None and pam_side is None:
        return base
    side = pam_side or base.pam_side
    if side not in {'3prime', '5prime'}:
        raise ValueError("pam_side must be '3prime' or '5prime'")
    effective_name = name if name in NUCLEASES and (pam is None or pam.upper() == base.pam) else 'Custom'
    return NucleaseSpec(
        name=effective_name,
        pam=(pam or base.pam).upper(),
        spacer_length=spacer_length or base.spacer_length,
        pam_side=side,  # type: ignore[arg-type]
        cut_offset=base.cut_offset,
    )
