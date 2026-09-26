from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Literal


Strand = Literal['+', '-']


@dataclass(frozen=True)
class NucleaseSpec:
    name: str
    pam: str
    spacer_length: int = 20
    pam_side: Literal['3prime', '5prime'] = '3prime'
    cut_offset: int = 3


@dataclass
class GuideCandidate:
    spacer: str
    expressed_spacer: str
    pam: str
    strand: Strand
    start: int
    end: int
    pam_start: int
    pam_end: int
    cut_site: int | None
    context: str | None = None
    gc_fraction: float = 0.0
    warnings: list[str] = field(default_factory=list)
    scores: dict[str, Any] = field(default_factory=dict)
    rank: int | None = None

    @property
    def id(self) -> str:
        return f"{self.strand}:{self.start}-{self.end}:{self.pam}"

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data['id'] = self.id
        return data


@dataclass
class OffTargetSummary:
    exact_other_sites: int = 0
    mismatches_1: int = 0
    mismatches_2: int = 0
    mismatches_3: int = 0
    mismatches_4: int = 0
    specificity_proxy: float = 1.0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class PrimeEdit:
    position: int
    ref: str
    alt: str


@dataclass
class NickingGuide:
    guide: GuideCandidate
    nick_distance: int
    mode: Literal['PE3', 'PE3b']

    def to_dict(self) -> dict[str, Any]:
        return {
            'guide': self.guide.to_dict(),
            'nick_distance': self.nick_distance,
            'mode': self.mode,
        }


@dataclass
class PrimeCandidate:
    spacer: GuideCandidate
    pbs_sequence: str
    pbs_length: int
    rtt_sequence: str
    rtt_length: int
    extension_sequence: str
    nick_to_edit: int
    pam_disrupted: bool
    nicking_guides: list[NickingGuide] = field(default_factory=list)
    scores: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    rank: int | None = None

    @property
    def id(self) -> str:
        return f"{self.spacer.id}:PBS{self.pbs_length}:RTT{self.rtt_length}"

    def to_dict(self) -> dict[str, Any]:
        return {
            'id': self.id,
            'spacer': self.spacer.to_dict(),
            'pbs_sequence': self.pbs_sequence,
            'pbs_length': self.pbs_length,
            'rtt_sequence': self.rtt_sequence,
            'rtt_length': self.rtt_length,
            'extension_sequence': self.extension_sequence,
            'nick_to_edit': self.nick_to_edit,
            'pam_disrupted': self.pam_disrupted,
            'nicking_guides': [x.to_dict() for x in self.nicking_guides],
            'scores': self.scores,
            'warnings': self.warnings,
            'rank': self.rank,
        }
