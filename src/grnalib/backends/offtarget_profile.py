from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from ..models import GuideCandidate, NucleaseSpec
from .crisprware import _guidescan_position


@dataclass(frozen=True)
class OffTargetHit:
    chromosome: str
    position: int
    strand: str
    mismatches: int

    def to_dict(self) -> dict:
        return {
            'chromosome': self.chromosome,
            'position': self.position,
            'strand': self.strand,
            'mismatches': self.mismatches,
        }


@dataclass
class WholeGenomeOffTargetProfile:
    guide_id: str
    spacer: str
    specificity: float | None
    max_mismatches: int
    mismatch_counts: dict[int, int]
    total_off_targets: int
    exact_duplicate_sites: int
    intended_target_found: bool
    representative_hits: list[OffTargetHit] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            'guide_id': self.guide_id,
            'spacer': self.spacer,
            'specificity': self.specificity,
            'max_mismatches': self.max_mismatches,
            'mismatch_counts': {
                str(mm): self.mismatch_counts.get(mm, 0)
                for mm in range(self.max_mismatches + 1)
            },
            'total_off_targets': self.total_off_targets,
            'exact_duplicate_sites': self.exact_duplicate_sites,
            'intended_target_found': self.intended_target_found,
            'representative_hits': [hit.to_dict() for hit in self.representative_hits],
            'notes': list(self.notes),
        }


@dataclass
class WholeGenomeProfileResult:
    profiles: dict[str, WholeGenomeOffTargetProfile]
    backend: str
    index: str
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            'backend': self.backend,
            'index': self.index,
            'profiles': [profile.to_dict() for profile in self.profiles.values()],
            'notes': list(self.notes),
        }


def build_crisprots_kmers(
    guides: Iterable[GuideCandidate],
    spec: NucleaseSpec,
    *,
    chromosome: str,
    reference_start: int = 0,
) -> str:
    """Build the GuideScan2-compatible one-column kmers file used by crispr-ots."""
    if reference_start < 0:
        raise ValueError('reference_start must be >= 0')
    if not chromosome or ',' in chromosome or '\t' in chromosome:
        raise ValueError('chromosome must be a non-empty name without commas/tabs')

    lines = ['id,sequence,pam,chromosome,position,sense']
    for guide in guides:
        position = _guidescan_position(guide, spec, reference_start)
        lines.append(
            ','.join([
                guide.id,
                guide.spacer,
                spec.pam,
                chromosome,
                str(position),
                guide.strand,
            ])
        )
    return '\n'.join(lines) + '\n'


def _resolve_crisprots_executable(executable: str) -> str:
    if '/' in executable or '\\' in executable:
        path = Path(executable)
        if not path.exists():
            raise FileNotFoundError(f'crispr-ots executable not found: {executable}')
        return str(path)
    resolved = shutil.which(executable)
    if resolved is None:
        raise FileNotFoundError(
            f'{executable!r} is not on PATH. Install the crispr-ots binary or pass its path explicitly.'
        )
    return resolved


def _expected_target(
    guide: GuideCandidate,
    spec: NucleaseSpec,
    *,
    chromosome: str,
    reference_start: int,
) -> tuple[str, int, str]:
    # crispr-ots CSV reports a 0-based genomic start. GuideScan-style input
    # stores the same anchor as a 1-based position.
    return (
        chromosome,
        _guidescan_position(guide, spec, reference_start) - 1,
        guide.strand,
    )


def parse_crisprots_csv(
    path: str | Path,
    guides: Iterable[GuideCandidate],
    spec: NucleaseSpec,
    *,
    chromosome: str,
    reference_start: int = 0,
    max_mismatches: int = 4,
    hit_limit: int = 100,
) -> dict[str, WholeGenomeOffTargetProfile]:
    """Summarize one-row-per-hit crispr-ots CSV output into per-guide profiles.

    The intended genomic hit is excluded only when chromosome, 0-based start,
    strand and a zero mismatch count all match the supplied target coordinates.
    Other zero-mismatch genomic copies remain classified as exact duplicates.
    """
    if hit_limit < 0:
        raise ValueError('hit_limit must be >= 0')

    guide_list = list(guides)
    by_id = {guide.id: guide for guide in guide_list}
    expected = {
        guide.id: _expected_target(
            guide,
            spec,
            chromosome=chromosome,
            reference_start=reference_start,
        )
        for guide in guide_list
    }
    counts = {guide.id: Counter() for guide in guide_list}
    specificities: dict[str, float | None] = {guide.id: None for guide in guide_list}
    intended_found = {guide.id: False for guide in guide_list}
    hits: dict[str, list[OffTargetHit]] = {guide.id: [] for guide in guide_list}
    rows_seen = {guide.id: 0 for guide in guide_list}

    with Path(path).open(newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh)
        required = {
            'id',
            'match_chrm',
            'match_position',
            'match_strand',
            'match_distance',
            'specificity',
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                'crispr-ots CSV is missing required columns: '
                + ', '.join(sorted(required))
            )

        for row in reader:
            guide_id = row.get('id', '')
            guide = by_id.get(guide_id)
            if guide is None:
                continue
            rows_seen[guide_id] += 1
            try:
                position = int(row['match_position'])
                mismatches = int(row['match_distance'])
            except (TypeError, ValueError) as exc:
                raise ValueError(f'invalid crispr-ots coordinate/mismatch row for {guide_id}') from exc

            try:
                score = float(row['specificity'])
                if specificities[guide_id] is None:
                    specificities[guide_id] = score
            except (TypeError, ValueError):
                pass

            locus = (
                row.get('match_chrm', ''),
                position,
                row.get('match_strand', ''),
            )
            if mismatches == 0 and locus == expected[guide_id]:
                intended_found[guide_id] = True
                continue

            if 0 <= mismatches <= max_mismatches:
                counts[guide_id][mismatches] += 1
            if hit_limit:
                hits[guide_id].append(
                    OffTargetHit(
                        chromosome=locus[0],
                        position=position,
                        strand=locus[2],
                        mismatches=mismatches,
                    )
                )
                # Bound memory for repetitive guides while keeping the nearest
                # mismatch classes. Trimming in chunks avoids sorting every row.
                if len(hits[guide_id]) > max(hit_limit * 4, 1000):
                    hits[guide_id] = sorted(
                        hits[guide_id],
                        key=lambda h: (h.mismatches, h.chromosome, h.position, h.strand),
                    )[:hit_limit]

    profiles: dict[str, WholeGenomeOffTargetProfile] = {}
    for guide in guide_list:
        guide_hits = sorted(
            hits[guide.id],
            key=lambda h: (h.mismatches, h.chromosome, h.position, h.strand),
        )[:hit_limit]
        mismatch_counts = {
            mm: int(counts[guide.id].get(mm, 0))
            for mm in range(max_mismatches + 1)
        }
        notes: list[str] = []
        if rows_seen[guide.id] == 0:
            notes.append('crispr-ots returned no rows for this guide')
        if rows_seen[guide.id] and not intended_found[guide.id]:
            notes.append(
                'intended target row was not identified by exact chromosome/start/strand match; '
                'zero-mismatch copies were therefore retained as off-target evidence'
            )
        total = sum(mismatch_counts.values())
        profiles[guide.id] = WholeGenomeOffTargetProfile(
            guide_id=guide.id,
            spacer=guide.spacer,
            specificity=specificities[guide.id],
            max_mismatches=max_mismatches,
            mismatch_counts=mismatch_counts,
            total_off_targets=total,
            exact_duplicate_sites=mismatch_counts.get(0, 0),
            intended_target_found=intended_found[guide.id],
            representative_hits=guide_hits,
            notes=notes,
        )
    return profiles


def profile_with_crispr_ots(
    guides: Iterable[GuideCandidate],
    spec: NucleaseSpec,
    *,
    index: str | Path,
    chromosome: str,
    reference_start: int = 0,
    executable: str = 'crispr-ots',
    threads: int = 4,
    mismatches: int = 4,
    hit_limit: int = 100,
) -> WholeGenomeProfileResult:
    """Enumerate whole-genome off-target loci and mismatch classes with crispr-ots.

    This uses crispr-ots' GuideScan2-compatible CSV output because it retains
    one row per genomic hit. Bulges are not requested here: current crispr-ots
    bin scanning supports mismatch enumeration but rejects non-zero bulges.
    """
    guide_list = list(guides)
    if not guide_list:
        return WholeGenomeProfileResult({}, 'crispr-ots', str(index), notes=['no guides supplied'])
    if mismatches < 0:
        raise ValueError('mismatches must be >= 0')

    exe = _resolve_crisprots_executable(executable)
    index = str(index)
    with tempfile.TemporaryDirectory(prefix='grnalib-crisprots-profile-') as tmp:
        tmpdir = Path(tmp)
        kmers = tmpdir / 'guides.kmers.tsv'
        output = tmpdir / 'offtargets.csv'
        kmers.write_text(
            build_crisprots_kmers(
                guide_list,
                spec,
                chromosome=chromosome,
                reference_start=reference_start,
            ),
            encoding='utf-8',
        )
        cmd = [
            exe,
            'enumerate',
            '--kmers-file',
            str(kmers),
            '--mismatches',
            str(mismatches),
            '--threads',
            str(threads),
            '--format',
            'csv',
            '--spec-convention',
            'guidescan',
            '--threshold',
            '-1',
            '--max-off-targets-per-bin',
            '-1',
            '--mode',
            'complete',
            '--output',
            str(output),
            index,
        ]
        proc = subprocess.run(cmd, text=True, capture_output=True)
        if proc.returncode != 0:
            detail = (proc.stderr or proc.stdout or '').strip()[-5000:]
            raise RuntimeError(f'crispr-ots profiling failed ({proc.returncode}): {detail}')
        if not output.exists():
            raise RuntimeError('crispr-ots completed but did not create the expected CSV output')
        profiles = parse_crisprots_csv(
            output,
            guide_list,
            spec,
            chromosome=chromosome,
            reference_start=reference_start,
            max_mismatches=mismatches,
            hit_limit=hit_limit,
        )

    return WholeGenomeProfileResult(
        profiles=profiles,
        backend='crispr-ots',
        index=index,
        notes=[
            'mismatch counts exclude only the coordinate-matched intended target',
            'representative_hits are ordered by mismatch count, then genomic coordinate',
            'current crispr-ots mismatch profiling does not support RNA/DNA bulges',
        ],
    )
