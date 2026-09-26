from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from ..models import GuideCandidate, NucleaseSpec


@dataclass
class IndexedSpecificityResult:
    scores: dict[str, float]
    backend: str
    index: str
    specificity_column: str | None = None
    notes: list[str] = field(default_factory=list)


def _guidescan_position(
    guide: GuideCandidate,
    spec: NucleaseSpec,
    reference_start: int,
) -> int:
    """Return the 1-based position CRISPRware/GuideScan uses in its composite ID."""
    if spec.pam_side == '5prime':
        anchor = guide.pam_start
    elif guide.strand == '+':
        anchor = guide.start
    else:
        anchor = guide.pam_start
    return reference_start + anchor + 1


def build_crisprware_bed(
    guides: Iterable[GuideCandidate],
    spec: NucleaseSpec,
    *,
    chromosome: str,
    reference_start: int = 0,
) -> str:
    """Build the six-column guide BED expected by CRISPRware score_guides.

    reference_start is the 0-based genomic coordinate of base 0 in the sequence
    originally passed to design_guides.
    """
    if reference_start < 0:
        raise ValueError('reference_start must be >= 0')
    if not chromosome or ',' in chromosome or '\t' in chromosome:
        raise ValueError('chromosome must be a non-empty name without commas/tabs')

    header = '#chr\tstart\tstop\tid,sequence,pam,chromosome,position,sense\tcontext\tstrand'
    lines = [header]
    for guide in guides:
        start = reference_start + guide.start
        stop = reference_start + guide.end
        position = _guidescan_position(guide, spec, reference_start)
        ident = f'{chromosome}:{position}:{guide.strand}'
        composite = ','.join([
            ident,
            guide.spacer,
            spec.pam,
            chromosome,
            str(position),
            guide.strand,
        ])
        context = guide.context or (guide.spacer + guide.pam)
        lines.append(
            f'{chromosome}\t{start}\t{stop}\t{composite}\t{context}\t{guide.strand}'
        )
    return '\n'.join(lines) + '\n'


def _resolve_executable(executable: str) -> str:
    if '/' in executable or '\\' in executable:
        path = Path(executable)
        if not path.exists():
            raise FileNotFoundError(f'CRISPRware executable not found: {executable}')
        return str(path)
    resolved = shutil.which(executable)
    if resolved is None:
        raise FileNotFoundError(
            f'{executable!r} is not on PATH. Install CRISPRware separately or pass --crisprware-executable.'
        )
    return resolved


def score_with_crisprware(
    guides: Iterable[GuideCandidate],
    spec: NucleaseSpec,
    *,
    index: str | Path,
    chromosome: str,
    reference_start: int = 0,
    executable: str = 'crisprware',
    threads: int = 4,
    mismatches: int = 3,
    rna_bulges: int = 0,
    dna_bulges: int = 0,
) -> IndexedSpecificityResult:
    """Score candidates against an installed crispr-ots or GuideScan2 index.

    CRISPRware is intentionally an external integration: no CRISPRware source,
    model, or index is vendored in this project.
    """
    guide_list = list(guides)
    if not guide_list:
        return IndexedSpecificityResult({}, 'crisprware', str(index), notes=['no guides supplied'])
    exe = _resolve_executable(executable)
    index = str(index)

    with tempfile.TemporaryDirectory(prefix='grnalib-crisprware-') as tmp:
        tmpdir = Path(tmp)
        bed_path = tmpdir / 'guides.bed'
        outdir = tmpdir / 'out'
        outdir.mkdir()
        bed_path.write_text(
            build_crisprware_bed(
                guide_list,
                spec,
                chromosome=chromosome,
                reference_start=reference_start,
            ),
            encoding='utf-8',
        )
        cmd = [
            exe,
            'score_guides',
            '-b',
            str(bed_path),
            '-i',
            index,
            '--skip_rs3',
            '--threshold',
            '-1',
            '--mismatches',
            str(mismatches),
            '--rna_bulges',
            str(rna_bulges),
            '--dna_bulges',
            str(dna_bulges),
            '--threads',
            str(threads),
            '--drop_duplicates',
            '-o',
            str(outdir),
        ]
        proc = subprocess.run(cmd, text=True, capture_output=True)
        if proc.returncode != 0:
            detail = (proc.stderr or proc.stdout or '').strip()[-5000:]
            raise RuntimeError(f'CRISPRware scoring failed ({proc.returncode}): {detail}')

        outputs = list(outdir.rglob('*scoredgRNA*.bed'))
        if not outputs:
            outputs = list(outdir.rglob('*.bed'))
        if not outputs:
            raise RuntimeError('CRISPRware completed but no scored BED was found')
        scored_path = max(outputs, key=lambda p: p.stat().st_mtime)

        with scored_path.open(newline='', encoding='utf-8') as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        if not rows:
            return IndexedSpecificityResult(
                {}, 'crisprware', index, notes=['CRISPRware returned an empty scored BED']
            )
        specificity_cols = [c for c in rows[0] if c.startswith('specificity_')]
        if not specificity_cols:
            raise RuntimeError('CRISPRware output did not contain a specificity_* column')
        column = specificity_cols[0]

        by_spacer: dict[str, float] = {}
        for row in rows:
            try:
                score = float(row[column])
            except (TypeError, ValueError):
                continue
            if score < 0:
                continue
            by_spacer[row.get('sequence', '').upper()] = score

    scores = {
        guide.id: by_spacer[guide.spacer]
        for guide in guide_list
        if guide.spacer in by_spacer
    }
    notes: list[str] = []
    missing = len(guide_list) - len(scores)
    if missing:
        notes.append(f'{missing} guide(s) had no indexed specificity result')
    return IndexedSpecificityResult(
        scores=scores,
        backend='crisprware',
        index=index,
        specificity_column=column,
        notes=notes,
    )


def build_crisprware_index(
    fasta: str | Path,
    *,
    pam: str = 'NGG',
    spacer_length: int = 20,
    pam_5_prime: bool = False,
    bin_width: int | None = None,
    output_directory: str | Path = '.',
    executable: str = 'crisprware',
) -> None:
    """Invoke the CRISPRware index_genome command for a FASTA."""
    exe = _resolve_executable(executable)
    cmd = [
        exe,
        'index_genome',
        '-f',
        str(fasta),
        '-p',
        pam,
        '-l',
        str(spacer_length),
        '-o',
        str(output_directory),
    ]
    if pam_5_prime:
        cmd.append('--pam_5_prime')
    if bin_width is not None:
        if not 1 <= bin_width <= 15:
            raise ValueError('bin_width must be between 1 and 15')
        cmd.extend(['--bin_width', str(bin_width)])
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or '').strip()[-5000:]
        raise RuntimeError(f'CRISPRware index build failed ({proc.returncode}): {detail}')
