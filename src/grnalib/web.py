from __future__ import annotations

from importlib.resources import files

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from pydantic import BaseModel

from . import (
    PrimeEdit,
    design_guides,
    design_prime_edit,
    rank_guides,
    rank_prime_candidates,
    resolve_nuclease,
    score_prime_candidates_deepprime,
)
from .backends import score_with_crisprware

app = FastAPI(title='gRNA Library')


class GuideRequest(BaseModel):
    sequence: str
    genome_sequence: str | None = None
    nuclease: str = 'SpCas9'
    pam: str | None = None
    spacer_length: int | None = None
    pam_side: str | None = None
    crisprware_index: str | None = None
    chromosome: str | None = None
    reference_start: int = 0
    threads: int = 4
    mismatches: int = 3
    rna_bulges: int = 0
    dna_bulges: int = 0


class PrimeRequest(BaseModel):
    sequence: str
    position: int
    ref: str = ''
    alt: str = ''
    nuclease: str = 'SpCas9'
    pam: str | None = None
    limit: int = 200
    scorer: str = 'structural'
    pe_system: str = 'PE2max'
    cell_type: str = 'HEK293T'


@app.get('/', response_class=HTMLResponse)
def index() -> str:
    return files('grnalib').joinpath('static/index.html').read_text(encoding='utf-8')


@app.post('/api/guides/design')
def api_design(request: GuideRequest):
    try:
        spec = resolve_nuclease(
            request.nuclease, pam=request.pam, spacer_length=request.spacer_length, pam_side=request.pam_side
        )
        guides = design_guides(request.sequence, spec)
        return {'nuclease': spec.__dict__, 'guides': [g.to_dict() for g in guides]}
    except (ValueError, RuntimeError, FileNotFoundError) as exc:
        raise HTTPException(400, str(exc)) from exc


@app.post('/api/guides/rank')
def api_rank(request: GuideRequest):
    try:
        spec = resolve_nuclease(
            request.nuclease, pam=request.pam, spacer_length=request.spacer_length, pam_side=request.pam_side
        )
        guides = design_guides(request.sequence, spec)
        indexed = None
        source = None
        if request.crisprware_index:
            if not request.chromosome:
                raise ValueError('chromosome is required when crisprware_index is supplied')
            result = score_with_crisprware(
                guides,
                spec,
                index=request.crisprware_index,
                chromosome=request.chromosome,
                reference_start=request.reference_start,
                threads=request.threads,
                mismatches=request.mismatches,
                rna_bulges=request.rna_bulges,
                dna_bulges=request.dna_bulges,
            )
            indexed = result.scores
            source = f'{result.backend}:{result.index}:{result.specificity_column}'
        ranked = rank_guides(
            guides,
            spec=spec,
            genome_sequence=request.genome_sequence,
            indexed_specificity=indexed,
            indexed_specificity_source=source,
        )
        return {'nuclease': spec.__dict__, 'guides': [g.to_dict() for g in ranked]}
    except (ValueError, RuntimeError, FileNotFoundError) as exc:
        raise HTTPException(400, str(exc)) from exc


@app.post('/api/prime/design')
def api_prime(request: PrimeRequest):
    try:
        spec = resolve_nuclease(request.nuclease, pam=request.pam)
        edit = PrimeEdit(request.position, request.ref, request.alt)
        candidates = design_prime_edit(request.sequence, edit, spec)
        if request.scorer.lower() == 'deepprime':
            candidates = score_prime_candidates_deepprime(
                request.sequence,
                edit,
                candidates,
                pe_system=request.pe_system,
                cell_type=request.cell_type,
                pam=spec.pam,
            )
        elif request.scorer.lower() != 'structural':
            raise ValueError('scorer must be structural or deepprime')
        candidates = rank_prime_candidates(candidates)
        return {'nuclease': spec.__dict__, 'candidates': [x.to_dict() for x in candidates[:request.limit]]}
    except (ValueError, RuntimeError, FileNotFoundError) as exc:
        raise HTTPException(400, str(exc)) from exc
