from __future__ import annotations

from importlib.resources import files

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from pydantic import BaseModel

from . import PrimeEdit, design_guides, design_prime_edit, rank_guides, resolve_nuclease

app = FastAPI(title='gRNA Library')


class GuideRequest(BaseModel):
    sequence: str
    genome_sequence: str | None = None
    nuclease: str = 'SpCas9'
    pam: str | None = None
    spacer_length: int | None = None
    pam_side: str | None = None


class PrimeRequest(BaseModel):
    sequence: str
    position: int
    ref: str = ''
    alt: str = ''
    nuclease: str = 'SpCas9'
    pam: str | None = None
    limit: int = 200


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
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc


@app.post('/api/guides/rank')
def api_rank(request: GuideRequest):
    try:
        spec = resolve_nuclease(
            request.nuclease, pam=request.pam, spacer_length=request.spacer_length, pam_side=request.pam_side
        )
        guides = design_guides(request.sequence, spec)
        ranked = rank_guides(guides, spec=spec, genome_sequence=request.genome_sequence)
        return {'nuclease': spec.__dict__, 'guides': [g.to_dict() for g in ranked]}
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc


@app.post('/api/prime/design')
def api_prime(request: PrimeRequest):
    try:
        spec = resolve_nuclease(request.nuclease, pam=request.pam)
        candidates = design_prime_edit(request.sequence, PrimeEdit(request.position, request.ref, request.alt), spec)
        return {'nuclease': spec.__dict__, 'candidates': [x.to_dict() for x in candidates[:request.limit]]}
    except ValueError as exc:
        raise HTTPException(400, str(exc)) from exc
