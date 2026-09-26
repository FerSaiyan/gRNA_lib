from __future__ import annotations

from importlib.resources import files

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from pydantic import BaseModel, Field

from . import (
    PrimeEdit,
    design_guides,
    design_prime_edit,
    rank_guides,
    rank_prime_candidates,
    resolve_nuclease,
    score_prime_candidates_deepprime,
)
from .backends import profile_with_crispr_ots, score_with_crisprware
from .genome_sources import download_ncbi_genome, fetch_ensembl_region, read_fasta_record

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


class GenomeProfileRequest(BaseModel):
    sequence: str
    crispr_ots_index: str
    chromosome: str
    reference_start: int = 0
    nuclease: str = 'SpCas9'
    pam: str | None = None
    spacer_length: int | None = None
    pam_side: str | None = None
    threads: int = 4
    mismatches: int = 4
    hit_limit: int = 100
    crispr_ots_executable: str = 'crispr-ots'


class NcbiSourceRequest(BaseModel):
    accession: str
    chromosomes: list[str] = Field(default_factory=list)
    record: str | None = None
    inline_limit: int = 10_000_000


class EnsemblSourceRequest(BaseModel):
    species: str
    region: str


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


@app.post('/api/guides/profile')
def api_profile(request: GenomeProfileRequest):
    try:
        spec = resolve_nuclease(
            request.nuclease,
            pam=request.pam,
            spacer_length=request.spacer_length,
            pam_side=request.pam_side,
        )
        guides = design_guides(request.sequence, spec)
        result = profile_with_crispr_ots(
            guides,
            spec,
            index=request.crispr_ots_index,
            chromosome=request.chromosome,
            reference_start=request.reference_start,
            executable=request.crispr_ots_executable,
            threads=request.threads,
            mismatches=request.mismatches,
            hit_limit=request.hit_limit,
        )
        return {'nuclease': spec.__dict__, **result.to_dict()}
    except (ValueError, RuntimeError, FileNotFoundError) as exc:
        raise HTTPException(400, str(exc)) from exc


@app.post('/api/sources/ncbi')
def api_source_ncbi(request: NcbiSourceRequest):
    try:
        result = download_ncbi_genome(
            request.accession,
            chromosomes=request.chromosomes,
        )
        payload = result.to_dict()
        selected = request.record
        if selected is None and len(result.records) == 1:
            selected = result.records[0].name
        if selected is not None:
            info = next((item for item in result.records if item.name == selected), None)
            if info is None:
                raise ValueError(f'FASTA record {selected!r} is not present in the downloaded package')
            if info.length <= request.inline_limit:
                _, sequence = read_fasta_record(result.fasta_path, selected)
                payload['selected_record'] = selected
                payload['sequence'] = sequence
            else:
                payload['selected_record'] = selected
                payload['sequence_omitted'] = (
                    f'record is {info.length:,} bp; it was cached on disk instead of copied into the browser'
                )
        return payload
    except (ValueError, RuntimeError, FileNotFoundError) as exc:
        raise HTTPException(400, str(exc)) from exc


@app.post('/api/sources/ensembl')
def api_source_ensembl(request: EnsemblSourceRequest):
    try:
        sequence = fetch_ensembl_region(request.species, request.region)
        return {
            'source': 'Ensembl REST',
            'species': request.species,
            'region': request.region,
            'length': len(sequence),
            'sequence': sequence,
        }
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
