from copy import deepcopy

import pytest

from grnalib import PrimeEdit, design_guides, design_prime_edit, rank_guides, resolve_nuclease
from grnalib.backends import crisprware
from grnalib.backends.crisprware import build_crisprware_bed
from grnalib.backends.offtarget_profile import build_crisprots_kmers, parse_crisprots_csv
from grnalib.prime_scoring import (
    apply_deepprime_rows,
    deepprime_input_sequence,
    rank_prime_candidates,
)


def _prime_candidates():
    seq = 'AAAA' + 'GCGCGCGCGCGCGCGCGCGC' + 'TGG' + 'ACGTACGTACGTACGTACGTACGTACGTACGT'
    edit = PrimeEdit(position=23, ref=seq[23], alt='A' if seq[23] != 'A' else 'C')
    candidates = design_prime_edit(
        seq,
        edit,
        resolve_nuclease(),
        pbs_lengths=[10],
        rtt_lengths=[10, 15],
    )
    assert candidates
    return candidates


def test_deepprime_input_sequence_uses_60_61_context():
    reference = 'A' * 60 + 'C' + 'G' * 61
    notation = deepprime_input_sequence(reference, PrimeEdit(60, 'C', 'T'))
    assert notation == 'A' * 60 + '(C/T)' + 'G' * 61


def test_deepprime_input_requires_context():
    with pytest.raises(ValueError, match='60 nt upstream'):
        deepprime_input_sequence('A' * 20 + 'C' + 'G' * 80, PrimeEdit(20, 'C', 'T'))


def test_deepprime_rows_match_our_candidate_schema():
    candidate = _prime_candidates()[0]
    rows = [{
        'Spacer': candidate.spacer.spacer,
        'PBS_len': candidate.pbs_length,
        'RTT_len': candidate.rtt_length,
        'RT-PBS': candidate.extension_sequence,
        'PE2max_score': 12.5,
    }]
    scored = apply_deepprime_rows(
        [candidate],
        rows,
        score_column='PE2max_score',
        pe_system='PE2max',
        cell_type='HEK293T',
        backend_version='test',
    )
    assert scored[0].scores['deepprime']['score'] == 12.5
    assert scored[0].scores['deepprime']['backend'] == 'genet.DeepPrime'


def test_prime_ranking_prefers_learned_score_when_present():
    candidates = _prime_candidates()
    a = deepcopy(candidates[0])
    b = deepcopy(candidates[-1])
    a.scores['structural_prime'] = 1.0
    b.scores['structural_prime'] = 0.1
    b.scores['deepprime'] = {
        'score': 2.0,
        'pe_system': 'PE2max',
        'cell_type': 'HEK293T',
        'backend': 'genet.DeepPrime',
        'version': 'test',
    }
    ranked = rank_prime_candidates([a, b])
    assert ranked[0].id == b.id


def test_crisprware_bed_uses_true_genomic_offset_for_plus_guide():
    seq = 'AAAA' + 'ACACACACACACACACACAC' + 'AGG' + 'AAAA'
    spec = resolve_nuclease()
    guide = next(g for g in design_guides(seq, spec) if g.strand == '+')
    bed = build_crisprware_bed([guide], spec, chromosome='chr7', reference_start=1000)
    line = bed.strip().splitlines()[1].split('\t')
    assert line[0] == 'chr7'
    assert int(line[1]) == 1000 + guide.start
    composite = line[3].split(',')
    assert int(composite[4]) == 1000 + guide.start + 1
    assert composite[5] == '+'


def test_crisprware_bed_minus_uses_pam_anchor():
    seq = 'CCATTTTTTTTTTTTTTTTTTT' + 'AAAAAA'
    spec = resolve_nuclease()
    guide = next(g for g in design_guides(seq, spec) if g.strand == '-')
    bed = build_crisprware_bed([guide], spec, chromosome='chr1', reference_start=500)
    composite = bed.strip().splitlines()[1].split('\t')[3].split(',')
    assert int(composite[4]) == 500 + guide.pam_start + 1
    assert composite[5] == '-'


def test_indexed_specificity_controls_ranking():
    seq = (
        'AAAA' + 'GCGCGCGCGCGCGCGCGCGC' + 'TGG' +
        'AAAA' + 'ACACACACACACACACACAC' + 'AGG' + 'AAAA'
    )
    spec = resolve_nuclease()
    guides = [g for g in design_guides(seq, spec) if g.strand == '+']
    assert len(guides) >= 2
    scores = {guides[0].id: 0.1, guides[1].id: 0.9}
    ranked = rank_guides(
        guides[:2],
        spec=spec,
        use_rs3=False,
        indexed_specificity=scores,
        indexed_specificity_source='test-index',
    )
    assert ranked[0].id == guides[1].id
    assert ranked[0].scores['specificity_used_for_ranking'] == 'indexed'


def test_missing_indexed_score_is_not_treated_as_perfect():
    seq = (
        'AAAA' + 'GCGCGCGCGCGCGCGCGCGC' + 'TGG' +
        'AAAA' + 'ACACACACACACACACACAC' + 'AGG' + 'AAAA'
    )
    spec = resolve_nuclease()
    guides = [g for g in design_guides(seq, spec) if g.strand == '+'][:2]
    assert len(guides) == 2
    scores = {guides[0].id: 0.2}
    ranked = rank_guides(
        guides,
        spec=spec,
        use_rs3=False,
        indexed_specificity=scores,
        indexed_specificity_source='test-index',
    )
    assert ranked[0].id == guides[0].id
    missing = next(g for g in ranked if g.id == guides[1].id)
    assert missing.scores['specificity_used_for_ranking'] == 'indexed-missing'


def test_crisprots_kmers_uses_guide_id_and_true_genomic_position():
    seq = 'AAAA' + 'ACACACACACACACACACAC' + 'AGG' + 'AAAA'
    spec = resolve_nuclease()
    guide = next(g for g in design_guides(seq, spec) if g.strand == '+')
    text = build_crisprots_kmers([guide], spec, chromosome='chr7', reference_start=1000)
    lines = text.strip().splitlines()
    assert lines[0] == 'id,sequence,pam,chromosome,position,sense'
    fields = lines[1].split(',')
    assert fields[0] == guide.id
    assert fields[1] == guide.spacer
    assert fields[2] == 'NGG'
    assert fields[3] == 'chr7'
    assert int(fields[4]) == 1000 + guide.start + 1
    assert fields[5] == '+'


def test_crisprots_profile_excludes_only_coordinate_matched_intended_target(tmp_path):
    seq = 'AAAA' + 'ACACACACACACACACACAC' + 'AGG' + 'AAAA'
    spec = resolve_nuclease()
    guide = next(g for g in design_guides(seq, spec) if g.strand == '+')
    intended = 1000 + guide.start
    csv_path = tmp_path / 'offtargets.csv'
    csv_path.write_text(
        'id,sequence,match_chrm,match_position,match_strand,match_distance,specificity\n'
        f'{guide.id},{guide.spacer}NGG,chr7,{intended},+,0,0.42\n'
        f'{guide.id},{guide.spacer}NGG,chr7,{intended + 100},+,0,0.42\n'
        f'{guide.id},{guide.spacer}NGG,chr8,2000,-,1,0.42\n'
        f'{guide.id},{guide.spacer}NGG,chr2,3000,+,2,0.42\n'
        f'{guide.id},{guide.spacer}NGG,chr3,4000,-,2,0.42\n',
        encoding='utf-8',
    )
    profiles = parse_crisprots_csv(
        csv_path,
        [guide],
        spec,
        chromosome='chr7',
        reference_start=1000,
        max_mismatches=3,
        hit_limit=2,
    )
    profile = profiles[guide.id]
    assert profile.intended_target_found is True
    assert profile.specificity == 0.42
    assert profile.mismatch_counts == {0: 1, 1: 1, 2: 2, 3: 0}
    assert profile.exact_duplicate_sites == 1
    assert profile.total_off_targets == 4
    assert len(profile.representative_hits) == 2
    assert [hit.mismatches for hit in profile.representative_hits] == [0, 1]


def test_crisprots_profile_keeps_zero_mismatch_when_intended_locus_not_found(tmp_path):
    seq = 'AAAA' + 'ACACACACACACACACACAC' + 'AGG' + 'AAAA'
    spec = resolve_nuclease()
    guide = next(g for g in design_guides(seq, spec) if g.strand == '+')
    csv_path = tmp_path / 'offtargets.csv'
    csv_path.write_text(
        'id,sequence,match_chrm,match_position,match_strand,match_distance,specificity\n'
        f'{guide.id},{guide.spacer}NGG,chrUn,999,+,0,0.2\n',
        encoding='utf-8',
    )
    profile = parse_crisprots_csv(
        csv_path,
        [guide],
        spec,
        chromosome='chr7',
        reference_start=1000,
        max_mismatches=2,
    )[guide.id]
    assert profile.intended_target_found is False
    assert profile.exact_duplicate_sites == 1
    assert profile.total_off_targets == 1
    assert any('intended target row was not identified' in note for note in profile.notes)


def test_crisprware_index_passes_5prime_pam_and_bin_width(tmp_path, monkeypatch):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    monkeypatch.setattr(crisprware.shutil, "which", lambda executable: "/usr/bin/crisprware")
    seen = {}

    class Result:
        returncode = 0
        stderr = ""
        stdout = ""

    def fake_run(cmd, text, capture_output):
        seen["cmd"] = cmd
        return Result()

    monkeypatch.setattr(crisprware.subprocess, "run", fake_run)
    crisprware.build_crisprware_index(
        fasta,
        pam="TTTV",
        spacer_length=23,
        pam_5_prime=True,
        bin_width=14,
        output_directory=tmp_path / "out",
    )
    assert "--pam_5_prime" in seen["cmd"]
    assert "--bin_width" in seen["cmd"]
    assert seen["cmd"][seen["cmd"].index("--bin_width") + 1] == "14"
