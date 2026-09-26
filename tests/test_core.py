from grnalib import PrimeEdit, design_guides, design_prime_edit, resolve_nuclease, reverse_complement
from grnalib.sequence import iupac_match


def test_reverse_complement():
    assert reverse_complement('ACGT') == 'ACGT'


def test_iupac_pam():
    assert iupac_match('AGG', 'NGG')
    assert iupac_match('TGA', 'NGR')
    assert not iupac_match('TGT', 'NGR')


def test_u6_g_is_prefixed_not_appended():
    seq = 'AAAA' + 'ACACACACACACACACACAC' + 'AGG' + 'AAAA'
    guides = design_guides(seq, resolve_nuclease())
    guide = next(g for g in guides if g.strand == '+' and g.spacer == 'ACACACACACACACACACAC')
    assert guide.expressed_spacer == 'GACACACACACACACACACAC'
    assert guide.spacer.endswith('C')


def test_poly_t_is_flagged():
    seq = 'AAAA' + 'GAAAATTTTCCCCGGGGAAA' + 'TGG' + 'AAAA'
    guides = design_guides(seq, resolve_nuclease())
    guide = next(g for g in guides if g.strand == '+')
    assert any('poly-T' in w for w in guide.warnings)


def test_minus_coordinates_are_reference_coordinates():
    # CCN on the reference creates an NGG site on the reverse complement.
    seq = 'CCATTTTTTTTTTTTTTTTTTT' + 'AAAAAA'
    guides = design_guides(seq, resolve_nuclease())
    minus = [g for g in guides if g.strand == '-']
    assert minus
    assert all(0 <= g.start < g.end <= len(seq) for g in minus)


def test_prime_design_makes_rtt_pbs_candidates_for_downstream_edit():
    # Protospacer at 4..24, PAM 24..27, nick at 21. Edit at 23 is downstream of nick.
    seq = 'AAAA' + 'GCGCGCGCGCGCGCGCGCGC' + 'TGG' + 'ACGTACGTACGTACGTACGTACGTACGTACGT'
    edit = PrimeEdit(position=23, ref=seq[23], alt='A' if seq[23] != 'A' else 'C')
    candidates = design_prime_edit(seq, edit, resolve_nuclease(), pbs_lengths=[10], rtt_lengths=[10, 15])
    assert candidates
    c = candidates[0]
    assert c.extension_sequence == c.rtt_sequence + c.pbs_sequence
    assert c.pbs_length == 10


def test_exact_intended_site_is_not_counted_as_offtarget():
    from grnalib.scoring import local_offtarget_summary
    seq = 'AAAA' + 'GCGCGCGCGCGCGCGCGCGC' + 'TGG' + 'AAAA'
    spec = resolve_nuclease()
    guide = next(g for g in design_guides(seq, spec) if g.strand == '+')
    summary = local_offtarget_summary(guide, seq, spec)
    assert summary.exact_other_sites == 0
