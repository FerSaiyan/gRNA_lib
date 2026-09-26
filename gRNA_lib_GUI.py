"""Compatibility API for notebooks/scripts written against the TCC-era module."""
from grnalib import design_guides, rank_guides, resolve_nuclease, reverse_complement


def complimentary_sequence(gene_sequence):
    return reverse_complement(gene_sequence)


def find_gRNA(gene_sequence, PAM='NGG'):
    spec = resolve_nuclease('SpCas9', pam=PAM)
    return [[g.spacer, g.start, g.pam] for g in design_guides(gene_sequence, spec) if g.strand == '+']


def gRNA_ranking(gRNA_list, genome_sequence, Sp_cas9=True):
    # Preserve the old return shape while using actual candidates from the supplied genome.
    spec = resolve_nuclease('SpCas9')
    wanted = {x[0] if isinstance(x, (list, tuple)) else str(x) for x in gRNA_list}
    guides = [g for g in design_guides(genome_sequence, spec) if g.spacer in wanted]
    ranked = rank_guides(guides, spec=spec, genome_sequence=genome_sequence or None)
    return [(g.spacer, g.scores) for g in ranked]


def main(gene_sequence_plus, PAM='NGG', show_sequence_minus=False):
    spec = resolve_nuclease('SpCas9', pam=PAM)
    guides = design_guides(gene_sequence_plus, spec)
    plus = [[g.spacer, g.start, g.pam] for g in guides if g.strand == '+']
    minus = [[g.spacer, g.start, g.pam] for g in guides if g.strand == '-']
    return plus, minus
