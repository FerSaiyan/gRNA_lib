from __future__ import annotations

import re

IUPAC = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'}, 'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'},
}

_RC = str.maketrans('ACGTRYKMSWBDHVN', 'TGCAYRMKSWVHDBN')


def normalize_dna(sequence: str) -> str:
    seq = re.sub(r'\s+', '', sequence.upper()).replace('U', 'T')
    invalid = sorted(set(seq) - set(IUPAC))
    if invalid:
        raise ValueError(f"Invalid DNA symbols: {', '.join(invalid)}")
    return seq


def reverse_complement(sequence: str) -> str:
    return normalize_dna(sequence).translate(_RC)[::-1]


def iupac_match(sequence: str, pattern: str) -> bool:
    sequence = normalize_dna(sequence)
    pattern = normalize_dna(pattern)
    if len(sequence) != len(pattern):
        return False
    return all(base in IUPAC[p] for base, p in zip(sequence, pattern))


def hamming_distance(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError('Hamming distance requires equal-length sequences')
    return sum(x != y for x, y in zip(a, b))
