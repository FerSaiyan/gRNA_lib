# Scoring and ranking

The project intentionally keeps **candidate generation**, **component scoring**, and **ranking policy** separate.

## Default guide policy

1. Hard sequence filters (currently `TTTT`/Pol III termination risk).
2. Specificity when a genome/background sequence is supplied.
3. Rule Set 3 sequence score when the optional `rs3` package is installed and a 30-nt SpCas9 context is available.
4. Otherwise, a transparent sequence-quality fallback based primarily on GC balance.

The built-in specificity value is named `specificity_proxy` on purpose. It counts PAM-compatible sites with up to four spacer mismatches and applies explicit mismatch weights. It is **not** labeled CFD and is not a substitute for indexed whole-genome GuideScan2/CRISPRware scoring.

## Prime editing

Prime candidates expose spacer, PBS, RTT, full extension, nick-to-edit distance, PAM-disruption status, and PE3/PE3b-style nicking-guide candidates. The default `structural_prime` score is a transparent fallback, not PRIDICT2 or DeepPrime.

A future learned scorer should populate a new named score field and provenance metadata rather than overwrite `structural_prime`; all frontends already pass the `scores` object through unchanged.
