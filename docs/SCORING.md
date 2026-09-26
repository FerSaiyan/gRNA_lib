# Scoring and ranking

The project intentionally keeps **candidate generation**, **component scoring**, and **ranking policy** separate.

## Default guide policy

1. Hard sequence filters (currently `TTTT`/Pol III termination risk).
2. Indexed specificity when an external index backend is configured.
3. Otherwise, the local specificity proxy when a background sequence is supplied.
4. Rule Set 3 sequence score when the optional `rs3` package is installed and a 30-nt SpCas9 context is available.
5. Otherwise, a transparent sequence-quality fallback based primarily on GC balance.

The built-in specificity value is named `specificity_proxy` on purpose. It counts PAM-compatible sites with up to four spacer mismatches and applies explicit mismatch weights. It is **not** labeled CFD and is not a substitute for indexed whole-genome scoring.

An indexed backend result is stored as `indexed_specificity` with its source. If an indexed backend was requested but omitted a candidate, that candidate is treated as missing indexed evidence rather than receiving neutral/perfect specificity.

## Prime editing

Prime candidates expose spacer, PBS, RTT, full extension, nick-to-edit distance, PAM-disruption status, and PE3/PE3b-style nicking-guide candidates.

The built-in `structural_prime` remains a transparent fallback. The optional GenET adapter adds a separate `deepprime` object containing:

- learned score;
- PE system;
- cell type;
- backend name;
- GenET package version.

The structural score is never overwritten. When a valid learned score is present, prime candidate ranking uses it before the structural fallback.

See `BACKENDS.md` for installation, model-context, whole-genome coordinate, and licensing details.
