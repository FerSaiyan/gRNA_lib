# gRNA Library

A local-first CRISPR guide-RNA and prime-edit design toolkit with a shared Python core and four frontends: the original-style CustomTkinter UI, a localhost web UI, JSON-friendly CLI commands, and an MCP server for AI agents.

## Repository cleanup

Historically, `main` contained ZIP releases while the inspectable TCC-era Python source lived on `master`. This branch is based on `main` and brings the source back into the default-branch lineage. The original `master` files are preserved under `legacy/master_snapshot/`; the top-level legacy filenames remain as compatibility shims.

After this branch is merged into `main`, `main` should be treated as the canonical development branch. The old `master` branch can remain as historical provenance or be archived/deleted later after confirming no external automation depends on it.

## What changed

- Correct IUPAC PAM matching instead of assuming only the first PAM base can be `N`.
- Candidate enumeration on both strands with coordinates normalized to the supplied reference sequence.
- U6 compatibility now prefixes an extra 5' G when needed; it never mutates the genomic spacer by appending G.
- `TTTT` is correctly flagged as a Pol III termination risk.
- GC and sequence warnings are annotations rather than a single opaque hand-written score.
- Optional Rule Set 3 sequence scoring (`pip install .[rs3]`) with a transparent fallback when RS3 is unavailable.
- Local off-target scanning counts 0–4 mismatch PAM-compatible sites and reports a clearly labeled specificity **proxy**. It does not claim to be CFD/GuideScan2.
- Prime-edit enumeration produces PBS/RTT combinations, pegRNA extension sequences, PAM-disruption annotation, and PE3/PE3b-style nicking-guide candidates.
- Tkinter, localhost web, CLI, and MCP all call the same core functions.

## Install

```bash
python -m pip install -e .
python -m pip install -e '.[ui,web,mcp,rs3]'
```

## Tkinter

```bash
grna-lib-tk
# or: python GUI_gRNA.py
```

The default guide-design flow intentionally retains the original visual hierarchy: 1200×1200 window, large `Guide RNA Library` title, large Helvetica controls, stacked input/output frames, and full-width Run/Rank buttons.

## Localhost web UI

```bash
uvicorn grnalib.web:app --reload
```

Open `http://127.0.0.1:8000`. The HTML/CSS mirrors the Tkinter hierarchy instead of redesigning it around a generic dashboard.

## CLI

All main commands emit JSON so agents do not have to scrape human-oriented terminal tables.

```bash
grna-lib guide design --sequence ACGT...
grna-lib guide rank --sequence ACGT... --genome genome.fa
grna-lib prime design --sequence ACGT... --position 42 --ref G --alt A
```

## MCP

```bash
grna-lib-mcp
```

Tools exposed:

- `design_grnas`
- `rank_grnas`
- `design_prime_candidates`
- `explain_candidate`

## Scoring policy

Ranking is deliberately hierarchical rather than mixing unrelated scores into a falsely precise weighted sum:

1. hard sequence filters;
2. specificity when a genome/background sequence is available;
3. RS3 sequence score when installed, otherwise a transparent sequence-quality fallback.

Prime-edit candidates currently use a transparent structural fallback that prefers PAM disruption and PE3b-compatible nicking guides. The API was designed so PRIDICT2/DeepPrime adapters can be added without changing the frontends or output schema.

## Important limitations in this first modernization patch

- The built-in off-target implementation is intended for local/small reference sequences. Whole-genome production workflows should use an indexed engine such as GuideScan2/CRISPRware and feed its specificity output through a future scorer adapter.
- The built-in prime score is **not** PRIDICT2 or DeepPrime. The candidate-generation data model is ready for those learned scorers, but they are intentionally not vendored or silently approximated.
- Prime extension generation currently targets Cas9-family 3'-PAM systems.

## Tests

```bash
python -m pip install -e '.[dev]'
pytest
```
