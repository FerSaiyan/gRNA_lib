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
- Optional Rule Set 3 sequence scoring with a transparent fallback when RS3 is unavailable.
- Local off-target scanning counts 0–4 mismatch PAM-compatible sites and reports a clearly labeled specificity **proxy**.
- Optional indexed whole-genome specificity through an externally installed CRISPRware/GuideScan2 or crispr-ots index.
- Prime-edit enumeration produces PBS/RTT combinations, pegRNA extension sequences, PAM-disruption annotation, and PE3/PE3b-style nicking-guide candidates.
- Optional learned DeepPrime efficiency scoring through the MIT-licensed GenET package.
- Tkinter, localhost web, CLI, and MCP all call the same core functions.

## Install

Basic package:

```bash
python -m pip install -e .
```

Common frontends and Rule Set 3:

```bash
python -m pip install -e '.[ui,web,mcp,rs3]'
```

DeepPrime/GenET currently fits best in a separate Python 3.10 environment:

```bash
python -m pip install -e '.[deepprime]'
```

CRISPRware is intentionally not a Python dependency of this project. Install it separately if you want indexed whole-genome specificity; see [docs/BACKENDS.md](docs/BACKENDS.md) for backend and licensing notes.

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

Open `http://127.0.0.1:8000`. The HTML/CSS mirrors the Tkinter hierarchy instead of redesigning it around a generic dashboard. Advanced indexed/learned scorer parameters are also available through the API, CLI and MCP without forcing them into the simple default UI.

## CLI

Commands emit JSON so agents do not have to scrape human-oriented terminal tables.

```bash
grna-lib backend status

grna-lib guide design --sequence ACGT...
grna-lib guide rank --sequence ACGT... --genome genome.fa

grna-lib guide rank \
  --sequence ACGT... \
  --crisprware-index /path/to/index \
  --chromosome chr7 \
  --reference-start 55019016

grna-lib prime design \
  --sequence ACGT... \
  --position 72 --ref G --alt A \
  --scorer deepprime --pe-system PE2max --cell-type HEK293T
```

## MCP

```bash
grna-lib-mcp
```

Tools exposed:

- `get_backend_status`
- `design_grnas`
- `rank_grnas`
- `build_offtarget_index`
- `design_prime_candidates`
- `explain_candidate`

## Scoring policy

Ranking is deliberately hierarchical rather than mixing unrelated scores into a falsely precise weighted sum:

1. hard sequence filters;
2. indexed specificity when configured, otherwise the local specificity proxy when supplied;
3. RS3 sequence score when installed, otherwise a transparent sequence-quality fallback.

If an indexed run was requested and no indexed score is returned for a guide, that missing value is not treated as perfect specificity.

Prime-edit candidates preserve the transparent `structural_prime` score. When GenET/DeepPrime is requested and returns a learned score, ranking prefers the learned score and keeps PE-system/cell-type/version provenance.

See [docs/SCORING.md](docs/SCORING.md) and [docs/BACKENDS.md](docs/BACKENDS.md).

## Important limitations

- The built-in off-target implementation is intended for local/small reference sequences. Use an indexed backend for production whole-genome specificity.
- CRISPRware is an optional external integration and is **not vendored**. Its current upstream license is noncommercial; review upstream terms for your intended use.
- DeepPrime is optional and its current GenET dependency stack is best isolated in Python 3.10.
- Prime extension generation currently targets Cas9-family 3'-PAM systems.
- PRIDICT2 is not yet integrated; no score is labeled as PRIDICT2 unless a genuine model adapter is added.

## Tests

```bash
python -m pip install -e '.[dev]'
pytest
```
