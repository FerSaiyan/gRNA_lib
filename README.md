# gRNA Library

Design and rank CRISPR guide RNAs and prime-editing candidates from a local Python tool.

gRNA Library started as an undergraduate project for finding SpCas9 guides. The current version keeps the original desktop interface but moves the design and scoring code into a shared Python package.

You can use it from the Tkinter app, the command line, a localhost web interface, or an MCP client.

[![Tests](https://github.com/FerSaiyan/gRNA_lib/actions/workflows/tests.yml/badge.svg)](https://github.com/FerSaiyan/gRNA_lib/actions/workflows/tests.yml)

## What it does

### CRISPR guide design

- Finds candidate guides on both DNA strands.
- Supports IUPAC PAM sequences and custom PAMs.
- Keeps genomic spacer sequences separate from expression changes such as an added 5′ G for U6.
- Flags poly-T and unusual GC content.
- Supports SpCas9, SpCas9-NG, SaCas9, AsCas12a, and custom nuclease settings.

### Guide ranking

For SpCas9, Rule Set 3 can be used for on-target activity scoring.

For small supplied reference sequences, the built-in scanner reports PAM-compatible off-targets with up to four mismatches. This score is deliberately reported as a local specificity proxy rather than CFD.

For genome-wide searches, gRNA Library can use an external CRISPRware installation with a crispr-ots or GuideScan2 index. A separate `guide profile` command reports mismatch counts and representative genomic off-target loci from a crispr-ots index.

See [Scoring and ranking](docs/SCORING.md) for the ranking rules.

### Prime editing

The prime-editing designer enumerates:

- pegRNA spacers
- PBS lengths
- RTT lengths
- pegRNA extensions
- nick-to-edit distances
- PAM-disrupting designs
- PE3 and PE3b nicking-guide candidates

A simple structural score is available without extra dependencies.

DeepPrime scoring is also supported through the GenET package when installed in a compatible Python environment.

### Interfaces

The same design code is used by:

- CustomTkinter desktop app
- command-line interface
- localhost FastAPI interface
- MCP server

This keeps guide generation and scoring consistent across interfaces.

## Install

Clone the repository:

```bash
git clone https://github.com/FerSaiyan/gRNA_lib.git
cd gRNA_lib
```

Install the core package:

```bash
python -m pip install -e .
```

For the desktop UI, web interface, MCP server, and Rule Set 3:

```bash
python -m pip install -e '.[ui,web,mcp,rs3]'
```

Python 3.10 or newer is required.

## Desktop app

Run:

```bash
grna-lib-tk
```

The old entry point still works:

```bash
python GUI_gRNA.py
```

The current UI stays close to the original TCC-era application rather than replacing it with a separate desktop design.

## Command line

Find guides:

```bash
grna-lib guide design \
  --sequence AAAAGCGCGCGCGCGCGCGCGCGCTGGAAAA
```

Rank guides against a supplied background sequence:

```bash
grna-lib guide rank \
  --sequence AAAAGCGCGCGCGCGCGCGCGCGCTGGAAAA \
  --genome genome.fa
```

Check which optional scoring backends are installed:

```bash
grna-lib backend status
```

### Genome-wide off-target scoring

If CRISPRware and a compatible genome index are installed:

```bash
grna-lib guide rank \
  --sequence ACGT... \
  --crisprware-index /path/to/index \
  --chromosome chr7 \
  --reference-start 55019016
```

`--reference-start` is the zero-based genomic coordinate corresponding to the first base of the supplied sequence.

To inspect the mismatch distribution and genomic loci behind the score:

```bash
grna-lib guide profile \
  --sequence ACGT... \
  --crispr-ots-index /path/to/index \
  --chromosome chr7 \
  --reference-start 55019016 \
  --mismatches 4
```

Setup, output fields, and licensing details are in [Optional scoring backends](docs/BACKENDS.md).

## Prime editing

Basic candidate generation:

```bash
grna-lib prime design \
  --sequence ACGT... \
  --position 72 \
  --ref G \
  --alt A
```

DeepPrime scoring:

```bash
grna-lib prime design \
  --sequence ACGT... \
  --position 72 \
  --ref G \
  --alt A \
  --scorer deepprime \
  --pe-system PE2max \
  --cell-type HEK293T
```

DeepPrime is provided through GenET. Its current dependency stack is best installed in a separate Python 3.10 environment:

```bash
python -m pip install -e '.[deepprime]'
```

## Web interface

Install the web dependencies and start the local server:

```bash
python -m pip install -e '.[web]'
uvicorn grnalib.web:app --reload
```

Then open:

```text
http://127.0.0.1:8000
```

## MCP

Install the MCP dependency:

```bash
python -m pip install -e '.[mcp]'
```

Start the server:

```bash
grna-lib-mcp
```

Available tools include guide design and ranking, prime-edit design, backend checks, off-target index construction, and candidate explanations.

## Scoring notes

Guide ranking follows this order:

1. sequence filters
2. indexed specificity, when available
3. local specificity proxy, when a background sequence is supplied
4. Rule Set 3 for SpCas9, when installed
5. built-in sequence score as a fallback

Prime-edit candidates keep the built-in structural score even when DeepPrime is used. DeepPrime results also record the PE system, cell type, backend, and package version used to produce the prediction.

More detail is available in [docs/SCORING.md](docs/SCORING.md).

## Current limitations

The built-in off-target scanner is meant for short sequences, not whole genomes. Use an indexed backend for genome-wide specificity analysis.

Prime-edit extension generation currently targets Cas9-family systems with a 3′ PAM.

PRIDICT2 is not integrated yet.

CRISPRware is an external dependency and is not distributed with this repository. Its upstream license currently restricts commercial use; see [docs/BACKENDS.md](docs/BACKENDS.md).

## Tests

```bash
python -m pip install -e '.[dev]'
pytest
```

The test suite currently runs on Python 3.11 and 3.12 in GitHub Actions.

## Project history

The original TCC-era implementation is preserved under [`legacy/master_snapshot`](legacy/master_snapshot/).

The current implementation lives in `src/grnalib`.
