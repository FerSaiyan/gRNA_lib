# Optional scoring backends

gRNA Library keeps external/scientific scoring engines behind explicit adapters. The core package still works without them, and every returned score includes enough provenance to tell which backend produced it.

## Check installed capabilities

```bash
grna-lib backend status
```

The same information is available to agents through the MCP tool `get_backend_status`.

## Rule Set 3

Rule Set 3 is the preferred optional SpCas9 sequence activity scorer.

```bash
python -m pip install -e '.[rs3]'
```

When a 30-nt SpCas9 context is available, guide results include `rs3_sequence`. Otherwise the package keeps the transparent sequence-quality fallback and explains why RS3 was unavailable.

## Indexed whole-genome specificity through CRISPRware

The local mismatch scanner is useful for short supplied sequences, but it is not a whole-genome CFD/GuideScan2 replacement.

For indexed specificity, gRNA Library can call an **externally installed** CRISPRware command. CRISPRware can score against crispr-ots or GuideScan2 indices. No CRISPRware source, models, binaries, or genome indices are bundled into this project.

Important licensing note: the CRISPRware repository currently states a UC Santa Cruz Noncommercial License and directs commercial users to UCSC Innovation Transfer. Review those terms before commercial use.

Example index build:

```bash
grna-lib genome index \
  --fasta hg38.fa \
  --pam NGG \
  --spacer-length 20 \
  --output-directory indexes
```

Example indexed ranking:

```bash
grna-lib guide rank \
  --sequence ACGT... \
  --crisprware-index indexes/hg38_crisprots/hg38_crisprots \
  --chromosome chr7 \
  --reference-start 55019016
```

`--reference-start` is the **0-based genomic coordinate corresponding to base 0 of the supplied sequence**. It and `--chromosome` matter: an indexed off-target engine needs the real intended locus so that the on-target perfect match is not treated as an off-target.

Optional controls:

```text
--threads
--mismatches
--rna-bulges
--dna-bulges
--crisprware-executable
```

The indexed result is stored separately from the local proxy and is preferred for ranking when available.

### Whole-genome mismatch profiles

For a per-guide view of the genomic matches behind the aggregate score, use the `crispr-ots` profiler directly:

```bash
grna-lib guide profile \
  --sequence ACGT... \
  --crispr-ots-index indexes/hg38_crisprots/hg38_crisprots \
  --chromosome chr7 \
  --reference-start 55019016 \
  --mismatches 4
```

The profile reports, for every candidate guide:

- the indexed specificity value returned by the search;
- off-target counts for each mismatch distance from 0 through the requested maximum;
- the number of additional exact genomic copies after excluding the coordinate-matched intended target;
- whether the intended target was found at the supplied chromosome/start/strand;
- a bounded list of representative genomic loci ordered by mismatch count.

`--hit-limit` limits only the loci retained in the JSON response. It does not truncate the mismatch counts.

The profiler uses the current GuideScan2-compatible per-hit CSV emitted by `crispr-ots enumerate`. The intended target is excluded only when its chromosome, 0-based genomic start, strand and zero-mismatch status agree with the coordinates supplied to gRNA Library. If that exact row cannot be identified, zero-mismatch hits are retained rather than guessed away.

Current `crispr-ots` bin scanning supports mismatch enumeration but rejects non-zero RNA/DNA bulges. Use the CRISPRware/GuideScan2 scoring path when bulge-aware specificity is required.

## DeepPrime through GenET

Learned prime-edit efficiency is available through the MIT-licensed GenET package maintained by the DeepPrime authors.

Because current GenET dependency metadata requires TensorFlow <2.10, the supported installation path for this optional backend is a Python 3.10 environment:

```bash
python -m pip install -e '.[deepprime]'
```

Then:

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

DeepPrime requires sufficient context around the edit and supports its own defined edit/PAM/model ranges. If the edit cannot be represented for GenET, gRNA Library raises an explicit error rather than silently substituting the structural score.

When DeepPrime is active, each candidate keeps both:

- `structural_prime`: the transparent built-in structural score.
- `deepprime`: learned score plus PE system, cell type, backend name, and installed GenET version.

Ranking prefers a valid learned DeepPrime score when one is present.

## Why these are adapters

Keeping these engines separate avoids three problems:

1. heavyweight ML/index dependencies becoming mandatory for the Tkinter/basic CLI install;
2. silently mixing scores trained for different editors or cell contexts;
3. redistributing third-party code under incompatible or more restrictive licenses.

The JSON schema stays stable whether a backend is installed or not, which is particularly useful for MCP/agent workflows.
