# Paper draft

This folder contains a working manuscript draft for the modernized gRNA Library project.

## Files

- `main.tex` — manuscript draft
- `references.bib` — literature and software references

Compile from this folder with:

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

## Framing

The draft deliberately does **not** claim that arbitrary genomes or custom PAMs are unique to gRNA Library. Current tools such as crisprVerse, GuideScan2, CHOPCHOP/CRISPOR local deployments, Cas-OFFinder, CRISPRitz, and CRISPRware already cover substantial parts of that space.

The proposed contribution is the combination of:

- local-first operation;
- one Python core across desktop, web, CLI, API, and MCP;
- configurable IUPAC PAM / spacer / PAM-side enumeration;
- conventional guide and prime-editing design in one package;
- optional scoring/index backends with provenance;
- whole-genome mismatch profiling;
- deterministic structured tools that AI agents can call directly.

## Before submission

1. Add authors, affiliations, contributions, acknowledgements, and funding.
2. Add a top-level software license. The repository currently has no `LICENSE` file, so the manuscript uses “public source” rather than “open source”.
3. Run the scientific benchmarks proposed in the Evaluation plan before turning them into Results.
4. Add figures, screenshots, quantitative tables, and runtime measurements.
5. Re-check the comparison table immediately before submission because web tools change.
