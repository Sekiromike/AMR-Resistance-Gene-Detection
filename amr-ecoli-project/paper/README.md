# Research paper and reproducibility release

**Curated resistance gene families outperform Evo 2 embeddings on unseen
variants in Escherichia coli** — Rushath Rajeev, independent researcher.

Prepared for arXiv biology (`q-bio.GN`; possible `q-bio.QM` cross-list).
No arXiv submission or acceptance is claimed.

## Read

- [Manuscript PDF](manuscript.pdf), [Markdown](manuscript.md), [LaTeX](manuscript.tex)
- [Supplement PDF](supplement.pdf) and [complete CSV tables](tables/)
- [Verification report](verification/audit.json) and [claim map](provenance/CLAIM_MAP.md)
- [Reproduction instructions](REPRODUCIBILITY.md)
- [Submission checklist and metadata](SUBMISSION_CHECKLIST.md)
- [Reporting limitations](EDITORIAL_NOTES.md)

The numerical audit reproduced 54 original evaluations, 438 paired comparisons,
192 source-exclusion summaries, the external variant subgroups, and the
withholding read-out. All 18 committed prediction locks match. The final paper
uses full-precision values from the verified machine-readable outputs.

`manuscript.pdf` and `supplement.pdf` are compiled from their actual LaTeX
sources using Tectonic 0.17.0. Vector PDF figures and resolved references are
included in the source upload archive. The HTML file is a convenience preview.

## Rebuild

```bash
uv run --no-project --with-requirements paper/render-requirements.txt \
  python paper/scripts/build_paper.py --tectonic /path/to/tectonic
python3 paper/scripts/build_paper.py --verify-only
```

Run these commands from the project directory. The build uses committed
verification exports; raw artifacts are needed only to repeat the numerical
audit, as described in the reproduction guide. `provenance/artifact.sha256`
checks the paper files and excludes bulky local artifacts and release bundles.

The current study is on the [`main` branch](https://github.com/Sekiromike/AMR-Resistance-Gene-Detection/tree/main/amr-ecoli-project). Original study-plan and prediction-lock records are retained in [study_chronology.json](provenance/study_chronology.json).
