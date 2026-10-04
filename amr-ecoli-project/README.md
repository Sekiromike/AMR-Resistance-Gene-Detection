# E. coli MIC prediction: curated resistance features versus Evo 2

This project compares frozen Evo 2 resistance-locus embeddings with curated
AMRFinderPlus features for predicting ceftriaxone, ciprofloxacin, and gentamicin
minimum inhibitory concentrations (MICs) in *Escherichia coli*.

**The experiments are complete.** Curated allele, gene-family, and drug-class
features outperformed the tested Evo 2 representation when CTX-M variants were
withheld from supervised training. Adding Evo 2 did not demonstrate improvement
in the original locked external evaluation.

## Paper and release

- [Manuscript PDF](paper/manuscript.pdf) and [editable source](paper/manuscript.md)
- [Supplementary tables](paper/supplement.pdf)
- [Versioned v1.0.0 release](https://github.com/Sekiromike/AMR-Resistance-Gene-Detection/releases/tag/v1.0.0)
- [Reproduction and artifact instructions](paper/REPRODUCIBILITY.md)
- [Verified results](paper/verification/audit.json)
- [Current status](docs/PROJECT_STATUS.md) and [scope closeout](docs/SCOPE_CLOSEOUT.md)
- [Resume description](docs/RESUME_DESCRIPTION.md)

Author: **Rushath Rajeev, independent researcher**. Computational resources:
UMass Unity HPC. Claude Code and Codex are acknowledged as AI assistance in the
paper. The submission package targets arXiv biology (`q-bio.GN`); it is not a
claim of submission, acceptance, or peer review.

## Study and outputs

The final cohort contains 29,619 isolate–antibiotic MIC rows from 11,520 isolates:
10,015 development isolates, including 3,849 human-clinical isolates, and 1,505
external Japanese JARBS isolates. The learners are interval-censored ridge
regression and XGBoost AFT. Nested development folds group lineages, genomic
neighbors, and identifier-linked duplicates. PCA and feature selection are
fitted within training folds.

The release provides cohort and exclusion tables, frozen splits, allele and
family features, sequence units, Evo 2 embeddings, frozen predictions, evaluation
and environment manifests, full comparison tables, manuscript sources, and
checksums. All 54 archived evaluations and 438 paired comparisons reproduce from
the saved predictions, and all 18 committed prediction locks match.

## Interpretation

The original external comparison was locked before reference attachment. The
family-feature follow-up was designed after external errors were inspected and
is exploratory. The subsequent variant-withholding experiment uses development
data and has not been replicated in a new external collection. The family model
includes allele, family, and drug-class indicators. The embedding representation
uses annotation-selected loci; it is not a full-genome or fine-tuned Evo 2 model.

EA measures agreement with an exact or censored MIC reference. Neither EA nor
the descriptive ceftriaxone >2 mg/L read-out is clinical classification accuracy.
Missing assay metadata prevent validation of the secondary S/I/R endpoint.
There is no clinically validated predictor or demonstrated new resistance
mechanism in this repository.

## Reproduce scores and paper

From this `amr-ecoli-project` directory, extract the release result assets as
specified in [the reproduction guide](paper/REPRODUCIBILITY.md), then run:

```bash
uv run --no-project --with-requirements paper/verify-requirements.txt \
  python paper/scripts/verify_results.py
uv run --no-project --with-requirements paper/render-requirements.txt \
  python paper/scripts/build_paper.py --tectonic /path/to/tectonic
python3 paper/scripts/build_paper.py --verify-only
python3 skills/manage-amr-research/scripts/status.py check
```

Scoring verification does not fit models or select on external outcomes.
Slurm execution scripts and environment definitions remain in `hpc/unity/` and
`workflow/envs/`. See [Unity instructions](docs/UNITY_HPC.md) for cluster setup;
keep user, allocation, and credential details outside Git.

Original code is licensed under [MIT](LICENSE). Third-party data, databases,
models, and software retain their own terms; see [third-party notices](THIRD_PARTY_NOTICES.md).

The current study is on the [`main` branch](https://github.com/Sekiromike/AMR-Resistance-Gene-Detection/tree/main/amr-ecoli-project). Its analysis plans and prediction-lock chronology are recorded in [study provenance](paper/provenance/study_chronology.json).
