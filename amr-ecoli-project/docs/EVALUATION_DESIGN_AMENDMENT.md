# Evaluation Design Amendment: Nested Cross-Validation and Precision Planning

Amendment: 007
Date: 2026-09-27
Status: ADOPTED 2026-09-27 by the decision owner, with additions 4 and 5
Decision owner: Rushath Rajeev
Depends on: [`EXTERNAL_COHORT_AMENDMENT.md`](EXTERNAL_COHORT_AMENDMENT.md) (002),
[`EXTERNAL_THRESHOLD_AMENDMENT.md`](EXTERNAL_THRESHOLD_AMENDMENT.md) (003),
[`GROUPING_VARIABLE_AMENDMENT.md`](GROUPING_VARIABLE_AMENDMENT.md) (004),
[`MIC_EVALUATION_SPEC.md`](MIC_EVALUATION_SPEC.md)

Everything here was measured on the validated cohort `64948707`
(`research_cohort.csv` SHA-256 `6b066cfd…`) before any genomic model was fitted.
Development MIC comparators were read; no external MIC value or comparator was
read. External figures use lineage-group sizes (genotype) only.

## What was measured

**Leakage units nest inside lineage groups.** Uniting SLV lineage groups,
99.9%-ANI genomic clusters and deduplication groups over development gives 485
units against 489 lineage groups. Genomic clusters (2,711; largest 717) sit
almost entirely inside SLV groups, unlike the external clustering that chained
(amendment 004).

**Effective sample size is small.** Human-clinical development has 218 lineage
groups but an inverse Simpson of 10.9; the SLV:ST131 group alone is 24.4%.
Resistance is strongly clonal. Taking the within-lineage intraclass correlation
(ICC) of the above-panel indicator as a planning value, with the design effect
DEFF = 1 + (m_w − 1)·ICC (m_w the size-weighted mean group size), the 95%
half-width for essential agreement at a planning EA of 0.90 is:

| Drug | Set | n | Groups | ICC | DEFF | n_eff | Half-width |
|---|---|---:|---:|---:|---:|---:|---:|
| ceftriaxone | dev human clinical | 2,822 | 188 | 0.276 | 66 | 43 | 0.090 |
| ciprofloxacin | dev human clinical | 3,809 | 218 | 0.496 | 175 | 22 | 0.126 |
| gentamicin | dev human clinical | 3,348 | 209 | 0.135 | 41 | 82 | 0.065 |
| ceftriaxone | dev all sources | 7,097 | 408 | 0.234 | 66 | 108 | 0.056 |
| ciprofloxacin | dev all sources | 8,493 | 437 | 0.446 | 161 | 53 | 0.081 |
| gentamicin | dev all sources | 9,514 | 484 | 0.141 | 49 | 194 | 0.042 |
| each drug | external | 1,505 | 106 | 0.10 (assumed) | 23 | 66 | 0.072 |

These are **conservative**: they assume model errors cluster by lineage as
strongly as the phenotype does. A model that captures lineage-associated
resistance leaves errors that cluster less, so the true precision will be
better, by an amount that cannot be known before a model exists.

## Decision (proposed)

### 1. Nested lineage-grouped cross-validation replaces the internal holdout

Amendment 002 planned an internal holdout carved from development as whole
genomic clusters. With about 11 effective lineages in the primary population,
a 20% holdout holds about two effective lineages, and whether it contains
SLV:ST131 would decide the result. Its estimate would be too imprecise to use,
and it would remove a fifth of an already small effective sample from training.

Instead:

- **Outer loop:** the five lineage-grouped folds report performance. Every
  outer test fold is sealed from all choices made for it.
- **Inner loop:** all model selection, feature selection, hyperparameters,
  calibration, thresholds and early stopping use grouped folds inside the
  outer training data only.
- **Untouched test:** the external-locked set remains the single evaluation
  that no development choice sees.

This gives an unbiased estimate after selection (Varma and Simon, 2006) using
every lineage.

### 2. The cross-validation unit also absorbs genomic clusters

The CV unit becomes the union of SLV lineage group, genomic cluster and
deduplication group. It costs 4 merges (489 to 485 units), and it closes the
remaining path by which 99.9%-ANI near-neighbours could sit in different
folds. Amendment 004 kept genomic clusters out of the CV unit because they
chained in the external set; in development they do not.

### 3. Precision is planned from development model errors before unsealing

The adequacy criterion is a 95% half-width of at most **0.05** (5 percentage
points) for essential agreement at a planning EA of 0.90.

The planning ICC is the within-lineage ICC of **essential-agreement errors**
of the first genomic model (the AMRFinderPlus model) in outer development
folds, computed before the external set is unsealed. The external half-width
is then computed per drug from that ICC and the external lineage-group sizes.
A drug that fails is reported as underpowered, not dropped or quietly
included. Reported intervals remain the whole-lineage bootstrap, which accounts
for clustering directly; this calculation only decides the underpowered label.

### 4. A nearest-relative baseline joins the no-genome baseline

Every model must also beat a genomic nearest-relative baseline on identical
rows: for each outer-test isolate and drug, predict from the training isolate
(outer training folds, same population, with that drug's MIC) of highest ANI
in the frozen development comparisons (ties: higher minimum aligned fraction,
then isolate ID). The prediction is that relative's reference set boundary:
its exact dilution, or the finite end of its censored range (`<= x` gives x,
`< x` gives x−1, `>= x` gives x, `> x` gives x+1). It has no tuning. A model
that cannot beat it has learned relatedness, not resistance.

### 5. The dominant lineage is reported separately

Every per-drug result is also reported for SLV:ST131 and for all other
lineages, each with its own whole-lineage bootstrap interval (the ST131 stratum
is one lineage group, so its interval resamples isolates within it and is
labelled as such). Ceftriaxone is reported as effectively dichotomous: 94% of
its human-clinical development references are censored.

## Consequences

- Any headline claim carries wide, honest intervals. Most AMR prediction
  studies report binomial intervals over isolates, which here would be
  roughly ten times too narrow.
- `build_split_manifest.py` gains the genomic-cluster union, and model code
  must accept an outer fold and run its own inner folds.
- The external-stress and lineage-novel strata of amendment 002 are unchanged.
