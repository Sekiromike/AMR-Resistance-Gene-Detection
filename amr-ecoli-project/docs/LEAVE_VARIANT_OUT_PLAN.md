# Leave-One-Variant-Out Experiment

Status: ADOPTED 2026-09-29 by the decision owner, before any model in it was run
Motivation: [`POSTHOC_FAMILY_ANALYSIS.md`](POSTHOC_FAMILY_ANALYSIS.md). Public data
hold no adequate independent set (21 candidate isolates carry rare or absent
cephalosporinase alleles; see the project ledger), so the unseen-variant
hypothesis is tested by controlled simulation on development data only. The
external set is not reused.

## Question

When every carrier of a resistance allele is removed from training, do models
still recognise those isolates as resistant? Compared: allele-level catalogue,
curated family-level catalogue, and Evo 2.

## Data and targets

Human-clinical development isolates with a ceftriaxone reference (2,822).
Targets, by rule: every cephalosporinase allele (AMRFinderPlus Subtype AMR,
Subclass containing CEPHALOSPORIN, intrinsic blaEC excluded) carried by ≥ 30
of these isolates. Counts are from genotypes only.

| Target | Carriers | Role |
|---|---:|---|
| blaCTX-M-15 | 410 | primary |
| blaCTX-M-27 | 170 | primary |
| blaCTX-M-14 | 47 | primary |
| blaCTX-M-55 | 39 | primary |
| blaCMY-2 | 78 | control: no other CMY allele with ≥ 5 carriers, so gene-family pooling cannot help; only the drug-class indicator could |
| blaOXA-1 | 218 | descriptive: not itself a ceftriaxone-resistance mechanism; usually co-carried with blaCTX-M-15 |

## Procedure, per target allele A

- **Test:** every isolate carrying A.
- **Training:** every other isolate. Settings are chosen on the five frozen
  development folds restricted to training rows, then the model is refitted
  on all training rows, exactly as `predict_external.py` does. A never
  appears in training, so it never enters the allele vocabulary.
- **Models:** ridge_censored and xgboost_aft on four feature sets:
  allele-level AMRFinderPlus; family-level (allele plus gene family plus
  drug class, `build_hierarchical_features.py`); Evo 2 alone; AMRFinderPlus
  plus Evo 2.

## Read-out, stated in advance

Primary: among test isolates with a right-censored (above-panel) ceftriaxone
reference, the share predicted on the resistant side (> 2 mg/L), per target
and model, pooled over the four CTX-M targets. Paired differences use a
whole-lineage bootstrap over the pooled test isolates.

- **H1**, curated grouping generalises at least as well as Evo 2: pooled
  family minus Evo 2-alone, with the 95% lower bound above −0.05
  (non-inferiority margin of 5 points), for both learners.
- **H2**, grouping helps over alleles: pooled family minus allele, with the
  95% lower bound above 0, for both learners.
- **Control:** for blaCMY-2, a family-level gain would come from the
  drug-class indicator, not from gene-family pooling. It is reported as such.

Secondary: essential agreement on all test isolates per target and model.

## Scope

A simulation of unseen alleles inside one population, not a test of other
regions. It strengthens or weakens the post-hoc external finding; it does not
replace an independent external validation, which current public data cannot
provide.

## Results (2026-09-30)

Array `65051616` (6 targets, zero RuntimeWarnings); read-out `65054505`
(checksums verified). The first read-out attempt, `65051617`, stopped on a
software defect before computing anything: isolates carrying two target
alleles appear once per target, and pooled scoring refused the duplicate.
Fixed by scoring per target, with a regression test. Predictions unchanged.

**Pre-stated hypotheses, pooled over the four CTX-M targets (568 held-out
carriers above the panel):**

| Learner | H1 family − Evo 2 alone | H2 family − allele |
|---|---|---|
| ridge | **+0.046 [+0.023, +0.121]**: supported | **+0.915 [+0.836, +0.943]**: supported |
| XGBoost | **+0.111 [+0.074, +0.221]**: supported | **+0.963 [+0.916, +0.981]**: supported |

Family-level grouping is not only non-inferior to Evo 2 but better, and it
restores detection that allele-level features lose almost completely.

Share of held-out carriers above the panel placed on the resistant side, per
target (ridge / XGBoost):

| Target (above panel / tested) | allele | family | Evo 2 alone | AMR + Evo 2 |
|---|---|---|---|---|
| blaCTX-M-15 (361/410) | 0.080 / 0.042 | 0.972 / 1.000 | 0.922 / 0.859 | 0.235 / 0.687 |
| blaCTX-M-27 (134/170) | 0.022 / 0.022 | 1.000 / 1.000 | 0.948 / 0.910 | 0.575 / 0.791 |
| blaCTX-M-14 (37/47) | 0.054 / 0.027 | 1.000 / 1.000 | 0.973 / 1.000 | 0.189 / 0.676 |
| blaCTX-M-55 (36/39) | 0.111 / 0.056 | 1.000 / 1.000 | 1.000 / 1.000 | 0.194 / 0.861 |
| blaCMY-2 control (22/78) | 0.182 / 0.182 | 0.773 / 0.727 | 0.227 / 0.591 | 0.182 / 0.318 |
| blaOXA-1 descriptive (174/218) | 0.994 for every model (co-carried blaCTX-M-15 remains in training) |

Essential agreement on all held-out CTX-M carriers: family 0.74–0.95, Evo 2
alone 0.44–0.87, allele 0.00–0.08.

**Reading.**
1. Evo 2 embeddings on their own do generalise to unseen CTX-M alleles
   (0.86–1.00 detection), which confirms the direction of the exploratory
   external finding under controlled conditions.
2. Curated gene-family and drug-class annotation generalises further
   (0.97–1.00) and is more accurate at the MIC level. In the CMY-2 control,
   where no family relative exists in training, only the curated drug-class
   indicator recovers most carriers (0.73–0.77, against Evo 2 0.23–0.59).
3. Combining allele-level catalogue features with Evo 2 features loses most
   of Evo 2's generalisation (AMR + Evo 2 0.19–0.86): the combined model
   relies on allele indicators and predicts susceptible when the allele is
   absent. How foundation-model features are combined with curated features
   decides whether their generalisation survives.

Scope: a simulation within one population. It corroborates the post-hoc
external finding; it does not replace an independent external validation.
