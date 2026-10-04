# Post-hoc Analysis: Curated Gene Families versus Evo 2

Status: specified 2026-09-29 by the decision owner, before this analysis's
predictions were made. **Post hoc**: defined after the locked external
evaluation was scored, prompted by its exploratory finding
([`EXTERNAL_RESULTS.md`](EXTERNAL_RESULTS.md)). It cannot change the
pre-registered verdict and supports hypotheses only.

## Question

In the external set, 102 isolates with ceftriaxone > 64 mg/L carry
resistance-gene variants absent or nearly absent from development (blaCTX-M-2,
-8, -3, -24, -65, blaIMP-1). Allele-level catalogue models placed 13–15% of
them on the resistant side; Evo 2 alone placed 78–81%. Is that advantage just
grouping related variants, which curated knowledge also provides?

## Rival model, fixed before prediction

Features are the AMRFinderPlus allele-level symbols (unchanged) **plus**:

1. **Gene family.** For genes (no `_` in the symbol): beta-lactamases drop the
   final allele number (`blaCTX-M-15` → `family:blaCTX-M`, `blaTEM-1` →
   `family:blaTEM`, `blaIMP-1` → `family:blaIMP`); other genes drop trailing
   allele digits and a final lower-case allele letter (`qnrS1` →
   `family:qnrS`, `aac(3)-IId` → `family:aac(3)-II`). Point mutations
   (`gyrA_S83L`) are not pooled.
2. **Drug class.** One indicator per AMRFinderPlus Subclass present
   (`subclass:CEPHALOSPORIN`, `subclass:CARBAPENEM`, ...).

The same in-fold frequency filter (five training isolates), learners
(ridge_censored and xgboost_aft), grids and training (human-clinical
development, settings chosen on the five frozen folds) are used.

## Procedure

1. Development nested CV with the family features, for reference.
2. External predictions without reading external MICs; hashes committed to
   `config/external_posthoc_lock.tsv` before scoring.
3. Scoring on the same external set.

## Read-out, stated in advance

- Primary: the share of the 102 unseen-variant carriers placed on the
  resistant side (> 2 mg/L), set against 0.13–0.15 (allele-level) and
  0.78–0.81 (Evo 2 alone).
- Secondary: overall external essential agreement per drug, and the paired
  difference from the allele-level model.

Interpretation: if the family model reaches the Evo 2-alone share, the Evo 2
advantage on these isolates is attributable to grouping that curation already
supplies. If it stays well below, Evo 2 representations carry generalisation
beyond curated families, which a pre-registered test on a new independent set
would then need to confirm.

## Results (2026-09-29)

Predictions: job `65022960` (development nested CV and external predictions);
external predictions locked in `config/external_posthoc_lock.tsv` (commit
`ac3a270`) before scoring; scoring: job `65023163` (checksums verified).

**Primary read-out.** Of the 102 external isolates with ceftriaxone > 64 mg/L
carrying unseen variants, the share placed on the resistant side was:

| Model | Unseen-variant carriers | Common-variant carriers (546) |
|---|---:|---:|
| allele-level ridge / XGBoost | 0.147 / 0.127 | 0.998 / 1.000 |
| **family-level ridge / XGBoost** | **0.990 / 0.990** | 1.000 / 1.000 |
| Evo 2 alone, ridge / XGBoost | 0.784 / 0.814 | 0.907 / 0.963 |
| AMR + Evo 2, ridge / XGBoost | 0.196 / 0.461 | 0.996 / 0.982 |

By the rule stated in advance, the Evo 2 advantage on these isolates is
attributable to grouping that curated knowledge already supplies, and curated
grouping supplies it more completely.

**Secondary read-out**, external essential agreement (95% whole-lineage
interval) and paired difference from the allele-level model:

| Model | Ceftriaxone | Ciprofloxacin | Gentamicin |
|---|---|---|---|
| family ridge | 0.771 [0.666, 0.817]; +0.092 [+0.034, +0.145] | 0.904; −0.028 [−0.105, +0.016] | 0.880; −0.007 [−0.031, +0.032] |
| family XGBoost | **0.854 [0.784, 0.892]; +0.257 [+0.138, +0.346]** | 0.938; +0.005 [−0.004, +0.012] | **0.901; +0.025 [+0.014, +0.048]** |

Family XGBoost under-calls resistant isolates in 10.6% (ceftriaxone), 2.2%
(ciprofloxacin) and 1.8% (gentamicin) of right-censored rows. Development
nested CV is unchanged or slightly higher (ridge 0.924, 0.914, 0.893; XGBoost
0.933, 0.917, 0.895), so the grouping costs nothing where training and target
share variants and recovers most of the geographic transfer failure.

**Status.** Post hoc: the analysis was motivated by the external errors,
although the grouping rule is generic and every model was trained on
development only. The family-level catalogue model is the leading hypothesis
for a pre-registered test on a new, independent external set.
