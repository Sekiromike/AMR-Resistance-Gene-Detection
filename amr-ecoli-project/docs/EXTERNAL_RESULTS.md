# Locked External Evaluation: Results

Date: 2026-09-29
Plan: [`EXTERNAL_EVALUATION_PLAN.md`](EXTERNAL_EVALUATION_PLAN.md) (adopted)
Predictions: job `65019273`, locked in `config/external_prediction_lock.tsv`
(commit `4c18aa3`, before any external MIC was read)
Scoring: job `65020832` (checksums verified). An earlier attempt, `65020797`,
stopped on a software defect before computing any statistic (empty uncertainty
column in point-prediction files); fixed, tested and rerun on the same locked
files, as the plan allows.

Set: 1,505 JARBS isolates × 3 drugs, one laboratory (MicroScan broth
microdilution, Japan), genomically disjoint from development. Primary models
were trained on human-clinical development. Essential agreement (EA) has 95%
whole-lineage bootstrap intervals.

## Pre-registered results

### Evo 2 credit: not awarded for any drug

Credit condition 2 (external gain of AMRFinderPlus + Evo 2 over AMRFinderPlus
alone, both learners) was tested for gentamicin, the only drug meeting
condition 1:

| Gentamicin, external | Difference | 95% interval |
|---|---:|---|
| ridge, AMR+Evo 2 − AMR | −0.033 | [−0.056, −0.015] |
| XGBoost, AMR+Evo 2 − AMR | −0.021 | [−0.046, +0.000] |

The development gain, traced to one source (PRJNA1297298), did not transfer.
It reversed. Ceftriaxone (+0.005, +0.015) and ciprofloxacin (−0.023, +0.002)
show no demonstrated gain. Evo 2 alone was worse than AMRFinderPlus for every
drug and learner (−0.010 to −0.100).

### Known-mechanism models transfer for ciprofloxacin and gentamicin

| Drug | Model | EA | Resistant under-call | vs constant | vs nearest relative |
|---|---|---|---|---|---|
| ciprofloxacin | ridge AMR | 0.932 [0.894, 0.950] | 0.012 [0.003, 0.026] | +0.532 | +0.082 |
| ciprofloxacin | XGBoost AMR | 0.933 [0.913, 0.954] | 0.025 [0.013, 0.064] | +0.533 | +0.084 |
| gentamicin | ridge AMR | 0.887 [0.841, 0.916] | 0.057 [0.031, 0.135] | +0.104 | +0.135 |
| gentamicin | XGBoost AMR | 0.876 [0.829, 0.909] | 0.121 [0.071, 0.269] | +0.093 | +0.124 |
| ceftriaxone | ridge AMR | 0.680 [0.539, 0.742] | 0.258 [0.181, 0.480] | +0.390 | +0.151 |
| ceftriaxone | XGBoost AMR | 0.597 [0.480, 0.692] | 0.513 [0.460, 0.669] | +0.307 | +0.068 |

Every AMRFinderPlus model beats the no-genome constant with intervals
excluding zero. Ciprofloxacin's resistant under-call (1.2% for ridge) is at
the ~1.5% very-major-error reference level. That is not a clinical claim: the
categorical endpoint is blocked and this is not a device evaluation.
Ceftriaxone carried the pre-declared "underpowered" label. Results trained on
all development sources were similar.

## Diagnosis of the ceftriaxone drop (after scoring)

Ceftriaxone EA fell from 0.92 in development to 0.60–0.68 here, for two
separate reasons:

1. **Geographic ESBL variants absent from development.** Of 710 external
   isolates with ceftriaxone > 64 mg/L, the ridge model placed 119 on the
   susceptible side. Every one carried a beta-lactamase, and the leading
   determinants were **blaCTX-M-2 (42 isolates; 0 in human-clinical
   development), blaCTX-M-8 (21; 2), blaCTX-M-3 (14; 2)**, CTX-M-24, CTX-M-65
   and blaIMP-1 (0). These CTX-M groups are known to be common in Japan and
   rare in the US/UK sources that dominate development. Catalogue features are
   allele symbols filtered at five training isolates, so these carriers look
   gene-negative to the model.
2. **Reporting scale.** External resistant results are ">64" where development
   mostly reported "≥64", so a prediction of 32 mg/L agrees in development but
   under-calls here (155 of XGBoost's 364 errors on ">64" rows). Separately,
   the additive ridge model extrapolates beyond any panel (up to 10⁴–10⁶ mg/L).

## Exploratory finding (post hoc; hypothesis-generating only)

For highly resistant external isolates carrying the unseen variants (102
isolates), the share predicted on the resistant side (> 2 mg/L) was:

| Model (human-clinical trained) | Unseen-variant carriers | Common-variant carriers (546) |
|---|---:|---:|
| ridge AMR | 0.147 | 0.998 |
| XGBoost AMR | 0.127 | 1.000 |
| ridge AMR+Evo 2 | 0.196 | 0.996 |
| XGBoost AMR+Evo 2 | 0.461 | 0.982 |
| ridge Evo 2 alone | **0.784** | 0.907 |
| XGBoost Evo 2 alone | **0.814** | 0.963 |

Evo 2 representations recognised resistance-gene variants that symbol-level
catalogue features could not. This is the property that would make a genomic
foundation model useful. It was defined after scoring, and it has not been
compared with the obvious curated alternative: catalogue features pooled to
gene family (for example "any CTX-M") or to AMRFinderPlus subclass, which may
do as well. It supports a pre-registered test on a new external set. It does
not change the verdict above.

## Conclusions

1. Known-mechanism models transfer from a development set dominated by one US
   hospital to a Japanese laboratory for ciprofloxacin and gentamicin.
2. Evo 2 added no robust accuracy. Its one development gain was
   source-specific and reversed externally.
3. The main external failure is epidemiological: resistance variants common in
   the target region were absent from training, and allele-level catalogue
   features cannot generalise to them.
4. Foundation-model embeddings may address exactly that failure (exploratory,
   above). The way they were combined here, averaged and PCA-reduced
   alongside catalogue features, did not realise it.
