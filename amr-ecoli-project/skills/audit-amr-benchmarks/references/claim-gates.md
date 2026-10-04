# AMR Claim Gates

## Dataset Gate

Require quantitative AST provenance, a frozen interpretation standard, duplicate
reconciliation, genome QC, species confirmation, lineage/near-neighbor groups,
and a prespecified untouched external source.

## Evaluation Gate

Require one task per antibiotic, identical cohorts across models, development-
only tuning, clinically justified locked thresholds, external evaluation once,
and cluster-bootstrap confidence intervals. Report sensitivity, specificity,
false-susceptible and false-resistant rates, PPV, NPV, MCC, balanced accuracy,
AUROC, AUPRC, Brier score, calibration, and coverage under abstention.

## Model Gate

Require prevalence, lineage/neighbor, AMRFinderPlus, AMRrules/ResFinder/
PointFinder/RGI, regularized curated-feature, and pangenome/unitig baselines before
crediting a foundation model. Match data, isolates, endpoints, and selection
budget. Audit pretraining overlap.

## Discovery Gate

Require lineage-aware bacterial GWAS, stability across folds/sites/lineages,
independent replication, genomic-context review, and targeted laboratory work
before mechanistic or causal language.
