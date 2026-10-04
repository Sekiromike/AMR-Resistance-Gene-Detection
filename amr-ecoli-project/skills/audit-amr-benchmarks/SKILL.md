---
name: audit-amr-benchmarks
description: Audit AMR genotype-to-phenotype datasets, models, metrics, and scientific claims for leakage, population structure, invalid AST semantics, external-test reuse, confounding, calibration, uncertainty, and reproducibility. Use for code reviews, result acceptance, manuscript claims, model comparisons, frontier-model evaluation, or deciding whether a benchmark is scientifically valid.
---

# Audit AMR Benchmarks

Read [claim-gates.md](references/claim-gates.md). Lead with findings ordered by
severity and cite file/line or artifact/hash evidence.

## Audit Order

1. Define the biological prediction target and intended-use population.
2. Trace every label to quantitative AST, method, comparator, standard, version,
   and independently derived interpretation.
3. Check isolate identity, duplicates, outbreak neighbors, lineage, site, time,
   and external membership before examining metrics.
4. Prove preprocessing, feature selection, calibration, thresholding, and early
   stopping use development data only and run inside folds where applicable.
5. Compare identical isolates and endpoints across prevalence, lineage,
   mechanistic, pangenome/unitig, hybrid, and foundation-model baselines.
6. Verify per-drug clinical errors, calibration, cluster-bootstrap intervals,
   subgroup/OOD performance, and abstention coverage.
7. Audit database/model versions, pretraining overlap, compute cost, hashes, and
   rerun instructions.
8. Reject causal, clinical, or discovery language beyond the evidence tier.

An impressive AUROC cannot rescue a contaminated cohort, random lineage split,
test-set tuning, unverified label, or unmatched baseline.
