# Current MIC study protocol

Consolidated on 4 October 2026 from the adopted analysis specifications and
completed experiments. This document summarizes the executed study; it is not a
new prospective registration. The dated amendments and analysis plans linked
below record the decisions that governed the experiments.

## Research question

Do frozen Evo 2 representations of annotation-selected resistance loci improve
prediction of ceftriaxone, ciprofloxacin, and gentamicin MICs in *Escherichia
coli* relative to curated AMRFinderPlus allele, family, and drug-class features?

## Endpoint and populations

Prefer quantitative MIC prediction with interval censoring when enough MIC data exist;
derive categorical performance as a downstream evaluation.

The executed primary endpoint is interval-censored log2 MIC, with essential
agreement as defined in [MIC_EVALUATION_SPEC.md](MIC_EVALUATION_SPEC.md).
Clinical S/I/R validation was not possible because required assay metadata and
independent breakpoint review were unavailable. See the adopted
[endpoint amendment](ENDPOINT_AMENDMENT.md).

The final cohort contains 29,619 isolate–antibiotic rows from 11,520 isolates:
10,015 development isolates, including 3,849 primary human-clinical isolates,
and 1,505 Japanese JARBS external isolates. All-source development is a
sensitivity population. Source snapshots, exclusions, metadata, genome QC,
split manifests, and hashes accompany the release.

The [population amendment](COUNTRY_AND_POPULATION_AMENDMENT.md),
[external-cohort amendment](EXTERNAL_COHORT_AMENDMENT.md),
[duplicate policy](DEDUPLICATION_POLICY.md), and
[grouping amendment](GROUPING_VARIABLE_AMENDMENT.md) define cohort construction.
Missing AST method, testing date, breakpoint version, site, and indication are
kept missing.

## Features, models, and evaluation

The executed model comparisons include constant and nearest-relative baselines,
interval-censored ridge regression, and XGBoost AFT. Feature sets include
AMRFinderPlus allele indicators, frozen Evo 2 locus embeddings, their
combination, and allele + family + drug-subclass indicators.

Nested development folds group lineages, genomic neighbors, and
identifier-linked duplicates. PCA, feature selection, and hyperparameter
selection use training data within the applicable folds. The adopted
[evaluation design](EVALUATION_DESIGN_AMENDMENT.md) and
[foundation-model protocol](FOUNDATION_MODEL_PROTOCOL.md) specify the methods.
Whole-lineage bootstrap intervals condition on the saved predictions.

## Analysis chronology and evidence tiers

1. **Original external comparison:** predictions were frozen and locked before
   reference attachment and scoring. The [external analysis
   plan](EXTERNAL_EVALUATION_PLAN.md) and
   [prediction lock](../config/external_prediction_lock.tsv) govern this result.
2. **Family-feature external follow-up:** designed after inspecting the original
   external errors. Its separately locked predictions do not make it
   confirmatory. See [the exploratory plan](POSTHOC_FAMILY_ANALYSIS.md) and
   [family prediction lock](../config/external_posthoc_lock.tsv).
3. **Variant withholding:** hypotheses were specified before these runs using
   development data. See [LEAVE_VARIANT_OUT_PLAN.md](LEAVE_VARIANT_OUT_PLAN.md).
   This is controlled development evidence, not independent external replication.
4. **Source sensitivity:** specified development predictions are rescored after
   excluding source rows; models are not retrained for these exclusions.

[Study chronology](../paper/provenance/study_chronology.json) preserves original
plan and lock contents, hashes, and local commit dates. These dates do not
constitute independent timestamps or public registration.

## Outputs and limits

The [paper](../paper/manuscript.md), [supplement](../paper/supplement.md),
[verified result exports](../paper/verification/), and
[reproduction guide](../paper/REPRODUCIBILITY.md) document the completed study.
The audit reproduced 54 evaluations, 438 paired comparisons, and all 18
prediction locks.

The representation uses selected resistance loci and frozen Evo 2 embeddings.
No full-genome or fine-tuned Evo 2 result, clinically validated diagnostic,
independent family-model replication, or novel causal resistance mechanism is
claimed. [SCOPE_CLOSEOUT.md](SCOPE_CLOSEOUT.md) records the remaining scientific
limits.
