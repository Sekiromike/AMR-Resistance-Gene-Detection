---
name: curate-amr-cohort
description: Acquire, normalize, reconcile, interpret, and validate isolate-level quantitative antimicrobial susceptibility data linked to E. coli genomes. Use for NCBI AST Browser or BigQuery exports, BioSample and assembly joins, MIC or disk-zone handling, breakpoint versioning, duplicate resolution, exclusion ledgers, external-cohort assignment, and cohort freeze decisions.
---

# Curate AMR Cohort

Read [cohort-contract.md](references/cohort-contract.md) before editing cohort
logic or data. Use official source documentation and record retrieval time,
query text, source release, byte size, and SHA-256.

## Workflow

1. Acquire quantitative AST and isolate metadata into immutable source files.
2. Run `scripts/normalize_ncbi_ast.py`; preserve submitted labels separately.
3. Join BioSample, assembly, project, date, site, clinical indication, and source
   provenance without inferring missing values. Specimen source is not a proxy
   for clinical indication.
4. Freeze one explicit breakpoint artifact and checksum. Independently derive
   categories with censoring-aware MIC or zone logic; retain discordances.
5. Reconcile duplicate isolate-drug tests with a prespecified rule and exclusion
   ledger. Never silently choose a convenient result.
6. Add lineage and genomic near-neighbor clusters from sequence analysis.
7. Prespecify external membership from an independent source/site/project.
8. Run `scripts/validate_research_cohort.py --study-config config/study.json`.
9. Freeze cohort and split hashes before feature extraction.

Stop if quantitative values, comparator, method, breakpoint artifact, or
external provenance is unavailable. A submitted S/I/R value alone is not a
publication-grade reference standard.
