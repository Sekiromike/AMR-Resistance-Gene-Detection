# Cohort Rebuild Readiness

Date: 2026-09-21
Endpoint: `mic_regression` (amendment 001)
External rule: putative-duplicate (amendment 003)

What must exist before `data/curated/cohort.csv` can be built, and the current
state of each. Nothing here may be substituted by a proxy.

## Ready

| Input | State |
|---|---|
| Normalized AST | `data/interim/ncbi_ast_normalized.csv`, 25,735 rows, hash-frozen |
| Quantitative measurements | 25,466 MIC rows, comparators preserved, 71.22% censored |
| Frozen external membership | 1,888 locked / 1,271 excluded / 536 lineage-novel, SHA-256 `b4c6baf7…` |
| Cohort builder | `scripts/construct_research_cohort.py --endpoint-profile mic_regression` |
| Metadata builder | `scripts/build_isolate_metadata.py` |
| Cohort validator | `scripts/validate_research_cohort.py --endpoint-profile mic_regression` |
| Assemblies | 13,493 built and checksum-verified on Unity |

Under `mic_regression` the builder no longer requires `ast_method`,
`ast_testing_date`, `clinical_indication`, `ast_category`, the breakpoint
fields, `site`, or `surveillance_network`. It consumes the **normalized** AST
directly; the breakpoint interpretation step is skipped.

## Blocking

| # | Blocker | Needed for | Status |
|---|---|---|---|
| 1 | Development genomic clusters | `genomic_cluster` | search done; normalization OOM-killed; streaming resume job 64793946 queued behind maintenance |
| 2 | Within-external clustering | `genomic_cluster` for the 1,888 locked externals | complete (job 64670770) |
| 3 | Species confirmation + contamination screen | `species_method`, `species_result` | genomics phase unimplemented |
| 4 | Assembly QC status | `genome_qc_status` | rules exist in `genomics.smk`, never executed |
| 5 | Lineage assignment (MLST/cgMLST) | `lineage_group` | not implemented |
| 6 | Deduplication rule | `deduplication_group` | adopted and implemented (amendment 005); needs duplicate pairs from the resumed development run |
| 7 | AMRFinderPlus database release | genomics DAG gate | pinned `2026-08-07.1` |

### Found on 2026-09-23

| # | Blocker | Status |
|---|---|---|
| 8 | `country` blank for all 25,735 rows | resolved by amendment 006: INSDC prefix, 25,403 rows filled, 332 absent |
| 9 | Specimen-source mapping | resolved: `build_source_attributes.py`; 126 development and 29 external absent |
| 10 | Stale normalized AST | resolved: regenerated (`653c89a7…`); a legacy view reproduces the audited `e72661c6…` exactly |
| 11 | Prevalence baseline and evaluator are still categorical (S/R, sensitivity); the MIC endpoint needs essential-agreement and interval-censored metrics — **in progress**: `scripts/evaluate_mic_predictions.py` implements the rules in `docs/MIC_EVALUATION_SPEC.md` (PROPOSED); the no-genome interval-censored baseline and the off-scale exclusion wait on its adoption | implementation |
| 12 | Intended-use population | resolved by amendment 006: primary human clinical (3,935 isolates), all-source sensitivity |
| 13 | ~~JARBS MICs are not yet normalized into the cohort schema~~ **Resolved**: `scripts/normalize_jarbs_ast.py` reshapes the frozen table (7,605 isolate-drug rows from 2,535 isolates, all MIC-eligible, country from the same INSDC rule; structural counts only, the external set stays sealed) and `construct_research_cohort.py` accepts repeated `--interpreted-ast`, recording each source's path, hash and row count | implementation |

### Blocker 2 was found by writing the builder

The comparison in job `63772662` was development-versus-external, and job
`64633773` is development-versus-development. **Neither clusters external
isolates against each other.** The cohort contract forbids a genomic cluster
spanning the development/external boundary and requires every isolate to carry
a cluster, so the 1,888 locked externals need a within-external clustering run.

It is cheap relative to what has already run: 3,159 queries against a
3,159-genome sketch, versus the 10,334-query job now executing. The existing
`hpc/unity/development-clusters.sbatch` pattern covers it with an external
sketch and `--query-cohort` handling.

`build_isolate_metadata.py` fails closed on this with reason code
`genomic_cluster_absent` rather than inventing an assignment, and a test pins
that behaviour.

### Blocker 6 resolved by amendment 005

Adopted 2026-09-23 in [`DEDUPLICATION_POLICY.md`](DEDUPLICATION_POLICY.md) and
implemented. Repeated MICs are reconciled by interval intersection, identifier
duplicates share a deduplication group, and cross-lineage genomic duplicates
are merged under a 10% cap or excluded. The candidate rule once written here --
treating a putative-duplicate component as one group -- was rejected because
that threshold chains 978 of 3,159 external isolates into one component.

## Order of operations once unblocked

1. Development clustering resume completes (`64793946`), then extract
   duplicate pairs (`development-duplicate-pairs.sbatch`).
2. Within-external clustering: done (`64670770`).
3. Freeze the AMRFinderPlus release and run the genomics DAG over the 13,493
   assemblies, producing species, QC, and lineage (blockers 3–5, 7).
4. Resolve blockers 8 (`country`) and 12 (intended-use population).
5. `build_isolate_metadata.py` → `data/curated/isolate_metadata.csv`.
6. `construct_research_cohort.py --endpoint-profile mic_regression`.
7. `validate_research_cohort.py --endpoint-profile mic_regression`.
8. Publish cohort flow, missingness, per-drug and per-split counts, censoring
   profile, and a per-drug precision calculation.
9. Only then: prevalence and AMRFinderPlus baselines — the first result.

## Not blocking, deliberately

Site, clinical indication, surveillance network, AST testing date, and
testing-standard version remain absent from every public source audited. They
are not required under `mic_regression` and must not be reconstructed. The
categorical endpoint stays blocked under amendment 001.
