# Deduplication and Repeat-Measurement Policy

Amendment: 005
Date: 2026-09-23
Status: ADOPTED (approved by the decision owner 2026-09-23)
Decision owner: Rushath Rajeev
Depends on: 001 (MIC endpoint), 003 (external duplicate rule), 004 (lineage grouping)
Resolves: blocker 6 in [`COHORT_REBUILD_READINESS.md`](COHORT_REBUILD_READINESS.md)

This policy is written before any phenotype is modelled and uses only
sequence identity, identifiers, and the structure of the measurements
themselves. It never selects a value because it is convenient.

It covers three separate problems that the cohort contract previously folded
into one field:

| Problem | Question | Governs |
|---|---|---|
| R | One isolate, one drug, several AST rows | the MIC value used |
| I | Several records that are the same physical isolate | `deduplication_group` |
| G | Distinct isolates that are genomically the same clone | fold assignment |

## Measured scope in the frozen development snapshot

| Check | Result |
|---|---|
| Isolate-drug groups with >1 AST row | **11 groups, 27 rows, 4 isolates** |
| …of which MIC-only repeats | 6 groups on 2 isolates |
| …of which one MIC + one disk-zone row | 5 groups on 2 isolates |
| MIC repeat groups whose intervals are mutually consistent | **6 of 6** |
| Isolates sharing a BioSample accession | 0 |
| Isolates sharing an assembly accession | 0 |
| Isolates sharing an SRA run (7,307 runs checked) | 0 |

Rule R therefore affects 2 isolates, and rule I currently merges nothing. The
substantive question is rule G.

## Rule R — repeated measurements for one isolate and drug

**R1. Endpoint filter first.** Only MIC-eligible rows (amendment 001) enter
reconciliation. Disk-zone rows are not part of the MIC endpoint. The 5 groups
that pair one MIC with one zone row therefore reduce to a single MIC row and are
not repeats.

**R2. Reconcile by interval intersection.** Each MIC row denotes an interval on
the concentration scale:

| Reported | Interval |
|---|---|
| `= v` | [v, v] |
| `<= v` | (0, v] |
| `< v` | (0, v) |
| `>= v` | [v, ∞) |
| `> v` | (v, ∞) |

The reconciled measurement is the **intersection of all intervals** in the
group: the most specific statement consistent with every submitted
measurement. This adds no information and discards none that is consistent.

| Group (observed) | Intervals | Reconciled |
|---|---|---|
| `>32`, `>=64` | (32, ∞) ∩ [64, ∞) | `>=64` |
| `>2`, `>=4`, `>2` | (2, ∞) ∩ [4, ∞) | `>=4` |
| `<=4`, `<=1`, `<=1` | (0, 4] ∩ (0, 1] | `<=1` |

**R3. Exclusions.** A group is excluded to the ledger, never resolved by
choosing a row, when:

- the intersection is empty (`discordant_repeat_measurements`), for example
  `= 2` against `= 8`;
- the intersection is a bounded two-sided interval that is not a single value
  (`repeat_not_representable`), for example `<= 4` with `>= 2`, because the
  endpoint schema carries one comparator and one value;
- units differ within the group (`repeat_unit_conflict`).

**R4. Provenance.** The reconciled row carries the lowest source row as its
record of origin, plus a `reconciled_from` field listing every source AST
record identifier in the group. The original rows stay in the normalized
source ledger unchanged.

Applied to the frozen snapshot, R keeps all 6 MIC repeat groups (2 isolates ×
3 drugs) that the current rule excludes, and excludes none.

## Rule I — identifier-level duplicates

Isolate records sharing a BioSample accession, an assembly accession, or any
SRA run accession are the same physical isolate and receive the same
`deduplication_group`. This rule currently merges nothing, since there are zero
shared identifiers, but it stays in force as a guard for future source versions.

Patient and outbreak identifiers do not exist in any public source used here.
They are **not** proxied by BioProject, submitter, collection date, or
geography. The only outbreak control available is genomic (rule G), and this is
stated as a limitation.

## Rule G — genomic duplicates between distinct isolates

A **genomic duplicate pair** is two development isolates meeting the frozen
putative-duplicate rule: ANI ≥ 99.99% and minimum reciprocal aligned fraction
≥ 0.95. These are the same values already used for external exclusion
(amendment 003), so no new threshold is introduced.

**Why rule G must not use connected components.** Grouping duplicate pairs by
transitive closure reproduces the chaining that amendment 004 removed. On the
external cohort, connected components at exactly this threshold put **978 of
3,159 isolates (31%)** into one group. Merging those groups into cross-validation
would collapse the folds again.

**G1. Within-lineage duplicates need no action.** At ≥ 99.99% ANI, two genomes
differ at roughly 500 bases or fewer across ~5 Mb and should share an MLST
sequence type. Lineage grouping (amendment 004) already keeps them in the same
fold.

**G2. Cross-lineage duplicates are anomalies, handled pair by pair.** A genomic
duplicate pair whose members have **different** `lineage_group` values points
to an MLST calling error, a novel or ambiguous ST, or contamination. Every such
pair is written to a ledger. Pairs are processed in a fixed order (descending
ANI, then ascending isolate identifiers) and each is resolved in one of two
ways:

- **Merge** — both members receive a shared `deduplication_group`, which
  unions their two lineages into one cross-validation group through the
  existing split code, **only if** the resulting group stays at or below
  **10% of development isolates**. That is half a fold at five folds, which
  leaves room to balance class counts.
- **Exclude** — otherwise, the member belonging to the smaller lineage is
  excluded as `cross_lineage_genomic_duplicate`. On a tie, the member with the
  lexicographically larger isolate identifier is excluded. The other member is
  retained.

This bounds the damage from anomalies. A few anomalous pairs cannot chain
lineages into a giant group, and every decision is deterministic and
ledgered.

**G3. Reporting.** Before any model is fitted, report the number of genomic
duplicate pairs, the within- versus cross-lineage split, merges performed,
isolates excluded, and the largest cross-validation group as a fraction of
development. If cross-lineage pairs are not rare (more than 1% of duplicate
pairs), the MLST-based lineage definition is suspect and amendment 004 is
reopened rather than patched.

## Consequence for evaluation: the bootstrap unit

`config/study.json` specifies cluster-bootstrap confidence intervals. The
external cohort has the same chaining problem: its largest ANI component holds
1,314 of 3,159 isolates. A bootstrap that resamples ANI components would
therefore resample one unit carrying 42% of the data, which makes the intervals
degenerate.

**The bootstrap resampling unit is `lineage_group`** for both development and
external evaluation, consistent with amendment 004. ANI components are kept for
leakage control only.

## What this proposal deliberately does not do

- It does not dereplicate the cohort to one representative per clone.
  Removing near-identical isolates would change the class balance and discard
  real repeat-measurement information, and dereplication's own clustering
  step would chain.
- It does not reconcile across endpoints. MIC and disk-zone values are never
  merged.
- It does not reconcile across isolates. Genomic duplicates keep their own MIC
  values; similar genomes with different MICs are measurement information, not
  errors.

## Noted while measuring

- `PDT000041778.1`, one of the two isolates with triplicate MIC submissions, is
  the same development isolate as the first confirmed cross-cohort collision
  (99.99% ANI to `JBABADA-19-0001`). Repeated submission of one isolate is
  typical of a reference or quality-control strain.
- `PDT001463387.1` reports a ciprofloxacin MIC of `= 2.5` mg/L, which is not on
  the standard doubling-dilution series. This is direct evidence of the
  method heterogeneity recorded in amendment 001. It is retained, and it is
  exactly the kind of value the heterogeneity analysis must stratify.
- `data/interim/ncbi_ast_normalized.csv` on disk predates the `mic_eligible`
  column and must be regenerated before the cohort rebuild.

## Implementation

| Rule | Where | Tests |
|---|---|---|
| R | `construct_research_cohort.py`: `prepare_mic_rows`, `reconcile_mic_repeats`; MIC endpoint only, the categorical secondary keeps its blanket exclusion | every R2/R3 case, plus all six real snapshot repeat groups |
| I | `build_isolate_metadata.py`: identifier union across BioSample, run/assembly accession | shared-BioSample merge; cross-split identifier refused |
| G | `build_genomic_disjointness.py extract-duplicate-pairs` (streaming) feeds `build_isolate_metadata.py` | within-lineage no-op, symmetric pairs counted once, merge under cap, exclusion over cap, tie-break, reopen trigger, and an anti-chaining test |
| Bootstrap | `evaluate_predictions.py` `BOOTSTRAP_UNIT = "lineage_group"`; predictions carry `lineage_group` | a single chained component no longer breaks resampling |

Latent defects fixed alongside, each now covered by a test:

- Under `mic_regression` the cohort builder copied categorical-only metadata
  (`site`) and a breakpoint hash, which would have raised `KeyError` on the
  first real MIC build.
- The metadata builder read `specimen_source` from the sequence-request table,
  which has no such column, so it would have excluded every isolate. It now
  takes an explicit specimen-source mapping.

The tiny end-to-end workflow passes all 8 steps and reports
`bootstrap_unit: lineage_group`.
