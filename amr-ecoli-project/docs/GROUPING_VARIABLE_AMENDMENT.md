# Evaluation Grouping Variable Amendment

Amendment: 004
Date: 2026-09-21
Status: ADOPTED
Decision owner: Rushath Rajeev
Depends on: [`ENDPOINT_AMENDMENT.md`](ENDPOINT_AMENDMENT.md) (001),
[`EXTERNAL_THRESHOLD_AMENDMENT.md`](EXTERNAL_THRESHOLD_AMENDMENT.md) (003)

## Decision

`config/study.json` `splitting.group` changes from `genomic_cluster` to
`lineage_group`.

The two fields now have separate, non-overlapping jobs:

| Field | Definition | Job |
|---|---|---|
| `genomic_cluster` | connected components over pairwise ANI at the frozen rule | **leakage control** — no component may span a split boundary |
| `lineage_group` | MLST sequence type, escalating to cgMLST where ST is uninformative | **cross-validation grouping** |

Both remain required for every isolate. Neither substitutes for the other.

## Why

Single-linkage connected components over pairwise ANI are not lineages. The
external clustering run (`64670770`, 3,159 isolates, 9,967,661 comparisons)
chains transitively: A~B and B~C merge A and C even when A and C are not
similar. The result is one giant component plus a long tail.

Measured on that run, sweeping the thresholds already frozen in
`disjointness_contract.json`:

| Rule | Clusters | Largest | % of cohort | Singletons | Clusters covering 50% |
|---|---:|---:|---:|---:|---:|
| ANI ≥ 99.9, AF ≥ 0.9 | 453 | 1,314 | 41.6% | 342 | **2** |
| ANI ≥ 99.99, AF ≥ 0.95 | 1,242 | 978 | 31.0% | 1,096 | 9 |
| ANI ≥ 99.999, AF ≥ 0.98 | 2,928 | 108 | 3.4% | 2,867 | 1,349 |

There is no usable operating point. Loosening produces an indivisible group
holding a third to a half of the cohort; tightening dissolves the structure
into 98% singletons, which groups nothing. A five-fold grouped
cross-validation on `genomic_cluster` would have to place roughly 40% of
isolates in a single fold.

Development clustering (`64633773`) is expected to chain at least as strongly,
because chaining opportunity grows with cohort size and that cohort is 10,334
genomes rather than 3,159.

**Outcome (recorded 2026-09-27, prediction not borne out).** Development
clustering `64633773` (106,729,566 comparisons) chains far less than the
external set: the largest component holds 739 of 10,334 isolates (7.2%),
against 1,314 of 3,159 (41.6%) externally. The external set is a single
national surveillance collection dominated by a few clones; the development
set pools many unrelated submitters. The decision stands regardless: the
external components are still degenerate, and a component is still not a
lineage.

## What this does not change

Over-merging is **correct** for leakage control. A conservative component that
groups more isolates than a strict lineage definition would is exactly what is
wanted when asking "could this isolate leak across a split boundary?" The
existing constraint stands unchanged and remains enforced:

- no `genomic_cluster` may span the development/external boundary;
- no `deduplication_group` may span an evaluation boundary;
- the putative-duplicate rule continues to define external exclusion
  (amendment 003).

The ANI components are therefore kept, and kept required. Only their use as the
grouping variable is withdrawn.

## Consequence: MLST becomes a hard prerequisite

`lineage_group` was previously a required field with no implementation.
Promoting it to the grouping variable makes MLST a blocking dependency for
every split, every fold, and therefore every model.

- `genomics_mlst_assembly` already exists in `workflow/rules/genomics.smk` and
  now runs in its own `genomics-mlst` environment, verified to solve.
- `lineage_group` is the MLST sequence type. Where MLST returns a novel or
  unresolved ST, the isolate is **not** given a placeholder; it is routed to
  the exclusion ledger, or the lineage definition is escalated to cgMLST for
  the whole cohort and the escalation recorded. Mixing ST and cgMLST
  definitions within one analysis is not permitted.
- Grouping quality must be reported before modelling: number of lineages,
  size distribution, largest lineage as a fraction of the cohort, and the
  number of lineages required to cover half the cohort. If MLST shows the same
  degenerate concentration as the ANI components, this amendment is reopened
  rather than worked around.

## Sensitivity analysis, prespecified

Because the grouping variable determines every performance estimate, the
locked evaluation reports a sensitivity analysis over grouping definitions:
MLST ST as primary, and at minimum one coarser and one finer alternative
(phylogroup, and cgMLST or ANI components) as secondary. Divergence between
them is reported, not reconciled silently.

## MLST adequacy criteria (ADOPTED 2026-09-27 by the decision owner)

Written before any MLST output exists, so the rule cannot be chosen to suit
the answer. Applied once, to the development isolates that pass genome QC,
separately for the primary (human clinical) and sensitivity (all sources)
populations. No phenotype is read. External isolates are not used to decide.

"Resolved" means an Achtman ST with MLST status PERFECT or OKAY.

| Check | Measure | Rule |
|---|---|---|
| 1. Typing coverage | resolved isolates / QC-passing isolates | **≥ 95%**: keep ST; unresolved isolates go to the exclusion ledger with reason `lineage_unresolved`, counted per drug. **< 95%**: escalate the whole cohort to cgMLST; do not exclude. |
| 2. Concentration | largest ST as a fraction of the population | **> 40%** (the degenerate level seen in the external ANI components): reopen this amendment; do not work around it. |
| 3. Effective diversity | inverse Simpson index 1/Σpᵢ² over STs | **< 10** (fewer than two effective lineages per fold at 5 folds): reopen this amendment. |

Always reported, never used to choose: number of STs, size distribution,
largest ST fraction, STs covering 50%, and resulting grouped-fold sizes. A
single ST above one fold's share (20%) is permitted, since a fold holding most
of ST131 is a legitimate test of generalisation to an unseen lineage, but it is
reported.

Escalation (check 1) means cgMLST for every isolate. It is not yet
implemented and would block all splits until it is.

Implemented in `scripts/assess_lineage_adequacy.py` with thresholds in
`config/study.json` `splitting.lineage_adequacy`. The most severe outcome over
the two populations wins (escalate, then reopen, then keep), and any outcome
other than keeping ST exits non-zero so no split can proceed past it.

## Voided application: wrong MLST scheme (2026-09-27)

Cohort build `64927962` applied the criteria to MLST calls from genome-analysis
array `64927454` and stopped with `escalate_cgmlst` (typing coverage 59% of
all sources, 72% of human clinical). That result is **void**: the array ran
`mlst --scheme ecoli`, which in mlst ≥ 2.23 is the **Pasteur** scheme (dinB,
icdA, pabB, polB, putP, trpA, trpB, uidA), not the Achtman scheme this
amendment specifies (adk, fumC, gyrB, icd, mdh, purA, recA; mlst name
`ecoli_achtman_4`). The finalizer's scheme check expected the pre-2.23 name
and so accepted it; the six-genome canary checked only that typing succeeded.

Recorded so the voided run is not hidden:

- Its most frequent Pasteur STs were 43 (Achtman ST131's Pasteur counterpart),
  2, 1, 3 and 4. Development genotypes only; no phenotype was read and no
  external lineage distribution was inspected.
- The criteria are **not** changed in response. They are re-applied, unchanged,
  to Achtman calls.

Remediation: `config/study.json` `genomics.mlst` now freezes the scheme and its
seven loci; `finalize_genome_analyses.py` stops the run if any report names a
different scheme or loci; `hpc/unity/mlst-array.sbatch` re-types every
verified assembly into `results/genomics/mlst/ecoli_achtman_4/`, leaving the
Pasteur outputs untouched as a record; the cohort build reads that tree.

## Reopened: near-identical genomes split across STs (2026-09-27)

With Achtman typing (re-typing arrays `64939241`/`64939277`), cohort build
`64939242` passed the adequacy criteria (`keep_st`; typing coverage 99%;
largest ST131 at 10.7% of all sources and 23.8% of human clinical; inverse
Simpson 43.4 and 12.4). It then stopped at amendment 005 rule G: **929 of
19,771** development putative-duplicate pairs (ANI ≥ 99.99, AF ≥ 0.95) join
isolates of different STs, above the 1% reopening trigger. As prespecified,
this amendment is reopened rather than patched.

Diagnosis (development genotypes only; no phenotype read):

- Of cross-ST pairs with both isolates typed, 970 of 1,007 differ at exactly
  one of the seven loci and 37 at two; none differ at more.
- Typical case: ST131 paired with ST9126, a single-locus variant of ST131 that
  has no clonal complex assigned in the scheme table.
- So the defect is ST granularity at its edges: one mutation in an MLST gene
  gives a near-clone a new ST number. It is not a failure of MLST as a
  lineage signal.

Candidate groupings, measured on typed, eligible development isolates:

| Grouping | Groups | Largest, all sources | Largest, human clinical | Inverse Simpson, human clinical |
|---|---:|---:|---:|---:|
| A. ST (current) | 1,222 | 10.7% | 23.8% | 12.4 |
| B. STs linked by single-locus variation | 492 | 11.0% | 24.6% | 10.9 |
| C. ST joined by duplicate pairs | 1,159 | 10.8% | 24.0% | 12.0 |

**Adopted 2026-09-27 by the decision owner: B.** Implemented in
`finalize_genome_analyses.assign_slv_groups`; groups are named `SLV:ST<n>`.
`lineage_group` becomes the single-locus-variant group: STs observed in
development are linked when their Achtman profiles differ at exactly one
locus, and each connected component is one lineage group, named by its most
frequent ST. This is the classical eBURST clonal-complex definition, needs no
ANI data and no tuning, and absorbs the near-clone splits at their cause.
Remaining cross-group duplicate pairs (the two-locus cases, at most 37 of
20,490 pairs, 0.18%) stay under rule G unchanged. The adequacy criteria are
re-applied to B unchanged; B passes all three.

C was not chosen because it defines lineage from the same ANI edges that rule
G audits, which amounts to lifting the 1% trigger.

## Open

- ~~Development clustering must still complete~~ — done; see Outcome above.
- ~~MLST adequacy criteria must be adopted before MLST results are read~~ — adopted 2026-09-27, before any cohort-level MLST output was read.
- Phylogroup, required for the prespecified sensitivity analysis, is not yet
  computed; the genome-analysis array produces ST only.
- ~~Frozen NCBI source files must be transferred to Unity~~ — done; the
  genome-analysis array runs on the 13,493 checksum-verified disjointness
  assemblies (`hpc/unity/genome-analysis-array.sbatch`).
