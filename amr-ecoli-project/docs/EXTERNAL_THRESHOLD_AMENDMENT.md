# External Exclusion Threshold Amendment

Amendment: 003
Date: 2026-09-21
Status: ADOPTED
Decision owner: Rushath Rajeev
Amends: [`EXTERNAL_COHORT_AMENDMENT.md`](EXTERNAL_COHORT_AMENDMENT.md) (002)
Depends on: [`ENDPOINT_AMENDMENT.md`](ENDPOINT_AMENDMENT.md) (001)

## Decision

The rule that removes an external isolate from the locked evaluation set changes
from the **near-neighbour** threshold to the **putative-duplicate** threshold.

| | ANI | min reciprocal AF | Externals removed | Externals retained |
|---|---:|---:|---:|---:|
| Amendment 002 (near-neighbour) | ≥ 99.9 | ≥ 0.90 | 2,623 | 536 |
| **Amendment 003 (duplicate)** | **≥ 99.99** | **≥ 0.95** | **1,271** | **1,888** |

Resulting cohort structure:

| Set | n | With all three target MICs | Role |
|---|---:|---:|---|
| **External-locked** | **1,888** | **1,555** | the one locked external evaluation |
| ↳ lineage-novel stratum | 536 | 454 | nested, strictest; isolates with no near neighbour at all |
| External-excluded | 1,271 | — | putative genomic duplicates of development isolates |

Both thresholds are unchanged values already frozen in
`disjointness_contract.json`. No new threshold is introduced; the amendment
changes only which of the two prespecified rules governs exclusion.

## Why the previous operating point was wrong

The `NEAR_NEIGHBOR` class in comparison `63772662` spans **ANI 99.90 to 100.00,
with a median of 99.93 and a 90th percentile of 99.97**. For *E. coli*, ANI at
that level is ordinary similarity *within a sequence type or clonal complex*.
It indicates shared lineage, not a shared isolate.

Using it as the exclusion rule therefore did not remove leakage. It removed
**lineage membership**, with three consequences:

1. **It discards the clinically dominant lineages.** Any Japanese isolate
   belonging to a globally disseminated lineage such as ST131 has a 99.9%-ANI
   relative somewhere among 10,334 globally sampled genomes. Those are the
   lineages that cause most human infection and are exactly the ones a
   susceptibility predictor must handle.
2. **It silently changes the research question.** A cohort of isolates whose
   lineages are absent from training answers "can the model predict on an
   unseen lineage?" That is a legitimate question and remains a stratum below,
   but it is not external validation and must not be reported as though it
   were.
3. **It costs the power needed for the primary endpoint.** See below.

The duplicate threshold is the one that means *the same clone*, which is the
leakage that actually invalidates an external evaluation.

## Power

No formal calculation has been performed, and the external set stays sealed, so
the resistance prevalence of either cohort is **not known and was deliberately
not inspected**. As an illustration only, at an assumed 20% resistance a
454-isolate cohort yields roughly 90 resistant isolates per drug and a 95%
confidence half-width near ±8 percentage points on sensitivity; 1,555 isolates
yield roughly 310 and a half-width near ±4. Lower-prevalence drugs are worse in
both cases.

A formal per-drug precision calculation is required before the locked
evaluation, using development prevalence as the planning assumption. Any drug
failing it is reported as underpowered rather than quietly included.

## Checks performed

**Assembly quality is not a confound.** A concern that collision-free isolates
might merely be fragmented assemblies — a fragmented genome has lower aligned
fraction and so escapes an AF-based rule — was tested against the frozen JARBS
assembly metrics and is not supported:

| Group | n | Median contigs | Median coverage | Median total bases |
|---|---:|---:|---:|---:|
| Collision-free | 536 | 188 | 53 | 5,121,504 |
| Colliding | 2,623 | 191 | 52 | 5,140,678 |

**Not yet checkable.** Whether the 536 are lineage-atypical requires MLST or
cgMLST, which the genomics phase has not yet produced. This must be reported
when it becomes available.

## Analysis consequence

The headline external analysis becomes **performance as a function of genomic
distance to the training set**, with the lineage-novel stratum as the most
distant band. This uses the disjointness computation as a covariate rather than
only as a filter, and it measures the out-of-distribution reliability question
that is this study's stated contribution, instead of reducing it to a single
small cohort.

## What does not change

- Putative genomic duplicates remain excluded. Clone-level leakage control is
  unchanged.
- The external set stays sealed from feature, model, threshold, calibration and
  epoch selection. No external phenotype has been read.
- Selection remains sequence-only and phenotype-blind, computed before any
  phenotype was available.
- JARBS sampling remains resistance-enriched, so its prevalence does not
  estimate routine prevalence and predictive values must be reported
  prevalence-adjusted.
- JARBS still publishes no AST testing date, hospital, or clinical indication,
  so leave-site-out evaluation remains impossible and must not be simulated.
- The categorical S/I/R endpoint remains blocked under amendment 001.

## Open

- The 166 collision components from `63772662` are still recorded as
  `UNRESOLVED`. This amendment excludes the duplicate-bearing isolates rather
  than adjudicating each component.
- Lineage composition of both strata, pending genomics.
- Per-drug precision calculation, pending the rebuilt cohort.
