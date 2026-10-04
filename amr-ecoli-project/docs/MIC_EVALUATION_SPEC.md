# MIC Evaluation Specification

Status: ADOPTED 2026-09-27 by the decision owner, before any model was fitted or external MIC read
Decision owner: Rushath Rajeev
Depends on: [`ENDPOINT_AMENDMENT.md`](ENDPOINT_AMENDMENT.md) (001),
[`GROUPING_VARIABLE_AMENDMENT.md`](GROUPING_VARIABLE_AMENDMENT.md) (004),
[`DEDUPLICATION_POLICY.md`](DEDUPLICATION_POLICY.md) (005)
Implements: readiness blocker 11
Code: `scripts/evaluate_mic_predictions.py`,
`tests/test_evaluate_mic_predictions.py`

Amendment 001 fixes essential agreement (EA, prediction within one doubling
dilution) as the primary metric, but does not say how to score a reference MIC
that is only a bound (`<= 0.25`, `> 4`). In the development snapshot that is
18,136 of 25,466 MIC rows (71%). This document fixes that, and the other
scoring details, before any model is fitted and before any external MIC is
read.

## Rules

**Scale.** Work in log2 mg/L. A reported value snaps to the nearest doubling
dilution when within 0.1 log2 of it, so conventional labels (0.12, 0.06, 0.03,
0.015, 0.008) are 2⁻³ … 2⁻⁷. A value further off the scale is not a doubling
dilution and is excluded from the MIC endpoint with reason
`mic_value_off_doubling_scale` (`construct_research_cohort.py`). In development this affects 3 of 25,466 rows
(values 0.2, 2.5 and 200 mg/L); the external count is not inspected.

**Reference set.** Each reference is the set of dilutions it allows:

| Reported | Allowed dilutions |
|---|---|
| `= y` | y |
| `<= y` | y and below |
| `< y` | y−1 and below |
| `>= y` | y and above |
| `> y` | y+1 and above |

**Error.** A prediction is rounded to the nearest dilution (halves up). Its
error is the number of dilution steps to the nearest allowed dilution, and 0
when it is allowed.

**Metrics, per antibiotic.**

| Metric | Definition |
|---|---|
| Essential agreement (**primary**) | error ≤ 1, over **all** rows |
| EA, exact references only | error ≤ 1, over rows with `= y` |
| Exact agreement | error = 0, both ways |
| Over-call ≥ 2 / under-call ≥ 2 | prediction ≥ 2 dilutions above / below the allowed set |
| Bias | mean (prediction − y) over exact references |
| Interval log-likelihood | when a predictive sd is given; dilution d means the true MIC is in (d−1, d] |

**Uncertainty.** 95% percentile intervals from 2,000 bootstrap resamples of
whole `lineage_group`s (config `evaluation`), seed 20260815.

**Reference comparator.** Every model is compared with a no-genome baseline on
identical rows: per antibiotic, the single dilution that maximises essential
agreement on the training folds (ties: higher exact-only agreement, then the
lower dilution), predicted for every held-out row. A model is credited only
for EA above this baseline, with the lineage-bootstrap interval of the
difference.

**Correction, 2026-09-27 (decision owner).** The adopted comparator was an
intercept-only interval-censored normal predicting its mean. Its first run
(`64940232`, development only) fitted degenerately: development MICs are
bimodal and censored at both panel limits, so a single normal flattened to sd
≈ 30 log2 with ceftriaxone mean ≈ 2⁻²⁴ mg/L, and its agreement came almost
entirely from left-censored rows. It was replaced by the best constant above,
the strongest no-genome prediction under the primary metric, before any
genomic model was fitted and without reading any external MIC. The baseline
reports no log-likelihood; log-likelihood comparisons are made only between
models that give predictive distributions.

**Specimen source.** Isolates whose specimen source is absent in source are
kept (the field stays blank, never inferred) and are excluded only from
specimen-type subgroup analyses (decision owner, 2026-09-27).

**Seal.** The evaluator refuses external rows unless run with
`--locked-external`, and a locked run must contain external rows only.

## Why all rows is primary, and what guards it

Scoring exact references only would keep 29% of rows, and not a random 29%:
an MIC is exact only when it falls inside the panel's range, so the clearly
susceptible and clearly resistant isolates, which are the clinically decisive
ones, are mostly the censored ones. EA on exact rows alone would measure
performance on the intermediate band.

Scoring all rows has one known weakness. A model that always predicts the
bottom of the range agrees with every `<= low` row for free. Two things guard
against that:

1. the exact-only EA is always reported beside it; and
2. credit is given only relative to the no-genome baseline, which collects the
   same free agreement on identical rows.

## Not yet specified

- Categorical agreement and very major and major errors: blocked with the
  categorical endpoint (no approved breakpoints).
- Error stratified by distance to the nearest training genome (amendment 001
  secondary): needs the development ANI table joined to fold assignments.
- Abstention coverage: no abstaining model exists yet; any model that
  abstains must report coverage and score only on covered rows, reported
  separately.

## Addendum A (adopted 2026-09-28, before all-sources Evo 2 results or any external MIC were read)

**Reporting by reference class.** Every model, in development and in the
locked external evaluation, is also reported separately for three classes of
reference MIC: left-censored (at or below the panel), exact, and
right-censored (above the panel). For the right-censored class, which holds
the clearly resistant isolates, two figures are reported beside overall
essential agreement: essential agreement within the class, and the
**under-call rate** (prediction at least two dilutions below the reference
set). The under-call is the MIC analogue of a very major error; clinical AST
devices must keep very major errors below about 1.5%. No result may be
described as clinically usable while that rate is above it.

Motivation, from human-clinical development results (job 64983779): the
no-genome baseline scores 100% on left-censored rows and 0% on right-censored
rows, so overall agreement is dominated by susceptible isolates, and the
AMRFinderPlus models under-call 4% (ciprofloxacin) to 19% (gentamicin) of
right-censored rows.

## Addendum B (adopted 2026-09-28, before any external MIC was read)

**Source robustness.** Public MIC data mix laboratories whose panels and
reporting conventions differ, and one submitting project can span many
lineages, so lineage-grouped cross-validation does not stop a model from
recognising a source. Every per-drug development result and every paired
model comparison is therefore also reported with each **major source**
(a BioProject contributing at least 10% of that drug's rows in that
population) excluded in turn, with the same paired whole-lineage bootstrap.

**Status and provenance, stated plainly.** This is a mandatory reported
sensitivity analysis, not a credit condition, and it was added **after** a
development diagnosis: the human-clinical gentamicin Evo 2 gain came entirely
from one phenotype-enriched collection (PRJNA1297298; excluding it, AMR and
AMR+Evo 2 both reach 0.965 essential agreement). The Evo 2 credit rule in
`FOUNDATION_MODEL_PROTOCOL.md` is left unchanged so that it is not tightened
after its outcome became predictable; the single-laboratory external
evaluation remains decisive. The primary cross-validation design is also
unchanged.
