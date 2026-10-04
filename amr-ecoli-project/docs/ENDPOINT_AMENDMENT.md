# Primary Endpoint Amendment

Amendment: 001
Date: 2026-09-20
Status: ADOPTED
Decision owner: Rushath Rajeev
Supersedes: the categorical-only endpoint implied by `config/study.json`
schema 1.0.0 and enforced by `scripts/validate_research_cohort.py`.

## Decision

The primary endpoint changes from a categorical EUCAST S/I/R call to an
**interval-censored minimum inhibitory concentration (MIC)** per isolate and
antibiotic. The categorical S/I/R endpoint is retained as a **prespecified
secondary endpoint** and remains blocked until independent EUCAST review and
source-level assay metadata exist.

This is an endpoint and evaluation-design change under
`skills/manage-amr-research/references/governance.md`. The `cohort` phase is
reopened accordingly.

## Why the previous endpoint cannot be reached

The categorical endpoint requires a defensible measurement-to-category
interpretation. The cohort contract therefore requires row-level AST method,
AST testing date, clinical indication, laboratory site, and breakpoint
authority *and version*. Frozen source audits in this repository show those
fields are absent, not merely incomplete:

| Required field | Development source (NCBI) | External candidate (JARBS) |
|---|---|---|
| AST method | absent 25,735 / 25,735 | broth microdilution, panel documented |
| AST testing date | absent 25,735 / 25,735 | not published at isolate level |
| Testing-standard version | 0 recovered by BioSample overlay | not published |
| Clinical indication | not in source schema | not published |
| Laboratory / hospital site | not in source schema | not published |
| Disk potency (zone rows) | absent for all 269 zone rows | not applicable |

Consequence, recorded in `docs/NCBI_SOURCE_NORMALIZATION_AUDIT.md`:
**0 of 25,735 source AST rows are breakpoint-eligible.** The two reason codes
responsible are `missing_ast_method` and `missing_ast_testing_date`.

No public genome-linked *E. coli* AST resource supplies these fields. Published
work that has them obtained them through direct institutional record access,
which is out of scope for a computational study built on public archives. The
categorical endpoint is therefore unsatisfiable from public data, and the
`cohort` phase could not exit under it.

`governance.md` states that a phase is blocked "only for a concrete external
dependency; record it and continue independent work." Independent EUCAST review
and absent source metadata are exactly such dependencies. They are recorded
here and below, and work continues on an endpoint that does not consume them.

## Why MIC is a valid primary endpoint here

MIC is the measured quantity. Predicting it requires no breakpoint table, no
breakpoint version, and no clinical indication, because no interpretive step is
performed. `docs/STUDY_PROTOCOL.md` already prefers this framing: "Prefer
quantitative MIC prediction with interval censoring when enough MIC data exist;
derive categorical performance as a downstream evaluation."

Available quantitative measurements in the frozen development snapshot:

| Antibiotic | MIC rows | Zone rows |
|---|---:|---:|
| Ceftriaxone | 7,198 | 25 |
| Ciprofloxacin | 8,627 | 165 |
| Gentamicin | 9,641 | 79 |
| **Total** | **25,466** | **269** |

Of 25,466 MIC rows, 18,136 (71.22%) are censored (`<`, `<=`, `>`, `>=`) and
7,330 are exact. Censoring is handled by interval-censored likelihood, not by
imputing point values.

## Primary analysis specification

- **Target:** log2 MIC, modelled as interval-censored. A row reported `<= x`
  contributes the interval (0, x]; a row reported `> x` contributes (x, inf).
- **Primary metric:** essential agreement, the proportion of predictions within
  one doubling dilution of the measured MIC. Exact agreement is reported
  alongside it.
- **Secondary metrics:** interval-censored log-likelihood, bias by dilution,
  and error stratified by distance to the nearest training genome.
- **Unit of analysis:** one reconciled row per isolate and antibiotic.
- **Evaluation structure:** unchanged. Lineage-grouped development folds,
  forward-time and geographic holdouts, and one locked external evaluation.

## Secondary endpoint, blocked

The categorical S/I/R endpoint remains specified and remains enforced by the
`categorical_sir` validator profile, which still fails closed. It is unblocked
only when both hold:

1. A second independent human reviewer completes the EUCAST v16.1
   transcription per `docs/EUCAST_V16_1_TRANSCRIPTION.md`, and an approved
   rules CSV exists with a recorded SHA-256.
2. Row-level assay method, testing date, and clinical indication are available
   from a documented source for the rows being interpreted.

Until then, no S/I/R category derived in this project is publication evidence.

## What this amendment does not change

- No integrity rule is relaxed. Absent fields stay absent and stay recorded.
- The `categorical_sir` profile keeps every field requirement it has today.
- All structural checks apply to both profiles: accession formats, species
  confirmation, genome QC, unit/measurement-type agreement, comparator
  validity, duplicate reconciliation, isolate-to-assembly uniqueness, and the
  prohibition on genomic near-neighbour clusters or deduplication groups
  crossing the development/external boundary.
- The external cohort stays sealed from feature, model, threshold, and epoch
  selection.
- Evidence tiers are unchanged.

## Limitation that must be reported

MICs in the development snapshot originate from heterogeneous laboratories with
**unlabelled assay methods and unknown dilution ranges**. This inflates
measurement error relative to a single-laboratory study and means observed
essential agreement is not directly comparable to a study using one documented
method. This limitation belongs in the abstract and the limitations section of
any output, not in a footnote.

Interval censoring absorbs part of the effect, because panel range limits are
expressed through the comparator. It does not remove it.

Absence of AST testing date also prevents any claim about breakpoint-era
consistency, and absence of site prevents leave-site-out evaluation. The
geographic and forward-time holdouts use `country` and `collection_date`, which
are source-derived and are not substitutes for site or AST date.

## Evidence tier path

`SMOKE_TEST_ONLY` (current)
→ `DEVELOPMENT_ONLY` after locked interpretable baselines on development folds
→ `LOCKED_EXTERNAL` after one frozen evaluation on the external cohort.

No clinical, causal, or diagnostic claim is authorised at any of these tiers.

## Affected artifacts

| Artifact | Change |
|---|---|
| `config/study.json` | adds `primary_endpoint`, `secondary_endpoint`, `mic` block; moves `clinical_indication` under the secondary endpoint |
| `scripts/validate_research_cohort.py` | required-column set and category checks become profile-conditional |
| `scripts/normalize_ncbi_ast.py` | adds `mic_eligible` and `mic_ineligibility_reasons` beside the unchanged `breakpoint_eligible` |
| `skills/curate-amr-cohort/references/cohort-contract.md` | documents both profiles |
| `docs/PROJECT_STATUS.json` | `cohort` phase reopened with this amendment recorded |

## Open items this amendment does not resolve

- Independent EUCAST v16.1 review (needs a second human).
- AMRFinderPlus database release is still `REQUIRED_BEFORE_COHORT_FREEZE`.
- JARBS remains a conditional external challenge until the cross-cohort
  collision ledger is resolved by the prespecified sequence-only rules.
- JARBS sampling is resistance-enriched; its prevalence does not estimate
  routine prevalence and must be reported as such.
