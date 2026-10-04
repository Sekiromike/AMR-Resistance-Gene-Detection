# Country Derivation and Intended-Use Population Amendment

Amendment: 006
Date: 2026-09-23
Status: ADOPTED (approved by the decision owner 2026-09-23)
Decision owner: Rushath Rajeev
Depends on: 001 (MIC endpoint), 005 (deduplication)
Resolves: blockers 8, 9, 10 and 12 in [`COHORT_REBUILD_READINESS.md`](COHORT_REBUILD_READINESS.md)

## Decision 1 — country from the INSDC `geo_loc_name` prefix

`country` was blank for all 25,735 development AST rows, and it is a required
cohort field, so every row would have failed validation. The August source
audit had ruled that `geo_loc_name` is free text and never fills `country`.

That ruling is revised. INSDC defines `/geo_loc_name` as
`<name>[:<region>][, <locality>]`, where `<name>` **must** come from a
controlled vocabulary. The prefix is therefore a structured field, not free
text, and reading it is extraction, not inference.

**Rule.** An explicit `country` always wins. Otherwise, `country` is the
`geo_loc_name` prefix before the first colon, **only if** that prefix is in the
frozen INSDC vocabulary. Anything else, including missingness tokens, leaves
`country` absent. Every row records its `country_source`: `explicit`,
`insdc_geo_loc_name_prefix`, or `absent`.

**Frozen vocabulary.** `config/insdc_geo_loc_name_vocabulary.tsv`, tracked in
git: 280 current and 17 historical names, from the INSDC page revised
July 08, 2024.

| Field | Value |
|---|---|
| Source | https://www.insdc.org/submitting-standards/geo_loc_name-qualifier-vocabulary/ |
| Retrieved | 2026-09-23 |
| Page SHA-256 | `c4760a193eaedbaf4de0076990b0aa65ab7b848d3dca915f1ed2e718f30b20be` |
| Vocabulary SHA-256 | `f8a9e6fa7338bf5241ec168bffd44d0a7fffd5a123b37c47d872d1e14db1bfcf` |

**Result on the frozen snapshot.** The data contains 31 distinct prefixes. All
30 real country names are valid INSDC names; the only mismatch is
`not collected` (17 isolates), which correctly stays absent. Country is filled
for 25,403 of 25,735 rows and absent for 332.

**Proof that nothing else changed.** The regenerated normalized table's hash
moved from `e72661c6…` to `653c89a7…` because three columns were added
(`mic_eligible`, `mic_ineligibility_reasons`, `country_source`) and countries
were filled. Regenerating with country derivation off and those columns
dropped reproduces the audited `e72661c63ecd96b2…` byte for byte.

## Decision 2 — the primary population is human clinical

The development set mixes human clinical isolates with food, animal and
environmental isolates (ground turkey 930, dogs about 1,361, chickens 1,127).
The external set is entirely human clinical. Training on retail meat and
testing on hospital patients would answer a muddier question than the study
intends.

**Rule.** The primary analysis is restricted to `human_clinical` isolates. The
all-source cohort is a prespecified **sensitivity analysis**, reported
alongside and never substituted for the primary.

**Classification (development).** It uses only source metadata, never a
phenotype, and never guesses:

| NCBI `epi_type` | `host` | Population |
|---|---|---|
| clinical | Homo sapiens | `human_clinical` |
| clinical | absent | `undetermined` |
| clinical | non-human | `non_human_or_environmental` |
| environmental/other | Homo sapiens | `undetermined` (conflict) |
| environmental/other | other or absent | `non_human_or_environmental` |
| absent | any | `undetermined` |

NCBI's `epi_type = clinical` alone is **not** sufficient: 259 clinical isolates
in the snapshot have animal hosts (veterinary clinical). A human host is
required. The snapshot has exactly one human spelling, `Homo sapiens`, and
every isolate's host, epi_type and isolation source agree across its rows.

**Classification (external).** JARBS isolates are `human_clinical` by
documented study design (hospital clinical-isolate surveillance). That basis is
recorded per isolate; it is not inferred from the specimen.

## Resulting cohort

| Population | Isolates | Ceftriaxone MICs | Ciprofloxacin MICs | Gentamicin MICs |
|---|---:|---:|---:|---:|
| **Human clinical (primary)** | **3,935** (3,888 with an MIC) | **2,835** | **3,853** | **3,376** |
| Non-human or environmental | 6,167 | 4,131 | 4,542 | 6,033 |
| Undetermined | 232 | 232 | 232 | 232 |

The primary development cohort remains larger than the best published
genome-linked *E. coli* MIC study (2,875 isolates).

## Consequences that must be reported

1. **Geography is narrow.** Human-clinical development isolates are 97% USA and
   UK (3,135 and 688). Any "leave a country out" evaluation is in practice USA
   versus UK, and is reported as such.
2. **Food and animal isolates are not wasted.** They enter the all-source
   sensitivity analysis, and whether they help or hurt human-clinical
   prediction is itself a reportable result.
3. **29 external isolates lack a specimen source** and are excluded by the
   specimen requirement. The count within the locked external set is reported
   in the cohort flow.

## Implementation

| Piece | Where | Tests |
|---|---|---|
| Country derivation | `normalize_ncbi_ast.py` `derive_country`, `load_country_vocabulary`; vocabulary is a tracked Snakemake input | prefix fill, missingness token, unknown prefix, explicit wins, historical names, vocabulary validation |
| Population and specimen source | new `build_source_attributes.py` | one test per classification row, cross-row conflict, cross-cohort identifier |
| Carried to cohort | `build_isolate_metadata.py`, `construct_research_cohort.py`, `validate_research_cohort.py` (allowed values enforced) | fixtures updated; tiny end-to-end workflow passes 8/8 |

Not yet implemented: restricting model fitting and evaluation to
`human_clinical`. That belongs with the MIC-endpoint evaluation layer (readiness
blocker 11), which does not exist yet.

## Implementation note (2026-09-27): absent values in the cohort

The first complete MIC cohort (build `64939994`) kept 59 development rows with
country absent in source (`not collected` or empty) and 3 with collection date
absent. Consistent with Decision 1, they are not filled. Under the MIC endpoint
the validator reports them as warnings: they remain in the lineage-grouped
evaluation and are excluded from the geographic and forward-time holdouts,
which need those fields. The categorical endpoint still requires both.

Collection dates are validated at the precision the submitter gave (INSDC
`YYYY`, `YYYY-MM` or `YYYY-MM-DD`); a missing month or day is never filled in.
