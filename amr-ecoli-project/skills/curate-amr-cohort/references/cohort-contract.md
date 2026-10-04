# Cohort Contract

The canonical table has one reconciled row per isolate and antibiotic. Required
fields are enforced in `scripts/validate_research_cohort.py` under one of two
endpoint profiles. `docs/ENDPOINT_AMENDMENT.md` records why.

## Endpoint Profiles

| Profile | Endpoint | Status |
|---|---|---|
| `mic_regression` | interval-censored MIC | primary |
| `categorical_sir` | EUCAST S/I/R | secondary, blocked |

`mic_regression` requires isolate identity, genome provenance, the quantitative
measurement with its comparator and unit, and the full evaluation structure. It
does not require assay method, AST testing date, clinical indication,
laboratory site, or breakpoint authority/version, because no measurement is
interpreted under it.

`categorical_sir` requires everything `mic_regression` requires plus every
interpretation field listed below. It still fails closed. An unspecified
profile resolves to `categorical_sir`.

Dropping a requirement from the primary profile never licenses inventing the
field. Absent values stay absent and stay recorded, and the validator emits a
warning naming every separation check it could not run.

## Required Fields

- isolate, BioSample, and versioned assembly accessions;
- source dataset, release, and AST record identifier;
- antibiotic, clinical indication, MIC or zone type, comparator, value, unit,
  method, and disk content where applicable;
- submitted category and independently derived category;
- breakpoint authority, version, and frozen artifact hash;
- collection date, site, country, and BioProject;
- lineage, genomic near-neighbor cluster, and prespecified evaluation split.

## Non-Negotiable Semantics

- Preserve `<`, `<=`, `>`, and `>=`; censored MICs are not point values.
- Do not infer breakpoint version from `CLSI` or `EUCAST` alone.
- Do not call BioProject, country, year, MLST, or submitter a laboratory site.
- Do not infer a clinical indication from specimen source. Apply indication-
  specific breakpoint rows only when indication is explicitly available.
- Do not treat submitter categories as independent reinterpretations.
- Do not collapse `I`, `SDD`, `NS`, or discordant repeats without a written rule.
- Do not allow a genomic near-neighbor cluster across development and external.

Primary sources:

- https://www.ncbi.nlm.nih.gov/pathogens/docs/ast/
- https://www.ncbi.nlm.nih.gov/pathogens/docs/ast_gcp/
- https://www.eucast.org/bacteria/clinical-breakpoints-and-interpretation/clinical-breakpoint-tables/
- https://clsi.org/shop/standards/m100/
