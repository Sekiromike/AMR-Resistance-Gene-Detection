# NCBI Source Normalization Audit

## Scope and evidence boundary

This audit records row-complete selected-field normalization of the frozen NCBI
Pathogen Detection AST snapshot. The immutable source remains authoritative for
fields outside the normalized schema. This is source-cohort evidence only, not phenotype
validation, a frozen analysis cohort, external validation, or evidence that the
submitted AST categories are clinically correct.

The normalizer preserves submitted categories and quantitative measurements,
but does not infer assay method, AST testing date, disk potency, clinical
indication, laboratory site, country, or breakpoint edition. A narrow
`breakpoint_eligible` flag checks only source-stage assay metadata: method and
testing date for all measurements, plus disk potency for zone measurements. It
does not establish overall cohort readiness.

## Frozen input provenance

| Artifact | Rows | SHA-256 |
|---|---:|---|
| `data/source/ncbi_ast_snapshot/ncbi_ast.csv` | 25,735 | `988c2ff295c19595f2e65ba898bb86890130d8d1031b64dc4234b950a494c4bd` |
| `data/source/ncbi_ast_snapshot/ncbi_isolates.csv` | 10,334 | `b62f023abe872d082205a5b531d7928501c788bf6424d794ddf8ee003880479c` |
| `data/source/ncbi_ast_snapshot/ncbi_acquisition_manifest.json` | - | `ab65e3b33943f6ae99ab7b664a5403dab577d877617a91e0987d7ec9e4044d16` |
| `data/source/ncbi_ast_snapshot/ncbi_acquisition_manifest_amendment.json` | - | `0989b7c8b56403a7ddf874f17a1ac008671b98d3d3dc23d67834870ceb5d76a0` |

The original acquisition manifest is unchanged. Its append-only amendment
corrects two provenance representations: the executed comma-separated
antibiotic parameter was originally recorded with pipe separators, and the
client identity is now represented by the portable basename rather than a
user-specific absolute path. The query text, source-table metadata, source
exports, row counts, and their hashes are unchanged.

Normalization schema `2.0.0` independently verified the two input hashes, both
query hashes, and the source-table-metadata hash before processing.

## Structural normalization result

| Result | Count |
|---|---:|
| Source AST rows | 25,735 |
| Structurally normalized rows | 25,735 |
| Structural exclusions | 0 |
| Native source IDs preserved | 25,735 |
| Source checksums preserved | 25,735 |
| Source-stage breakpoint-eligible rows | 0 |

Derived artifacts:

| Artifact | Rows | SHA-256 |
|---|---:|---|
| `data/interim/ncbi_ast_normalized.csv` | 25,735 | `e72661c63ecd96b29313161632b32dc89ad4b1c937f91e5cc0d5aec7d66ff3cf` |
| `data/interim/ncbi_ast_exclusions.csv` | 0 | `aabbcfeb213c9a9256d125d870d0633bc4bf50fde0602888acf02b8ce11128a8` |
| `data/interim/ncbi_ast_manifest.json` | - | `6109fd45930c6b385c708204f1e737b0bbca528445b11c322a76ac1595a1c6c0` |

`data/interim/ncbi_ast_manifest.json` authenticates both outputs and records the
exact input/output schemas, byte counts, field missingness, eligibility reason
counts, mapping version, and acquisition-manifest references.

## Endpoint and measurement profile

| Antibiotic | Rows | Isolates | MIC | Zone |
|---|---:|---:|---:|---:|
| Ceftriaxone | 7,223 | 7,219 | 7,198 | 25 |
| Ciprofloxacin | 8,792 | 8,786 | 8,627 | 165 |
| Gentamicin | 9,720 | 9,714 | 9,641 | 79 |
| **Total** | **25,735** | **10,334 unique overall** | **25,466** | **269** |

Comparator semantics are retained. Of 25,466 MIC rows, 18,136 (71.22%) are
censored (`<`, `<=`, `>`, or `>=`) and 7,330 are exact. Of 269 zone rows, 48
are censored and 221 are exact. Censored MICs must not be treated as point
measurements.

Submitted categories are retained separately from any future independent
interpretation: 18,423 `S`, 4,263 `R`, 273 `I`, 17 `NS`, and 2,759
`NOT DEFINED`. Submitted standard names are 23,182 `CLSI`, 2,481 `EUCAST`, 50
`SFM`, and 22 blank. A standard name does not identify a breakpoint edition.

## Missingness and unresolved semantics

| Source-stage issue | Affected AST rows |
|---|---:|
| AST method absent | 25,735 |
| AST testing date absent | 25,735 |
| Disk potency absent | 25,735 total; required for all 269 zone rows |
| Explicit normalized country absent | 25,735 |
| Assembly accession absent | 1,030 |
| Collection date absent | 341 |
| Raw `geo_loc_name` absent | 297 |

The raw `geo_loc_name` field is preserved and is not promoted to country or
laboratory site. Explicit row-level clinical indication and AST laboratory site
are not present in the frozen source schemas. `epi_type`, specimen/isolation
source, host, BioProject, platform, vendor, and reagent are preserved as source
metadata but cannot substitute for assay method, indication, or site.

The AST table's `creation_date` and isolate table's `target_creation_date` are
preserved separately as `ast_record_creation_date` and `target_creation_date`.
They disagree for all 10,334 targets by exactly four or five hours, including
1,674 targets whose calendar date changes. Their upstream timezone/source
semantics require clarification; neither is used as collection or AST testing
date.

The 1,030 rows without an assembly represent 348 isolates. All 25,735 AST rows
otherwise link to exactly one of 10,334 unique frozen isolate targets and
BioSamples, with no observed BioSample, BioProject, collection-date, or raw
geography disagreement between the two exports.

## Versioned isolate-enrichment overlay

An authenticated BigQuery `FOR SYSTEM_TIME AS OF` query was run against the
exact 10,334 sorted target accessions from the frozen isolate export. It did not
reselect AST records. The query used the parent retrieval timestamp
`2026-08-21T01:26:04.521087Z`, an explicit 10,000,000,000-byte safety cap, and
the confirmed unbilled sandbox project. BigQuery reported 909,529,512 bytes
processed and 910,163,968 bytes billed against sandbox quota; no billing
account was attached.

| Overlay artifact | Rows | SHA-256 |
|---|---:|---|
| `ncbi_isolate_enrichment.csv` | 10,334 | `90881fbc1b61b8190692d48fdae607afba666f6f834a3f899a0153fa5c12cbdf` |
| `ncbi_isolate_enrichment_query.sql` | - | `d7d8cc05a3f17ed289c0a0a98cf37af3eaf0763dd974eb8610b3cd50e358ed99` |
| `ncbi_isolate_enrichment_job_metadata.json` | - | `e5d8b3348e875bce1abc81a2ff49b7cf6d19fbc4a811f93526a08dbe90eaf6a9` |
| `ncbi_isolate_enrichment_source_metadata.json` | - | `44c34b4f7e46120c589ce944c8dd7b1c1c3450bf5d8e02e50122f19e8127c033` |
| `ncbi_isolate_enrichment_manifest.json` | - | `aa330009765bcafb8425f49b81f8a4a7d0e622a9cb5dc122e9c08caef2473bfd` |

All 10,334 source matches were one-to-one, with exact target, BioSample,
scientific-name, taxonomy, BioProject, collection-date, and raw-geography
reconciliation. The 9,986 existing assembly accessions were unchanged. The
source still had no assembly accession for 348 isolates, but every one of those
348 had an SRA run accession. Thus every frozen isolate has a sequence-retrieval
route (9,986 assemblies plus 348 read-only fallbacks), subject to subsequent
download, checksum, species, contamination, and assembly-QC gates.

The overlay also records raw collection organization, SRA center, host-disease,
epidemiological type, submitted assembly method/statistics, and sequence
platform where present. These fields remain provenance or triage context. They
do not establish AST method, AST testing date, laboratory/site, clinical
indication, surveillance membership, genome QC, lineage, or evaluation split.

## Repeats and reconciliation gate

There are 11 repeated isolate-antibiotic groups containing 27 records (maximum
three records per group). These are not exact duplicate rows. Some differ by
quantitative value or measurement type, five groups mix MIC and zone records,
and one group has discordant submitted categories. No record was collapsed.
A prespecified, audited repeat-reconciliation policy is required before cohort
freeze.

## Scientific decision

The snapshot is suitable as an immutable, quantitatively characterized source
ledger. It is not yet suitable as the publication phenotype cohort because no
row passes the source-stage assay-metadata gate and clinical indication, AST
site, repeat reconciliation, independently reviewed breakpoint rules, and an
untouched external source remain unresolved.

The next defensible work is to validate whether direct BioSample/NARMS records
recover source-level AST method metadata, identify a source that supplies the
remaining phenotype provenance, review the EUCAST transcription independently,
and prespecify repeat handling and the external cohort before any modeling.
