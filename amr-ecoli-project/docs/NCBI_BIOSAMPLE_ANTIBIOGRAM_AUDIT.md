# NCBI BioSample Antibiogram Overlay Audit

Checkpoint date: 2026-08-24 (America/New_York)

## Evidence boundary

This is an immutable source/provenance overlay for the frozen NCBI Pathogen
Detection AST snapshot. It is not an analysis-ready phenotype cohort, an
external validation set, or evidence of model performance. Missing method
detail, testing-standard version, AST testing date, disk potency, clinical
indication, and laboratory site remain missing.

## Acquisition identity and integrity

The live EFetch acquisition completed in
`data/source/ncbi_biosample_antibiograms/20260825T020258Z-v1`. The request list
was derived only from the hash-verified parent isolate file; the mutable NCBI
AST table was not reselected.

| Artifact | Rows or batches | SHA-256 |
|---|---:|---|
| Acquisition manifest | 1 | `38a76b8b44305b51504b86f9039f0578d554d4b1825dfc3afa22e8f20dc4fe8d` |
| Request plan | 104 batches | `a41dcece1a39ff4660cae04e42d616f119edfcc045a42e6da2f224b4b262ac71` |
| Antibiogram overlay | 165,462 | `ae9556c78c20f18e1025f684c42b545897d6d2fc7ed6011b0351b2b566f2b93a` |
| AST reconciliation ledger | 25,735 | `64baa271fe86e346aaa5ccf2dd377a40264d878f793f0ac7ea1b5e8db4316120` |

All 104 raw XML batches were independently rehashed after installation: zero
files were missing or mismatched. Their combined size is 100,795,808 bytes.
The complete installed directory contains 108 files and 217,784,646 bytes.
All 10,334 expected BioSample accessions appeared exactly once, every record
contained an antibiogram, and no unknown antibiogram headers were observed.
The runtime contact email and API key were not persisted; an exact contact scan
over every installed artifact returned no matches.

The first attempt used 200-accession batches and NCBI interrupted a chunked XML
response after 1,018,739 bytes. Atomic installation left no partial dataset.
The transport failure was made retryable, covered by a regression test, and the
successful acquisition used 100-accession batches with bounded retries.

## Coverage for the three prespecified drugs

The overlay contains exactly 25,735 rows for ciprofloxacin, ceftriaxone, and
gentamicin, matching the parent AST row count:

| Antibiotic | BioSample rows |
|---|---:|
| Ceftriaxone | 7,223 |
| Ciprofloxacin | 8,792 |
| Gentamicin | 9,720 |

Field coverage across those 25,735 target rows is:

| Submitted BioSample field | Present | Missing | Interpretation limit |
|---|---:|---:|---|
| Laboratory typing method | 25,735 | 0 | 25,466 `MIC`; 269 `disk diffusion`; modality alone is not a complete physical method |
| Laboratory typing platform | 22,875 | 2,860 | Platform is not laboratory site or AST date |
| Vendor | 19,231 | 6,504 | Values such as `NA` remain raw and are not promoted to present metadata |
| Method version or reagent | 21,393 | 4,342 | A panel/reagent is not a testing-standard version |
| Testing standard | 25,713 | 22 | Authority name alone does not identify a version |
| Testing-standard version | 0 | 25,735 | Not provided by the BioSample antibiogram schema |
| AST testing date | 0 | 25,735 | Record and collection timestamps are not substituted |

The 269 disk-diffusion rows still have zero explicit disk-potency values. Their
method-version/reagent cells are blank or `missing` for 261 rows; the remaining
eight describe Mueller-Hinton agar rather than disk content.

## Cross-source reconciliation

Every parent AST record had at least one exact BioSample match on BioSample
accession, antibiotic, submitted phenotype, comparator, decimal-equivalent
measurement, measurement unit, and MIC/disk modality.

- 25,725 parent records had one unique core match.
- 10 parent records remained ambiguous, forming five isolate-drug groups in
  two BioSamples with two identical core candidates apiece.
- The broader parent audit still contains 11 repeated isolate-drug groups.
  Matching a source row does not resolve which repeat belongs in a one-row-per-
  isolate-drug cohort; a prespecified reconciliation policy is still required.
- Among unique matches, submitted standard agreed when present for 25,703 rows,
  platform for 22,865, vendor for 19,221, and reagent for 21,393. Missing parent
  values were not treated as agreements.

## Scientific conclusion

BioSample closes the accession-complete provenance-overlay gate and recovers
useful platform, vendor, and reagent context. It does not close the reference-
standard gate. There is still no explicit AST testing date or testing-standard
version, the method label remains generic, disk potency is absent, and clinical
indication and laboratory site are unresolved. Consequently, zero rows may be
promoted to breakpoint-eligible phenotype evidence from this overlay alone.

The next cohort milestone is to freeze and audit the JARBS supplement and
sequence manifests, prove development/external disjointness, independently
review the EUCAST transcription, and resolve all repeated isolate-drug groups
before any modeling.
