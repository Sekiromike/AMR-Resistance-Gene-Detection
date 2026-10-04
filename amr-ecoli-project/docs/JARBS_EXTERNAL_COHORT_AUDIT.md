# JARBS External-Cohort Candidate Audit

Checkpoint date: 2026-08-24 (America/New_York)

## Decision

JARBS is now an immutable, accession-linked **conditional external challenge
candidate** for the three-drug *E. coli* task. It is not yet a locked external
cohort. Exact identifier checks passed for BioSample, SRA run, registered
assembly, and public isolate identifiers, but raw-read identity and
whole-genome near-neighbor disjointness remain untested. No model may inspect
JARBS outcomes while those gates remain open.

This is source/cohort evidence only. It is not model-performance, clinical,
causal, replication, or laboratory-validation evidence.

## Source identity, license, and integrity

The source is Kayama et al., *Nature Communications* 14, 8046 (2023),
[doi:10.1038/s41467-023-43516-4](https://doi.org/10.1038/s41467-023-43516-4),
PMCID `PMC10698200`, and BioProject `PRJDB10842`. The article and PMC Open
Access record identify the work as CC BY 4.0. The article's standard caveat for
separately credited third-party material still applies; no separate license or
credit line was found for Supplementary Data 6.

The publisher file host returned a JavaScript client-challenge document under
the `.xlsx` filename. It was rejected by file-signature and size checks. The
workbook was instead recovered from the official Europe PMC supplementary-file
API. The installed workbook is 1,493,480 bytes, has publisher-declared MD5
`00659b377a4830745048ba1643656276`, and SHA-256
`1e0ae22795d31efd496ccd80fb68b2a64a10bee9056d37cc45ba782006188bf8`.
The complete Europe PMC supplementary archive has SHA-256
`44ab0f47f19e189331901a888a652fbdf8cc12f5adc58837cee35e9c9fc9c98e`.

The January 2024 author correction changes a Supplementary Information figure
count and two article-text counts. It does not report a replacement of
Supplementary Data 6, but remains part of the source provenance:
[doi:10.1038/s41467-024-45349-1](https://doi.org/10.1038/s41467-024-45349-1).

All source and derived artifacts are frozen under
`data/source/jarbs/20260825T021722Z-v1`. The runtime NCBI contact email was not
persisted; an exact repository-artifact scan returned no matches.

## Workbook audit

Supplementary Data 6 contains 5,143 unique isolate rows plus one footnote row.
There are 3,159 *E. coli* isolates. The workbook explicitly contains isolate
ID, assembly-QC values, taxonomic assignment, quantitative MIC fields, MLST,
and resistance genes. It does not contain BioSample or SRA accessions,
collection date, AST testing date, hospital/site, specimen, clinical
indication, row-level AST method, or breakpoint-standard version.

Exactly 2,535 *E. coli* isolates have all three prespecified MIC fields;
the same 624 isolates are blank for ceftriaxone, ciprofloxacin, and gentamicin.
Missing values remain missing.

| Drug | Present | Missing | Censored | Comparator distribution |
|---|---:|---:|---:|---|
| Ceftriaxone (`CTRX`) | 2,535 | 624 | 2,041 | 613 `<=`; 494 exact; 1,428 `>` |
| Ciprofloxacin (`CPFX`) | 2,535 | 624 | 2,232 | 801 `<=`; 303 exact; 1,431 `>` |
| Gentamicin (`GM`) | 2,535 | 624 | 2,235 | 1,745 `<=`; 300 exact; 490 `>` |

The heavy left- and right-censoring must be retained as intervals. Values such
as `>64` or `<=0.25` are not exact MIC measurements and must not be silently
converted to their displayed numeric endpoint.

## Method and sampling boundary

The paper reports centralized broth microdilution at the National Institute of
Infectious Diseases using MicroScan WalkAway, NEG MIC 3.31E, and NEG MIC NF 1J
panels. It reports CLSI 2021 interpretation, except the piperacillin/tazobactam
change described for 2022. These are paper-level method facts, not row-level
AST dates or per-row standard-version fields.

The isolates were collected in 2019-2020 from 175 hospitals across 45 Japanese
prefectures after enrichment for third-generation-cephalosporin resistance in
*E. coli*/*K. pneumoniae* or reduced meropenem susceptibility in
Enterobacterales. JARBS therefore cannot estimate routine clinical prevalence.
SRA collection dates, geographic labels, and isolation-source values are
preserved as raw source metadata, but none is substituted for AST testing date,
hospital identity, laboratory site, or clinical indication.

## Sequence-accession reconciliation

The frozen BioProject query returned 5,387 runs representing 5,144 public
isolate identifiers and BioSamples: 5,145 Illumina runs and 242 Oxford Nanopore
runs. The current BioProject is larger than the paper's 5,143-isolate table, so
BioProject membership alone is not an isolate manifest.

For the target *E. coli* subset, all 3,159 Supplementary Data 6 isolate IDs map
exactly to one Illumina run and one BioSample. Across all species, 5,140 rows
map to exactly one Illumina run. The three non-target anomalies are retained,
not guessed away:

- `JBCDACJ-19-0094` and `JBEAAEE-19-0003` occur in the supplement but not under
  those exact SRA isolate names; SRA instead contains two older numeric-style
  *Klebsiella aerogenes* names (`23029-19-0094`, `40044-19-0003`). No alias
  equivalence is asserted without an authoritative crosswalk.
- `JBBDAGI-19-0084` is an *Acinetobacter baumannii* SRA isolate absent from
  Supplementary Data 6.
- `JBBDAAF-19-0029` has two Illumina runs and one Nanopore run. It is an
  *Enterobacter* isolate, not part of the *E. coli* target set.

NCBI Datasets currently reports ten assembly accessions under the BioProject,
representing nine strain names: seven GenBank accessions and three RefSeq
accessions. Two paired accessions represent one *E. coli* strain. This is a
registered-complete-assembly manifest only; it does not cover every draft
genome represented in the supplement.

## Development/external overlap audit

The external manifests were compared with all 10,334 rows in the frozen NCBI
development isolate and enrichment snapshots. Exact intersections were zero
for each tested identifier class:

| Exact identifier class | Overlaps |
|---|---:|
| BioSample accession | 0 |
| SRA run accession | 0 |
| Registered assembly accession | 0 |
| Public isolate identifier | 0 |

This does not prove genomic independence. The development overlay has no raw
read hashes, and neither cohort yet has the sequence sketches/cluster ledger
needed to exclude same-strain or near-neighbor genomes under the prespecified
threshold. Consequently the external-lock status is
`blocked_genomic_disjointness_pending`.

## Frozen audit artifacts

| Artifact | Rows | SHA-256 |
|---|---:|---|
| SRA run manifest | 5,387 | `96d266eb48750bf6b941f475d66fce143a4e9507001803eff7a4abba7cce7d50` |
| Isolate-sequence manifest | 5,143 | `992a8a37c11d0b983fb90357b4ab61840892b11406f06f11a35170dc7f3227f1` |
| *E. coli* three-drug audit table | 3,159 | `fff5c979e8cc1014b077c4c4abab0727242e9a8a364163573568943003bf4b45` |
| Registered assembly manifest | 10 | `47348ba56b7472099b11f7c269c6f500c2c49ee2af37539c4df302aac2a3697b` |
| Exact-overlap ledger | 0 | `62483db0a7f4dc858c5f3b11e9cb711f7ffd1c78e813d487073da0b69e23ed83` |
| Machine-readable audit | 1 | `6d3732720d816a2cb9e96941d4f42325498b4ca22cba0aa73e339472e7143f3d` |

`scripts/audit_jarbs_external.py` recreates the derived tables from the frozen
OOXML, SRA XML/RunInfo, NCBI Datasets report, and development manifests. The
workbook reader is read-only and uses the Python standard library because the
workspace spreadsheet runtime was unavailable; the source workbook itself was
never rewritten.

## Required next gates

1. Build compatible sequence sketches or genome clusters for development and
   all 3,159 JARBS *E. coli* isolates, freeze the threshold, and exclude exact
   genomes and near neighbors before external designation.
2. Decide whether absent row-level AST testing date, hospital/site, and clinical
   indication are acceptable for a challenge-only endpoint. Do not invent
   replacements.
3. Independently review the EUCAST v16.1 transcription and predeclare the
   handling of interval-censored MICs and intermediate categories.
4. Lock the 2,535 three-drug phenotype rows and all exclusions before any model
   run can access JARBS outcomes.
