# Phenotype Source Candidate Audit

## Decision boundary

This is a prespecification aid, not an external-cohort freeze. No candidate may
enter model development or evaluation until its raw artifacts, license or terms,
identifiers, quantitative AST fields, method/version provenance, genome links,
selection mechanism, and overlap exclusions are frozen by hash.

No public source reviewed so far supplies all desired fields together:
genome-linked quantitative AST, row-level assay method, AST testing date,
clinical indication, and laboratory/site. Missing fields must remain missing;
collection date, specimen source, submitting center, BioProject, geography, and
host disease are not substitutes.

## Prespecified roles

| Candidate | Intended role | Strength | Disqualifying or limiting gap |
|---|---|---|---|
| FDA/USDA NARMS plus NCBI BioSample | Development-source provenance enrichment | Public isolate-level quantitative MIC data and genome identifiers; all three target drugs are part of the *E. coli* panel | Overlaps the current NCBI source; spreadsheet-level method and AST-test date are absent or unproven; food/animal source is not clinical indication |
| Japan JARBS-GNR (`PRJDB10842`) | Conditional untouched clinical external challenge | Centralized quantitative MIC remeasurement, explicit broth microdilution platform/panels, three target drugs, multi-hospital clinical surveillance, public sequences | Deliberately enriched for third-generation-cephalosporin resistance or reduced meropenem susceptibility; cannot estimate routine prevalence; exact public site/specimen/AST-date fields and overlap must be audited |
| CDC/FDA AR Isolate Bank MuGSI-ESBL | Small secondary stress/replication set | Curated BioSample links, specimen source, quantitative modal MIC, explicit broth microdilution | Small, multi-species, deliberately ESBL/cephalosporin-resistant panel; interpretation display can change over time; not a primary performance cohort |
| Vet-LIRN/NAHLN animal-pathogen AMR | Optional cross-host out-of-distribution analysis | Large veterinary surveillance resource with AST and a substantial sequenced subset | Laboratory platforms vary and complete three-drug genome-linked quantitative coverage is not yet established |

CDC BEAM, public CIPARS dashboards/downloads, and FSIS summary reports are not
eligible isolate-level phenotype sources for this project because their public
forms are aggregate rather than genome-linked quantitative AST records.

## Evidence supporting the roles

The FDA describes NARMS as integrated surveillance of human clinical, animal
slaughter, and retail-meat samples, provides downloadable data and a data
dictionary, states that the tables are non-confidential, and requests NARMS
attribution. The CDC's current *E. coli* table includes gentamicin,
ceftriaxone, and ciprofloxacin as MIC-tested drugs. These properties make NARMS
a strong source-provenance target, but not an independent external cohort when
the same BioSamples already occur in the frozen NCBI snapshot.

- [FDA NARMS data downloads and data dictionary](https://www.fda.gov/animal-veterinary/national-antimicrobial-resistance-monitoring-system/integrated-reportssummaries)
- [CDC NARMS antibiotics tested](https://www.cdc.gov/narms/about/antibiotics-tested.html)
- [NCBI BioSample antibiogram package](https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/)
- [FDA NARMS BioProject PRJNA290865](https://www.ncbi.nlm.nih.gov/bioproject/290865)
- [CVM NARMS *E. coli* BioProject PRJNA292663](https://www.ncbi.nlm.nih.gov/bioproject/292663)

A live Entrez pilot on [BioSample SAMN12912155](https://www.ncbi.nlm.nih.gov/biosample/SAMN12912155)
confirmed an antibiogram table containing all three target drugs, comparator,
MIC, unit, submitted category, platform `Sensititre`, vendor `Trek`, reagent
`GM-NEG`, and standard `CLSI`. Its method cell is only `MIC`, not a complete
physical assay-method description, and the table has no AST testing-date
column. This one accession proves recoverability of additional submitted
metadata, not its cohort-wide completeness or fitness for breakpoint use.

The JARBS study sequenced 5,143 Enterobacterales isolates and centrally
remeasured quantitative MICs for 4,195 using broth microdilution on MicroScan
WalkAway with NEG MIC 3.31E and NEG MIC NF 1J panels. The measured panel
included ceftriaxone, ciprofloxacin, and gentamicin. The sequenced collection
included 3,159 *E. coli* and drew from 175 hospitals across 45 Japanese
prefectures in 2019-2020. Its selection was resistance-enriched, which is why it
is appropriate only as a difficult external challenge with that sampling frame
reported explicitly.

Version `20260825T021722Z-v1` now freezes Supplementary Data 6, the official
Europe PMC archive, 5,387 current SRA run records, and the registered-assembly
report. The workbook contains 3,159 *E. coli* rows; 2,535 have all three target
MIC fields and 624 have none of the three. All 3,159 *E. coli* isolate IDs map
to exactly one Illumina run and BioSample. Exact overlap with the 10,334-isolate
development snapshot is zero for BioSample, SRA run, registered assembly, and
public isolate identifiers. This is not yet a disjointness proof: development
read hashes are unavailable and genome-level near-neighbor clustering remains
pending. See [`JARBS_EXTERNAL_COHORT_AUDIT.md`](JARBS_EXTERNAL_COHORT_AUDIT.md).

- [JARBS primary article and supplementary files](https://www.nature.com/articles/s41467-023-43516-4)
- [JARBS Supplementary Data 6](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-43516-4/MediaObjects/41467_2023_43516_MOESM9_ESM.xlsx)
- [JARBS BioProject PRJDB10842](https://www.ncbi.nlm.nih.gov/bioproject/1262179)

The CDC/FDA MuGSI-ESBL panel documents a resistance-enriched panel, BioSample
accessions, collection year, country, urine/blood source, and modal MICs from
broth microdilution. It is suitable for stress testing only.

- [CDC/FDA AR Isolate Bank MuGSI-ESBL panel](https://wwwn.cdc.gov/arisolatebank/Panel/PanelDetail?ID=1162)
- [AR Isolate Bank quality and method documentation](https://wwwn.cdc.gov/arisolatebank/QA)

## BioSample overlay implementation checkpoint

`scripts/acquire_ncbi_biosample_antibiograms.py` verified the frozen AST,
isolate, query, metadata, manifest, and amendment hashes before completing live
version `20260825T020258Z-v1`. The retrieval fixed all 10,334 BioSample
accessions in 104 POST batches of 100, archived 165,462 antibiogram rows, failed
closed against identity drift, and wrote one reconciliation record for every
frozen native AST record. Of those records, 25,725 had one unique core match and
10 remained ambiguous across five isolate-drug groups. See
[`NCBI_BIOSAMPLE_ANTIBIOGRAM_AUDIT.md`](NCBI_BIOSAMPLE_ANTIBIOGRAM_AUDIT.md).

The parser preserves unknown headers and includes dedicated
`testing_standard_version_raw` and `ast_testing_date_raw` fields even though the
documented BioSample antibiogram schema does not provide those columns. This is
intentional: record submission/update dates and isolate collection dates cannot
fill AST testing date, and `CLSI` or `EUCAST` alone cannot fill a version.

## Next acquisition gates

1. Before designating JARBS external, exclude any shared BioSample, SRA,
   assembly, read hash, patient/outbreak component, or genomic near-neighbor
   from development. If disjointness cannot be proven, redesign the development
   source or reject JARBS as external.
2. Resolve whether a challenge-only cohort may proceed without public
   isolate-level AST testing date, hospital/site, and clinical indication.
   Missing fields remain missing.
3. Freeze the external eligibility rule, censoring policy, and 2,535-row
   three-drug endpoint set before inspecting model performance.

Until those gates pass, the project remains at source/cohort evidence and must
not report phenotype model performance.
