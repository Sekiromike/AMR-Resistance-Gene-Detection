# Third-party resources

The MIT license covers original project code and associated software
documentation. It does not relicense upstream datasets, database contents,
model weights, dependencies, or third-party publication material.

- **JARBS:** Kayama, Yahara, Sugawara and colleagues (2023), *National genomic
  surveillance integrating standardized quantitative susceptibility testing
  clarifies antimicrobial resistance in Enterobacterales*, Nature Communications
  14, 8046, [DOI](https://doi.org/10.1038/s41467-023-43516-4).
  The article and accompanying material carry
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), subject to any
  separately marked third-party exceptions. Our derived tables select E. coli,
  normalize identifiers/MIC records, and apply the documented exclusions.
  Source records and attribution are retained. Read accessions belong to
  BioProject PRJDB10842.
- **NCBI:** Pathogen Detection, BioSample, SRA, and assembly resources retain
  their upstream terms. See [NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/).
  NCBI places no restrictions on use/distribution of molecular database data,
  while noting that submitters may retain rights. Original accession and source
  identifiers are preserved in the derived research tables.
- **AMRFinderPlus:** annotations were generated with 4.2.7 and database release
  2026-08-07.1. The full database and executable are not redistributed here;
  use the [official repository](https://github.com/ncbi/amr) and upstream terms.
- **Evo 2:** the analysis used frozen `evo2_20b` representations through Evo 2
  0.5.5. We distribute derived embeddings and selected public bacterial sequence
  units, not pretrained weights. Obtain weights and implementation from the
  [official repository](https://github.com/ArcInstitute/evo2) under its terms.
- **OpenGenome2:** the pretraining overlap audit reports matches against the
  IMG/PR component. The large training corpus is not included in the release.
- **Other software:** pinned environment exports identify NumPy, SciPy, pandas,
  scikit-learn, XGBoost, PyTorch and other dependencies. Their licenses remain
  unchanged. Rendering tools and TeX distributions are not bundled.

No EUCAST workbook or rule-table content is redistributed in this release.
The primary endpoint is quantitative MIC; secondary S/I/R validation remains
uncompleted.

Third-party source snapshots were preserved locally. The release includes the
analysis inputs and acquisition metadata needed for score/model reproduction;
large public read collections and peripheral source records are retrieved
through the recorded upstream accessions and manifests.
