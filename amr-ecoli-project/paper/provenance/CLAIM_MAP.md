# Claim-to-source map

The CSV tables contain original machine-readable source file names. The immutable
ledger excerpt is `evidence_snapshot.json`; `source_manifest.json` records
SHA-256 hashes of the report and implementation files used by the build.

| Manuscript claim | Source | Evidence boundary |
|---|---|---|
| 29,619 rows, 11,520 isolates; 3,849 human-clinical development isolates | Ledger entry 2026-09-27T15:34:27.949526+00:00, cohort job 64948707 | Cohort file retrieved and counts recomputed in verification/audit.json |
| 485 development CV units; 106 external lineage groups | Ledger entry 2026-09-27T18:10:44.312763+00:00 and amendment 007 | Structural cohort analysis |
| Gentamicin gain +0.027/+0.024 in development | Ledger entry 2026-09-28T19:14:59.128254+00:00; bounded rerun reproduces it | Development paired comparison |
| Gain disappears when PRJNA1297298 rows are removed | Ledger entry 2026-09-29T02:22:07.702935+00:00; source-robustness script | Re-evaluation of existing predictions, not refitting |
| No demonstrated external Evo 2 gain; allele-feature EA by drug | `docs/EXTERNAL_RESULTS.md` | Locked external, score job 65020832 |
| Family features reach 0.990 on the 102-isolate subgroup | `docs/POSTHOC_FAMILY_ANALYSIS.md` | Post-hoc external subset and model design |
| Family XGBoost ceftriaxone EA 0.854, paired difference 0.257 | `docs/POSTHOC_FAMILY_ANALYSIS.md` | Post-hoc external comparison |
| Family beats Evo 2 by 0.046/0.111 on pooled withheld CTX-M observations | `docs/LEAVE_VARIANT_OUT_PLAN.md` | Controlled development; hypotheses fixed before this experiment |
| 568 pooled observations represent 565 unique isolates in 45 lineages | `scripts/leave_variant_out.py`, `resistant_side` and `paired_bootstrap` | Unit is target–isolate; bootstrap is by lineage |
| Evo 2-only features still use annotation-selected units | `scripts/extract_resistance_loci.py`, `embedding_features.py`, protocol revision 0.2 | Limits scope of representation comparison |
| 14,482/22,767 acquired units have near-identical pretraining-component matches | Ledger entry 2026-09-28T17:11:18.327358+00:00 | Partial sequence-overlap audit; no demonstrated MIC-label leakage |
| Approximately 3.6 GPU-hours for extraction | Ledger entry 2026-09-28T19:14:59.128254+00:00 | Extraction only; not full study cost |

## Verification of original results

`verification/audit.json` records score, interval, and hash verification.
`verification/all_evaluations.csv` identifies original evaluation files by run
path; `paired_comparisons.csv` and `source_sensitivity.csv` carry the complete
comparisons. The source reports in the table above preserve analysis chronology;
the rendered paper now uses the verified full-precision original results.
`release_manifest.json` retains original and distributed hashes for public copies.

## Literature verification

Primary publications and official software documentation were retrieved with
Exa on 30 September 2026. Five literature searches returned 25 candidate results; the
selected sources below were opened and checked. Search results and secondary
summaries were not treated as the final citation authority.
The separate arXiv category lookup reviewed three results and confirmed the
category names against the official taxonomy and quantitative-biology archive.

| Key | Verified page | Use in paper |
|---|---|---|
| `brixi2026` | https://www.nature.com/articles/s41586-026-10176-5 | Evo 2 foundation model; published Nature citation |
| `evo2software` | https://github.com/ArcInstitute/evo2 | Official 20B checkpoint and inference implementation |
| `feldgarden2021` | https://www.nature.com/articles/s41598-021-91456-0 | Curated gene, mutation, and functional-class representation |
| `kayama2023` | https://www.nature.com/articles/s41467-023-43516-4 | JARBS source identity and standardized susceptibility testing |
| `moradigaravand2018` | https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006258 | Prior genomic E. coli resistance prediction |
| `yu2025` | https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3003539 | Population structure and biased-sampling confounding |
| `varma2006` | https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-91 | Model selection and nested cross-validation |
| `xgboostaft` | https://xgboost.readthedocs.io/en/stable/tutorials/aft_survival_analysis.html | Censored regression objective and bounded labels |

The opening literature coverage is targeted, not a systematic review. The
draft makes no priority claim such as "first" and does not infer novelty from
an absence of search results.
