# Editorial notes for the first preprint

The user confirmed that the experiments are finished. The items below concern
accurate reporting, release of existing artifacts, and submission details.
They do not propose new experiments as a prerequisite for this draft.

## Main conclusion and its scope

The defensible headline is that curated allele, gene-family, and drug-class
features outperform the tested frozen Evo 2 embedding representation on
withheld resistance variants. The locked external evaluation also found no
demonstrated gain from adding Evo 2 to the allele-level model. Preserve both
findings; a universal claim about all uses of genomic foundation models would
go beyond the experiments.

The family model includes **allele + family + subclass** indicators. Calling
it a pure family-only model would misdescribe the implementation. Similarly,
"Evo 2 alone" means no binary catalogue features in the learner: AMRFinderPlus
still located the acquired sequence units, and the embedding representation
also includes missingness indicators and an acquired-unit count.

## Analysis chronology

| Analysis | Design timing | Reporting status |
|---|---|---|
| Main MIC endpoint and grouped model evaluation | Adopted before genomic models and external scoring | Development and locked external |
| Rule for awarding an Evo 2 gain | Adopted after human-clinical development results, before all-source Evo 2 results and external scoring | Prospectively fixed for external scoring; not prospectively fixed for all development results |
| Source robustness | Added after diagnosis of a development gain, before external scoring | Development sensitivity; evaluation rows excluded without retraining |
| Gene-family comparison on Japan | Designed after inspection of locked external errors | Exploratory/post hoc, despite locking the new predictions |
| Leave-one-variant-out comparison | Specified before its own runs, motivated by the external family analysis | Controlled development experiment; not independent external replication |

Use "prespecified before external scoring" or "specified before this
experiment" with the appropriate qualifier. The records establish dated local
protocols and original local commit metadata; they do not establish registration
in a public registry. `provenance/study_chronology.json` preserves six plan and
prediction-lock records exported before the owner-requested history replacement.

## Original artifact verification completed

The completed Unity artifacts were retrieved for release preparation. All 388
recorded artifact digests verify, with no missing hash targets, including 64
embedding shards and 67 extraction tables. The 18 committed prediction locks
match exactly. The verifier reproduces 54 evaluations, 438 paired comparisons,
192 source-exclusion summaries, the external variant subgroups, and the complete
withholding read-out. No training or external model selection was repeated.

The supplementary CSV files contain the complete original result matrix,
reference and lineage strata, selected settings, and convergence information.
The original hashes and public-copy hashes are separately retained when private
runtime paths are redacted in release copies. Prediction CSVs are unchanged.

## Specific issues resolved in this draft

- **Pooled denominator:** 568 is the number of above-panel target–isolate
  observations over four CTX-M withholding experiments. The code permits an
  isolate once per target. The verified membership table gives 565 unique isolates in 45 lineages.
- **Withholding scope:** target carriers are removed from supervised training.
  The experiment does not remove all related lineages, and pretrained
  sequence exposure remains possible.
- **Bootstrap scope:** reported intervals resample lineages conditional on
  predictions. They do not incorporate the complete model-selection process.
- **Source exclusion:** job 65016028 re-evaluates existing predictions after
  removing a source's evaluation rows. It is not leave-source-out training.
- **Quantitative endpoint:** EA against a censored reference measures distance
  to the allowed dilution set. The 2 mg/L ceftriaxone read-out is descriptive.
  Clinical S/I/R sensitivity, specificity, and very major error are not
  validated endpoints here.
- **Selection:** the external set excludes putative genomic duplicates at
  99.99% ANI and 0.95 minimum reciprocal aligned fraction. The old 454-isolate
  near-neighbor-free set is not the final 1,505-isolate evaluation population.
- **Rounding:** family ridge ceftriaxone 0.771 minus allele ridge 0.680 is
  0.091 when rounded first; the reported paired estimate is 0.092. Retain the
  original comparison rather than recomputing it from rounded values.
- **No mechanistic discovery:** recovering an unfamiliar allele or a source
  association does not establish a new biological resistance mechanism.

## Numerical and execution disclosures to retain

The original Evo 2 development round logged 14 numerical warnings in ridge
fits. A bound on residual scale removed those warnings in rerun 65006995,
with the conclusions reproduced. The interim human-clinical rerun record
reports 20 of 1,275 fits not formally converged. The external prediction
record reports 7 nonconverged fits for human-clinical AMR + Evo 2, 16 for
all-source AMR + Evo 2, and 3 for all-source Evo 2. Supplementary Table S4 and `verification/selected_settings.csv` give the complete counts, including 57/1,275 nonconverged all-source combined development fits and 15/1,275 all-source Evo 2 fits.

The successful embedding canary recorded flash-attn 2.8.3 above the declared
transformer_engine support range. Finite embeddings and checksum verification
are documented; this does not establish a separate numerical-equivalence
benchmark. Report the actual environment and checkpoint metadata.

Failed attempts are retained in the source ledger: mismatched initial MLST
scheme; empty uncertainty-column defect in external scoring 65020797; and
duplicate target–isolate handling in the first withholding read-out 65051617.
The final cited runs supersede them. The external and withholding fixes did
not change the locked predictions.

The approximately 3.6 GPU-hours describe production embedding extraction.
They exclude preceding engineering, CPU modelling, pretraining, genome
processing, and data transfer. Do not present them as the total experiment
budget or as a measured end-to-end speedup.

## Author and submission details

Rushath Rajeev confirmed independent-researcher affiliation, correspondence at
rushathrajeevav@gmail.com, UMass Unity HPC acknowledgement, no funding, and no
competing interests. The original code uses the author-selected MIT license.
Claude Code and Codex are named in the AI-assistance disclosure. Under
[arXiv's AI policy](https://info.arxiv.org/help/moderation/index.html#policy-for-authors-use-of-generative-ai-language-tools),
AI tools are not listed as authors.

The actual manuscript and supplement LaTeX sources have been compiled using
Tectonic 0.17.0. The arXiv source archive contains resolved references and vector
figures. The author must review the scientific text and the PDF produced by
arXiv before submitting. No account action, endorsement, submission, acceptance,
or peer review is asserted.

## Remaining provenance limits

Production manifests recorded `evo2_20b` and software versions, but no weight
content digest. The existing cache reference recovered during release is
`8b0f0a9a70c66367ed181a17d049b95699a28fed`; this is not represented as an
independently logged production-time hash. Published frozen embeddings have
original shard checksums, supporting exact downstream reproduction.

Production logs also record flash-attn 2.8.3 above transformer-engine's stated
support range. Checksums and finite outputs do not establish numerical
interchangeability across inference environments.

## Project ledger

The canonical status reflects the completed MIC comparison.
See `docs/SCOPE_CLOSEOUT.md` for
unperformed categorical, prospective, pangenome, abstention, and laboratory work.
