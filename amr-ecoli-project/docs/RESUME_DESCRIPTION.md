# Resume Description

## Project

**Predicting Antibiotic Resistance in E. coli from Genomic Data**

Developed a reproducible bioinformatics and machine learning pipeline to
predict minimum inhibitory concentrations (MICs) for ceftriaxone,
ciprofloxacin, and gentamicin from E. coli genomic features. Compared
resistance-gene annotations, Evo 2 DNA embeddings, and combined features,
and evaluated how predictions transfer across bacterial lineages and
geographic populations.

**Technologies:** Python, pandas, NumPy, SciPy, scikit-learn, XGBoost,
PyTorch, Evo 2, AMRFinderPlus, Snakemake, Slurm, Linux, Git.

## Resume bullets

- Built a reproducible genomics and machine learning workflow, integrating
  29,619 antibiotic susceptibility measurements across 11,520 E. coli
  isolates into a curated research cohort.
- Implemented regression models that account for MIC measurement bounds
  using regularized regression and XGBoost; benchmarked resistance-gene
  features against Evo 2 embeddings with nested cross-validation grouped
  by bacterial lineage and bootstrap confidence intervals.
- Evaluated frozen models on 1,505 independent Japanese isolates, achieving
  93.3% MIC essential agreement for ciprofloxacin with XGBoost and 88.7%
  for gentamicin with regularized regression.
- Investigated errors caused by source differences and unseen resistance
  variants; engineered gene-family features that increased ceftriaxone
  MIC essential agreement from 59.7% to 85.4% in an exploratory follow-up
  analysis on the same external cohort.

## Project outputs

- Curated isolate and quantitative susceptibility tables, exclusion ledgers,
  genome annotations, and reproducible split manifests.
- Resistance-gene feature matrices and Evo 2 sequence embeddings.
- Per-isolate, per-antibiotic MIC predictions, prediction hashes, and model
  and evaluation manifests.
- Reports comparing models on development folds and the locked external
  cohort, with lineage-bootstrap uncertainty and error analyses.
- Snakemake workflows, Slurm batch scripts, and automated software tests.

## Metric and evidence notes

Essential agreement means a prediction is within one doubling dilution of
the reference MIC or its allowed interval when the measurement is reported
as a bound. These percentages measure MIC agreement. Clinical use has not
been validated.

The gene-family result is post hoc and requires confirmation on a new
independent cohort. Evo 2 did not demonstrate the prespecified improvement
over AMRFinderPlus features in the locked external evaluation.

This description summarizes existing repository evidence; it does not
represent a new model run or independent reproduction of the HPC results.

### Evidence sources

- [Project ledger](PROJECT_STATUS.json): cohort build 64948707 records
  29,619 rows, 11,520 isolates, 10,015 development isolates, and 1,505
  external isolates; genome processing and feature extraction are recorded
  separately.
- [External results](EXTERNAL_RESULTS.md): locked external cohort and
  per-drug essential agreement.
- [Post-hoc family analysis](POSTHOC_FAMILY_ANALYSIS.md): exploratory
  ceftriaxone comparison, 0.597 to 0.854.
- [MIC evaluation specification](MIC_EVALUATION_SPEC.md): essential
  agreement definition and whole-lineage bootstrap.
- [Model implementation](../scripts/run_mic_models.py),
  [embedding implementation](../scripts/embed_resistance_units.py), and
  [workflow](../workflow/Snakefile): implemented models and tools.
