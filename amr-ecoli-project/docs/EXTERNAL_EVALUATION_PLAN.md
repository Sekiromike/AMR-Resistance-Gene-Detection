# Locked External Evaluation Plan

Status: ADOPTED 2026-09-29 by the decision owner
Decision owner: Rushath Rajeev
Governs: the single evaluation of the external-locked set (amendments 002, 003,
007; `MIC_EVALUATION_SPEC.md` with addenda A and B;
`FOUNDATION_MODEL_PROTOCOL.md` revision 0.2 and its credit rule)

The external set is evaluated **once**. Everything below is fixed before any
external MIC value, comparator or derived statistic is read.

## Evaluation set

External-locked JARBS isolates in validated cohort `64948707`: 1,505 isolates,
each with ceftriaxone, ciprofloxacin and gentamicin MICs (4,515 rows). All are
human clinical by study design. They are genomically disjoint from development
at the putative-duplicate rule (amendment 003), from one laboratory and one
method (central broth microdilution, MicroScan WalkAway panels), in Japan.

## Frozen models

Every model is fitted on development only, with hyperparameters chosen by the
same grouped inner cross-validation used during development (five folds over
the amendment 007 CV units of the whole development set), then refitted on all
of development. No choice uses the external set.

| # | Model | Role |
|---|---|---|
| 1 | No-genome constant (best constant dilution per drug) | baseline |
| 2 | Nearest development relative by ANI (frozen external-vs-development comparisons `63772662`) | baseline |
| 3 | ridge_censored on AMRFinderPlus | known mechanisms |
| 4 | xgboost_aft on AMRFinderPlus | known mechanisms |
| 5 | ridge_censored on AMRFinderPlus + Evo 2 | Evo 2 credit test |
| 6 | xgboost_aft on AMRFinderPlus + Evo 2 | Evo 2 credit test |
| 7 | ridge_censored on Evo 2 alone | descriptive |
| 8 | xgboost_aft on Evo 2 alone | descriptive |

Primary training population: human clinical. Sensitivity: the same eight
models trained on all development sources, evaluated on the same rows.

## Seal protocol

1. **Predict without reading labels.** The prediction step loads external rows
   by isolate and antibiotic with genotype features only; it never reads the
   MIC value or comparator columns.
2. **Lock.** Prediction files are checksummed, and the hashes are recorded in the
   project ledger and committed **before** any evaluation runs.
3. **Evaluate once.** `evaluate_mic_predictions.py --locked-external` and
   `compare_mic_models.py` run once on the locked files. A rerun is allowed only
   for a verified software defect; then both results are reported.

## Reported, for every model and drug

- Essential agreement with the whole-lineage bootstrap over external lineage
  groups; exact-reference agreement.
- Reference classes (addendum A), with the right-censored under-call rate
  beside the ~1.5% very-major-error threshold.
- SLV:ST131 and other lineages separately; isolates in lineage groups absent
  from development (18) listed descriptively only.
- Paired differences: every genomic model versus both baselines; models 5 and
  6 versus 3 and 4.
- Precision labels fixed by amendment 007: ceftriaxone in the primary
  population is labelled underpowered (planned half-width 0.052 > 0.05).

## Decisions this evaluation makes

- **Evo 2 credit**, per the pre-registered rule: condition 2 is tested for
  gentamicin, the only drug that met condition 1; the same comparison is
  reported for ceftriaxone and ciprofloxacin. The development diagnosis that
  the gentamicin gain came from one source (PRJNA1297298) is reported beside it.
- **Transfer of known-mechanism models** from a development population
  dominated by one US hospital to one Japanese laboratory.

## What the external set cannot answer

It is not representative. Its isolates were selected as the JARBS isolates
most distinct from development, and JARBS enriched for resistance. Results are
a test of transfer under strong shift, not an estimate of deployed accuracy.
No clinical claim follows, whatever the result.
