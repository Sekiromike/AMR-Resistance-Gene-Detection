#!/usr/bin/env bash
# docs/POSTHOC_FAMILY_ANALYSIS.md step 3: score the LOCKED family-model
# external predictions (config/external_posthoc_lock.tsv) and read out the
# unseen-variant shares beside the locked allele-level and Evo 2 models.
#SBATCH --job-name=amr-posthoc-score
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/posthoc-score-%j.out
#SBATCH --error=logs/posthoc-score-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_POSTHOC_JOB:?Set AMR_POSTHOC_JOB to the posthoc-family-predict job}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PY="${WORK_ROOT}/conda/envs/amr-models/bin/python"
COHORT="results/cohort/${AMR_COHORT_JOB:-64948707}/research_cohort.csv"
PH="results/posthoc/family/${AMR_POSTHOC_JOB}"
MAIN="results/external/evaluation/${AMR_MAIN_SCORING_JOB:-65020832}"
LOCK="config/external_posthoc_lock.tsv"
OUT="${PH}/scoring-${SLURM_JOB_ID}"

cd "${PROJECT_ROOT}"
[[ -s "${LOCK}" ]] || { printf 'No committed post-hoc lock\n' >&2; exit 2; }
(cd "${MAIN}" && sha256sum --check --quiet artifact.sha256)
mkdir -p "${OUT}"
for learner in ridge_censored xgboost_aft; do
    f="predictions_${learner}_family_human_clinical.csv"
    "${PY}" scripts/attach_external_references.py --predictions "${PH}/${f}" --cohort "${COHORT}" --lock "${LOCK}" \
        --output "${OUT}/scored_${f}"
    "${PY}" scripts/evaluate_mic_predictions.py --predictions "${OUT}/scored_${f}" --locked-external \
        --output "${OUT}/evaluation_${learner}_family_human_clinical.json"
    "${PY}" scripts/compare_mic_models.py --model "${learner}_family=${OUT}/scored_${f}" \
        --baseline "${learner}_allele=${MAIN}/scored_predictions_${learner}_amr_human_clinical.csv" \
        --baseline "constant=${MAIN}/scored_predictions_constant_human_clinical.csv" \
        --cohort "${COHORT}" --output "${OUT}/comparison_${learner}_family_human_clinical.json"
done
"${PY}" scripts/unseen_variant_analysis.py --features "${PH}/amr_family_features.tsv" --output "${OUT}/unseen_variants.json" \
    --scored "ridge_allele=${MAIN}/scored_predictions_ridge_censored_amr_human_clinical.csv" \
             "xgboost_allele=${MAIN}/scored_predictions_xgboost_aft_amr_human_clinical.csv" \
             "ridge_family=${OUT}/scored_predictions_ridge_censored_family_human_clinical.csv" \
             "xgboost_family=${OUT}/scored_predictions_xgboost_aft_family_human_clinical.csv" \
             "ridge_evo2=${MAIN}/scored_predictions_ridge_censored_evo2_human_clinical.csv" \
             "xgboost_evo2=${MAIN}/scored_predictions_xgboost_aft_evo2_human_clinical.csv" \
             "ridge_amr_evo2=${MAIN}/scored_predictions_ridge_censored_amr_evo2_human_clinical.csv" \
             "xgboost_amr_evo2=${MAIN}/scored_predictions_xgboost_aft_amr_evo2_human_clinical.csv"
(cd "${OUT}" && sha256sum ./*.csv ./*.json > artifact.sha256)
touch "${OUT}/COMPLETE"
