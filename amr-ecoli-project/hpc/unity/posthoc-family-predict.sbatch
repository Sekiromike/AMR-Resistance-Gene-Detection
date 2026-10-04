#!/usr/bin/env bash
# docs/POSTHOC_FAMILY_ANALYSIS.md steps 1-2: family-level catalogue features,
# development nested CV (reference) and external predictions without reading
# external MICs. Scoring waits for config/external_posthoc_lock.tsv.
#SBATCH --job-name=amr-posthoc-family
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/posthoc-family-%j.out
#SBATCH --error=logs/posthoc-family-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PY="${WORK_ROOT}/conda/envs/amr-models/bin/python"
COHORT="results/cohort/${AMR_COHORT_JOB:-64948707}/research_cohort.csv"
BASE="results/models/${AMR_MODELS_JOB:-64960056}"
OUT="results/posthoc/family/${SLURM_JOB_ID}"

cd "${PROJECT_ROOT}"
(cd "${BASE}" && sha256sum --check --quiet artifact.sha256)
mkdir -p "${OUT}"
"${PY}" scripts/build_hierarchical_features.py --features "${BASE}/amr_features.tsv" --output "${OUT}/amr_family_features.tsv"
for learner in ridge_censored xgboost_aft; do
    "${PY}" scripts/run_mic_models.py --learner "${learner}" --cohort "${COHORT}" --splits "${BASE}/isolate_splits.csv" \
        --features "${OUT}/amr_family_features.tsv" --population human_clinical \
        --output "${OUT}/dev_predictions_${learner}_family_human_clinical.csv" \
        --manifest "${OUT}/dev_manifest_${learner}_family_human_clinical.json"
    "${PY}" scripts/evaluate_mic_predictions.py --predictions "${OUT}/dev_predictions_${learner}_family_human_clinical.csv" \
        --output "${OUT}/dev_evaluation_${learner}_family_human_clinical.json"
    "${PY}" scripts/predict_external.py --model "${learner}:amr" --population human_clinical --cohort "${COHORT}" \
        --splits "${BASE}/isolate_splits.csv" --features "${OUT}/amr_family_features.tsv" \
        --output "${OUT}/predictions_${learner}_family_human_clinical.csv" \
        --manifest "${OUT}/manifest_${learner}_family_human_clinical.json"
done
(cd "${OUT}" && sha256sum ./*.tsv ./*.csv ./*.json > artifact.sha256)
touch "${OUT}/COMPLETE"
