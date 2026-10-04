#!/usr/bin/env bash
# Seal step 1 (docs/EXTERNAL_EVALUATION_PLAN.md): frozen-model predictions for
# the external-locked set. One array task per model; each writes predictions
# for both training populations. External MIC columns are never read
# (predict_external.py drops them while streaming the cohort).
#SBATCH --job-name=amr-external-predict
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-7
#SBATCH --output=logs/external-predict-%A_%a.out
#SBATCH --error=logs/external-predict-%A_%a.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${SLURM_ARRAY_JOB_ID:?Submit this script as a Slurm array}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PY="${WORK_ROOT}/conda/envs/amr-models/bin/python"
COHORT_JOB="${AMR_COHORT_JOB:-64948707}"
BASE="results/models/${AMR_MODELS_JOB:-64960056}"
COHORT="results/cohort/${COHORT_JOB}/research_cohort.csv"
EMB="results/genomics/embeddings/loci/${COHORT_JOB}/evo2_20b"
MAPPING="results/genomics/loci/${COHORT_JOB}/unit_mapping.tsv"
COMPARISONS="results/genomics/disjointness/comparisons/${AMR_EXTERNAL_COMPARISONS_JOB:-63772662}/comparisons.tsv"
OUT="results/external/predictions/${SLURM_ARRAY_JOB_ID}"
MODELS=(constant neighbor ridge_censored:amr xgboost_aft:amr ridge_censored:amr+evo2 xgboost_aft:amr+evo2
        ridge_censored:evo2 xgboost_aft:evo2)
MODEL="${MODELS[${SLURM_ARRAY_TASK_ID}]}"

cd "${PROJECT_ROOT}"
for root in "results/cohort/${COHORT_JOB}" "${BASE}" "results/genomics/loci/${COHORT_JOB}"; do
    [[ -f "${root}/COMPLETE" ]] || { printf 'Not COMPLETE: %s\n' "${root}" >&2; exit 2; }
    (cd "${root}" && sha256sum --check --quiet artifact.sha256)
done
expected="$(awk '$2 ~ /comparisons.tsv$/ {print $1}' "$(dirname "${COMPARISONS}")/artifact.sha256")"
[[ "$(sha256sum "${COMPARISONS}" | awk '{print $1}')" == "${expected}" ]] || { printf 'comparisons changed\n' >&2; exit 2; }
mkdir -p "${OUT}"
tag="${MODEL/:/_}"; tag="${tag/+/_}"
for population in human_clinical all_sources; do
    "${PY}" scripts/predict_external.py --model "${MODEL}" --population "${population}" \
        --cohort "${COHORT}" --splits "${BASE}/isolate_splits.csv" --features "${BASE}/amr_features.tsv" \
        --embeddings "${EMB}" --unit-mapping "${MAPPING}" --comparisons "${COMPARISONS}" \
        --output "${OUT}/predictions_${tag}_${population}.csv" \
        --manifest "${OUT}/manifest_${tag}_${population}.json"
done
printf '%s\n' "$(sha256sum scripts/predict_external.py | awk '{print $1}')" > "${OUT}/predict_external.sha256.task${SLURM_ARRAY_TASK_ID}"
