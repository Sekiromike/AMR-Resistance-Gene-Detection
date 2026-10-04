#!/usr/bin/env bash
# MIC_EVALUATION_SPEC addendum B: re-run the key development comparisons with
# each major source (BioProject >= 10% of a drug's rows) excluded in turn.
# Reads only COMPLETE, checksum-verified prediction sets. Development only.
#SBATCH --job-name=amr-source-robust
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/source-robustness-%j.out
#SBATCH --error=logs/source-robustness-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PY="${WORK_ROOT}/conda/envs/amr-models/bin/python"
COHORT="results/cohort/${AMR_COHORT_JOB:-64948707}/research_cohort.csv"
R1="results/models/${AMR_ROUND1_JOB:-64960056}"      # constant, neighbor, ridge round 1
R2="results/models/${AMR_ROUND2_JOB:-64965891}"      # XGBoost AFT on AMR
EVO="results/models/${AMR_EVO_JOB:-64983779}"       # XGBoost Evo 2 sets
RERUN="results/models/${AMR_RERUN_JOB:-65006995}"   # ridge sets, bounded optimiser
ABL="results/models/${AMR_ABLATION_JOB:-65006996}"  # ablation, human clinical
OUT="results/comparisons/${SLURM_JOB_ID}"

cd "${PROJECT_ROOT}"
for root in "${R1}" "${R2}" "${EVO}" "${RERUN}" "${ABL}"; do
    [[ -f "${root}/COMPLETE" ]] || { printf 'Not COMPLETE: %s\n' "${root}" >&2; exit 2; }
    (cd "${root}" && sha256sum --check --quiet artifact.sha256)
done
[[ -e "${OUT}" ]] && { printf 'Refusing to overwrite: %s\n' "${OUT}" >&2; exit 2; }
mkdir -p "${OUT}"
compare() {  # name model_csv baseline_name baseline_csv [baseline_name baseline_csv ...]
    local name="$1" model="$2"; shift 2
    local args=()
    while (( $# )); do args+=(--baseline "$1=$2"); shift 2; done
    "${PY}" scripts/compare_mic_models.py --model "${name}=${model}" "${args[@]}" \
        --cohort "${COHORT}" --output "${OUT}/${name}.json"
}
for pop in human_clinical all_sources; do
    compare "ridge_amr_${pop}" "${RERUN}/predictions_ridge_censored_${pop}.csv" \
        constant "${R1}/predictions_constant_${pop}.csv"
    compare "xgboost_amr_${pop}" "${R2}/predictions_xgboost_aft_${pop}.csv" \
        constant "${R1}/predictions_constant_${pop}.csv"
    compare "ridge_amr_evo2_${pop}" "${RERUN}/predictions_ridge_censored_amr_evo2_${pop}.csv" \
        ridge_amr "${RERUN}/predictions_ridge_censored_${pop}.csv"
    compare "xgboost_amr_evo2_${pop}" "${EVO}/predictions_xgboost_aft_amr_evo2_${pop}.csv" \
        xgboost_amr "${EVO}/predictions_xgboost_aft_${pop}.csv"
done
compare "ridge_amr_evo2_acquired_human_clinical" "${ABL}/predictions_ridge_censored_amr_evo2_acquired_human_clinical.csv" \
    ridge_amr "${ABL}/predictions_ridge_censored_human_clinical.csv"
compare "xgboost_amr_evo2_acquired_human_clinical" "${ABL}/predictions_xgboost_aft_amr_evo2_acquired_human_clinical.csv" \
    xgboost_amr "${ABL}/predictions_xgboost_aft_human_clinical.csv"
{
    printf 'job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'compare_mic_models_sha256=%s\n' "$(sha256sum scripts/compare_mic_models.py | awk '{print $1}')"
} > "${OUT}/provenance.txt"
(cd "${OUT}" && sha256sum ./*.json provenance.txt > artifact.sha256)
touch "${OUT}/COMPLETE"
