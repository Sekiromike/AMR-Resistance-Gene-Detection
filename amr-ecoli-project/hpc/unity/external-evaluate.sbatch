#!/usr/bin/env bash
# Seal step 3 (docs/EXTERNAL_EVALUATION_PLAN.md): evaluate the LOCKED external
# predictions once. attach_external_references.py refuses any file whose
# SHA-256 differs from config/external_prediction_lock.tsv.
#SBATCH --job-name=amr-external-evaluate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/external-evaluate-%j.out
#SBATCH --error=logs/external-evaluate-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_PREDICTIONS_JOB:?Set AMR_PREDICTIONS_JOB to the locked external-predict job}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PY="${WORK_ROOT}/conda/envs/amr-models/bin/python"
COHORT="results/cohort/${AMR_COHORT_JOB:-64948707}/research_cohort.csv"
PRED="results/external/predictions/${AMR_PREDICTIONS_JOB}"
LOCK="config/external_prediction_lock.tsv"
OUT="results/external/evaluation/${SLURM_JOB_ID}"

cd "${PROJECT_ROOT}"
[[ -s "${LOCK}" ]] || { printf 'No committed lock file\n' >&2; exit 2; }
[[ -e "${OUT}" ]] && { printf 'Refusing to overwrite: %s\n' "${OUT}" >&2; exit 2; }
mkdir -p "${OUT}"
while IFS=$'\t' read -r file sha; do
    [[ "${file}" == "file" ]] && continue
    "${PY}" scripts/attach_external_references.py --predictions "${PRED}/${file}" --cohort "${COHORT}" \
        --lock "${LOCK}" --output "${OUT}/scored_${file}"
    "${PY}" scripts/evaluate_mic_predictions.py --predictions "${OUT}/scored_${file}" --locked-external \
        --output "${OUT}/evaluation_${file%.csv}.json"
done < "${LOCK}"
for population in human_clinical all_sources; do
    s() { printf '%s/scored_predictions_%s_%s.csv' "${OUT}" "$1" "${population}"; }
    for m in ridge_censored_amr xgboost_aft_amr ridge_censored_amr_evo2 xgboost_aft_amr_evo2 \
             ridge_censored_evo2 xgboost_aft_evo2; do
        against=(--baseline "constant=$(s constant)" --baseline "neighbor=$(s neighbor)")
        case "${m}" in
            ridge_censored_amr_evo2|ridge_censored_evo2) against+=(--baseline "ridge_amr=$(s ridge_censored_amr)") ;;
            xgboost_aft_amr_evo2|xgboost_aft_evo2) against+=(--baseline "xgboost_amr=$(s xgboost_aft_amr)") ;;
        esac
        "${PY}" scripts/compare_mic_models.py --model "${m}=$(s "${m}")" "${against[@]}" \
            --cohort "${COHORT}" --output "${OUT}/comparison_${m}_${population}.json"
    done
done
{
    printf 'job_id=%s\npredictions_job=%s\nlock_sha256=%s\n' "${SLURM_JOB_ID}" "${AMR_PREDICTIONS_JOB}" \
        "$(sha256sum "${LOCK}" | awk '{print $1}')"
    for script in attach_external_references evaluate_mic_predictions compare_mic_models; do
        printf '%s_sha256=%s\n' "${script}" "$(sha256sum "scripts/${script}.py" | awk '{print $1}')"
    done
} > "${OUT}/provenance.txt"
(cd "${OUT}" && find . -maxdepth 1 -type f ! -name artifact.sha256 -printf '%f\n' | LC_ALL=C sort | xargs sha256sum > artifact.sha256)
touch "${OUT}/COMPLETE"
