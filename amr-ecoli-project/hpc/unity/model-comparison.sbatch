#!/usr/bin/env bash
# Paired model-versus-baseline comparison and external precision planning
# (amendment 007) for a COMPLETE development modelling round. Reads external
# lineage-group sizes only; no external MIC value or comparator.
#SBATCH --job-name=amr-model-compare
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/model-comparison-%j.out
#SBATCH --error=logs/model-comparison-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_MODELS_JOB:?Set AMR_MODELS_JOB to the COMPLETE development-models job}"
: "${AMR_COHORT_JOB:?Set AMR_COHORT_JOB to the COMPLETE cohort-build job}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PY="${AMR_CONTROLLER_ENV:-${WORK_ROOT}/conda/envs/amr-ecoli-controller}/bin/python"
MODELS_ROOT="results/models/${AMR_MODELS_JOB}"
COHORT_ROOT="results/cohort/${AMR_COHORT_JOB}"
RESULT_ROOT="results/comparisons/${SLURM_JOB_ID}"

resolved_work="$(realpath -m "${WORK_ROOT}")"
resolved_project="$(realpath -m "${PROJECT_ROOT}")"
case "${resolved_project}" in
    "${resolved_work}"/*) ;;
    *) printf 'Project root is outside AMR_WORK_ROOT: %s\n' "${resolved_project}" >&2; exit 2 ;;
esac
cd "${resolved_project}"
for root in "${MODELS_ROOT}" "${COHORT_ROOT}"; do
    [[ -f "${root}/COMPLETE" ]] || { printf 'Source is not COMPLETE: %s\n' "${root}" >&2; exit 2; }
    (cd "${root}" && sha256sum --check --quiet artifact.sha256)
done
[[ -e "${RESULT_ROOT}" ]] && { printf 'Refusing to overwrite: %s\n' "${RESULT_ROOT}" >&2; exit 2; }
mkdir -p "${RESULT_ROOT}"

for population in human_clinical all_sources; do
    "${PY}" scripts/compare_mic_models.py \
        --model "amrfinder=${MODELS_ROOT}/predictions_amrfinder_${population}.csv" \
        --baseline "constant=${MODELS_ROOT}/predictions_constant_${population}.csv" \
        --baseline "neighbor=${MODELS_ROOT}/predictions_neighbor_${population}.csv" \
        --cohort "${COHORT_ROOT}/research_cohort.csv" \
        --output "${RESULT_ROOT}/comparison_amrfinder_${population}.json"
done
{
    printf 'job_id=%s\nmodels_job=%s\ncohort_job=%s\n' "${SLURM_JOB_ID}" "${AMR_MODELS_JOB}" "${AMR_COHORT_JOB}"
    printf 'compare_mic_models_sha256=%s\n' "$(sha256sum scripts/compare_mic_models.py | awk '{print $1}')"
} > "${RESULT_ROOT}/provenance.txt"
(cd "${RESULT_ROOT}" && find . -maxdepth 1 -type f ! -name artifact.sha256 -printf '%f\n' \
    | LC_ALL=C sort | xargs sha256sum > artifact.sha256)
touch "${RESULT_ROOT}/COMPLETE"
