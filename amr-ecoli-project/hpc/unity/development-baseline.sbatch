#!/usr/bin/env bash
# Development folds and the no-genome MIC baseline on a COMPLETE cohort:
#   1. lineage-grouped development folds (build_split_manifest.py);
#   2. per population, the intercept-only interval-censored baseline
#      (run_mic_baseline.py) and its evaluation (evaluate_mic_predictions.py).
# Development rows only: the external set is not read by steps 2-3 and the
# split manifest reports external rows as counts per antibiotic.
#SBATCH --job-name=amr-dev-baseline
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/development-baseline-%j.out
#SBATCH --error=logs/development-baseline-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_COHORT_JOB:?Set AMR_COHORT_JOB to the COMPLETE cohort-build job}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PY="${AMR_CONTROLLER_ENV:-${WORK_ROOT}/conda/envs/amr-ecoli-controller}/bin/python"
COHORT_ROOT="results/cohort/${AMR_COHORT_JOB}"
RESULT_ROOT="results/baselines/${SLURM_JOB_ID}"

resolved_work="$(realpath -m "${WORK_ROOT}")"
resolved_project="$(realpath -m "${PROJECT_ROOT}")"
case "${resolved_project}" in
    "${resolved_work}"/*) ;;
    *)
        printf 'Project root is outside AMR_WORK_ROOT: %s\n' "${resolved_project}" >&2
        exit 2
        ;;
esac
cd "${resolved_project}"
[[ -f "${COHORT_ROOT}/COMPLETE" ]] || { printf 'Cohort is not COMPLETE: %s\n' "${COHORT_ROOT}" >&2; exit 2; }
(cd "${COHORT_ROOT}" && sha256sum --check --quiet artifact.sha256)
if [[ -e "${RESULT_ROOT}" ]]; then
    printf 'Refusing to overwrite baseline result: %s\n' "${RESULT_ROOT}" >&2
    exit 2
fi
mkdir -p "${RESULT_ROOT}"
read -r FOLDS SEED < <("${PY}" -c '
import json
s = json.load(open("config/study.json"))["splitting"]
print(s["folds"], s["seed"])')

"${PY}" scripts/build_split_manifest.py \
    --cohort "${COHORT_ROOT}/research_cohort.csv" \
    --output "${RESULT_ROOT}/isolate_splits.csv" \
    --manifest "${RESULT_ROOT}/split_manifest.json" \
    --folds "${FOLDS}" --seed "${SEED}"

for population in human_clinical all_sources; do
    "${PY}" scripts/run_mic_baseline.py \
        --cohort "${COHORT_ROOT}/research_cohort.csv" \
        --splits "${RESULT_ROOT}/isolate_splits.csv" \
        --population "${population}" \
        --output "${RESULT_ROOT}/baseline_predictions_${population}.csv" \
        --manifest "${RESULT_ROOT}/baseline_manifest_${population}.json"
    "${PY}" scripts/evaluate_mic_predictions.py \
        --predictions "${RESULT_ROOT}/baseline_predictions_${population}.csv" \
        --output "${RESULT_ROOT}/baseline_evaluation_${population}.json"
done

{
    printf 'job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'cohort_job=%s\n' "${AMR_COHORT_JOB}"
    for script in build_split_manifest run_mic_baseline evaluate_mic_predictions; do
        printf '%s_sha256=%s\n' "${script}" "$(sha256sum "scripts/${script}.py" | awk '{print $1}')"
    done
    printf 'study_config_sha256=%s\n' "$(sha256sum config/study.json | awk '{print $1}')"
} > "${RESULT_ROOT}/provenance.txt"
(cd "${RESULT_ROOT}" && find . -maxdepth 1 -type f ! -name artifact.sha256 -printf '%f\n' \
    | LC_ALL=C sort | xargs sha256sum > artifact.sha256)
touch "${RESULT_ROOT}/COMPLETE"
