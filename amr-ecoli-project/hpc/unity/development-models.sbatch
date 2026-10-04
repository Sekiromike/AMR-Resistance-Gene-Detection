#!/usr/bin/env bash
# Development modelling round on a COMPLETE cohort (amendment 007):
#   1. lineage/cluster/duplicate-grouped outer folds (build_split_manifest.py);
#   2. no-genome best-constant baseline (run_mic_baseline.py);
#   3. nearest-relative baseline over the frozen dev-vs-dev comparisons;
#   4. AMRFinderPlus AMR determinants (build_amr_features.py);
#   5. known-mechanism censored regression under nested CV;
#   6. evaluation of every model on identical rows, per population.
# Development only: no external MIC value or comparator is read.
#SBATCH --job-name=amr-dev-models
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/development-models-%j.out
#SBATCH --error=logs/development-models-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_COHORT_JOB:?Set AMR_COHORT_JOB to the COMPLETE cohort-build job}"
: "${AMR_DEV_CLUSTERS_JOB:?Set AMR_DEV_CLUSTERS_JOB to the COMPLETE development-clusters job}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PY="${AMR_CONTROLLER_ENV:-${WORK_ROOT}/conda/envs/amr-ecoli-controller}/bin/python"
COHORT_ROOT="results/cohort/${AMR_COHORT_JOB}"
CLUSTERS_ROOT="results/genomics/disjointness/development-clusters/${AMR_DEV_CLUSTERS_JOB}"
RESULT_ROOT="results/models/${SLURM_JOB_ID}"

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
for root in "${COHORT_ROOT}" "${CLUSTERS_ROOT}"; do
    [[ -f "${root}/COMPLETE" ]] || { printf 'Source is not COMPLETE: %s\n' "${root}" >&2; exit 2; }
done
(cd "${COHORT_ROOT}" && sha256sum --check --quiet artifact.sha256)
expected="$(awk '$2 ~ /comparisons.tsv$/ {print $1}' "${CLUSTERS_ROOT}/artifact.sha256")"
[[ -n "${expected}" ]] || { printf 'comparisons.tsv not recorded in %s\n' "${CLUSTERS_ROOT}" >&2; exit 2; }
[[ "$(sha256sum "${CLUSTERS_ROOT}/comparisons.tsv" | awk '{print $1}')" == "${expected}" ]] \
    || { printf 'comparisons.tsv changed since clustering\n' >&2; exit 2; }
if [[ -e "${RESULT_ROOT}" ]]; then
    printf 'Refusing to overwrite model result: %s\n' "${RESULT_ROOT}" >&2
    exit 2
fi
mkdir -p "${RESULT_ROOT}"
COHORT="${COHORT_ROOT}/research_cohort.csv"
read -r FOLDS SEED < <("${PY}" -c '
import json
s = json.load(open("config/study.json"))["splitting"]
print(s["folds"], s["seed"])')

"${PY}" scripts/build_split_manifest.py --cohort "${COHORT}" \
    --output "${RESULT_ROOT}/isolate_splits.csv" --manifest "${RESULT_ROOT}/split_manifest.json" \
    --folds "${FOLDS}" --seed "${SEED}"
"${PY}" scripts/build_amr_features.py --genome-manifest "${COHORT_ROOT}/genome_manifest.tsv" \
    --output "${RESULT_ROOT}/amr_features.tsv" --manifest "${RESULT_ROOT}/amr_features.json"
"${PY}" scripts/run_neighbor_baseline.py --cohort "${COHORT}" --splits "${RESULT_ROOT}/isolate_splits.csv" \
    --comparisons "${CLUSTERS_ROOT}/comparisons.tsv" \
    --output-prefix "${RESULT_ROOT}/predictions_neighbor" --manifest "${RESULT_ROOT}/neighbor_manifest.json"

for population in human_clinical all_sources; do
    "${PY}" scripts/run_mic_baseline.py --cohort "${COHORT}" --splits "${RESULT_ROOT}/isolate_splits.csv" \
        --population "${population}" --output "${RESULT_ROOT}/predictions_constant_${population}.csv" \
        --manifest "${RESULT_ROOT}/constant_manifest_${population}.json"
    "${PY}" scripts/run_censored_regression.py --cohort "${COHORT}" --splits "${RESULT_ROOT}/isolate_splits.csv" \
        --features "${RESULT_ROOT}/amr_features.tsv" --population "${population}" \
        --output "${RESULT_ROOT}/predictions_amrfinder_${population}.csv" \
        --manifest "${RESULT_ROOT}/amrfinder_manifest_${population}.json"
    for model in constant neighbor amrfinder; do
        "${PY}" scripts/evaluate_mic_predictions.py \
            --predictions "${RESULT_ROOT}/predictions_${model}_${population}.csv" \
            --output "${RESULT_ROOT}/evaluation_${model}_${population}.json"
    done
done

{
    printf 'job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'cohort_job=%s\n' "${AMR_COHORT_JOB}"
    printf 'development_clusters_job=%s\n' "${AMR_DEV_CLUSTERS_JOB}"
    for script in build_split_manifest build_amr_features run_neighbor_baseline run_mic_baseline \
                  run_censored_regression evaluate_mic_predictions; do
        printf '%s_sha256=%s\n' "${script}" "$(sha256sum "scripts/${script}.py" | awk '{print $1}')"
    done
    printf 'study_config_sha256=%s\n' "$(sha256sum config/study.json | awk '{print $1}')"
} > "${RESULT_ROOT}/provenance.txt"
(cd "${RESULT_ROOT}" && find . -maxdepth 1 -type f ! -name artifact.sha256 -printf '%f\n' \
    | LC_ALL=C sort | xargs sha256sum > artifact.sha256)
touch "${RESULT_ROOT}/COMPLETE"
