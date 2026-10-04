#!/usr/bin/env bash
# A further development modelling round on the SAME outer folds, features and
# baseline predictions as a COMPLETE round (AMR_MODELS_JOB), so every model is
# compared on identical rows. Runs each learner in AMR_LEARNERS under nested
# CV for both populations, evaluates it, and compares it with the no-genome,
# nearest-relative and round-1 ridge predictions. AMR_FEATURE_SETS adds Evo 2
# feature sets (amr+evo2, evo2; needs AMR_EMBEDDINGS and AMR_UNIT_MAPPING); each
# Evo 2 model is also compared with the same learner on amr features alone,
# which decides whether the foundation model is credited. Development only.
#SBATCH --job-name=amr-model-round
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/model-round-%j.out
#SBATCH --error=logs/model-round-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_MODELS_JOB:?Set AMR_MODELS_JOB to the COMPLETE round whose folds and features are reused}"
: "${AMR_COHORT_JOB:?Set AMR_COHORT_JOB to the COMPLETE cohort-build job}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
LEARNERS="${AMR_LEARNERS:-ridge_censored xgboost_aft}"
FEATURE_SETS="${AMR_FEATURE_SETS:-amr}"
POPULATIONS="${AMR_POPULATIONS:-human_clinical all_sources}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PY="${AMR_MODELS_ENV:-${WORK_ROOT}/conda/envs/amr-models}/bin/python"
BASE="results/models/${AMR_MODELS_JOB}"
COHORT_ROOT="results/cohort/${AMR_COHORT_JOB}"
RESULT_ROOT="results/models/${SLURM_JOB_ID}"

resolved_work="$(realpath -m "${WORK_ROOT}")"
resolved_project="$(realpath -m "${PROJECT_ROOT}")"
case "${resolved_project}" in
    "${resolved_work}"/*) ;;
    *) printf 'Project root is outside AMR_WORK_ROOT: %s\n' "${resolved_project}" >&2; exit 2 ;;
esac
cd "${resolved_project}"
[[ -x "${PY}" ]] || { printf 'Modelling environment missing: %s\n' "${PY}" >&2; exit 2; }
for root in "${BASE}" "${COHORT_ROOT}"; do
    [[ -f "${root}/COMPLETE" ]] || { printf 'Source is not COMPLETE: %s\n' "${root}" >&2; exit 2; }
    (cd "${root}" && sha256sum --check --quiet artifact.sha256)
done
[[ -e "${RESULT_ROOT}" ]] && { printf 'Refusing to overwrite: %s\n' "${RESULT_ROOT}" >&2; exit 2; }
mkdir -p "${RESULT_ROOT}"
COHORT="${COHORT_ROOT}/research_cohort.csv"
evo_args=()
if [[ "${FEATURE_SETS}" == *evo2* ]]; then
    : "${AMR_EMBEDDINGS:?Set AMR_EMBEDDINGS to the embedding shard directory}"
    : "${AMR_UNIT_MAPPING:?Set AMR_UNIT_MAPPING to unit_mapping.tsv}"
    for shard in "${AMR_EMBEDDINGS}"/shard-*.npz; do
        (cd "$(dirname "${shard}")" && sha256sum --check --quiet "$(basename "${shard%.npz}").sha256")
    done
    evo_args=(--embeddings "${AMR_EMBEDDINGS}" --unit-mapping "${AMR_UNIT_MAPPING}")
fi

for population in ${POPULATIONS}; do
    for learner in ${LEARNERS}; do
      for fs in ${FEATURE_SETS}; do
        tag="${learner}_${fs/+/_}"
        [[ "${fs}" == "amr" ]] && tag="${learner}"
        extra=()
        [[ "${fs}" != "amr" ]] && extra=("${evo_args[@]}")
        "${PY}" scripts/run_mic_models.py --learner "${learner}" --feature-set "${fs}" --cohort "${COHORT}" \
            --splits "${BASE}/isolate_splits.csv" --features "${BASE}/amr_features.tsv" "${extra[@]}" \
            --population "${population}" --output "${RESULT_ROOT}/predictions_${tag}_${population}.csv" \
            --manifest "${RESULT_ROOT}/manifest_${tag}_${population}.json"
        "${PY}" scripts/evaluate_mic_predictions.py \
            --predictions "${RESULT_ROOT}/predictions_${tag}_${population}.csv" \
            --output "${RESULT_ROOT}/evaluation_${tag}_${population}.json"
        against=(--baseline "constant=${BASE}/predictions_constant_${population}.csv"
                 --baseline "neighbor=${BASE}/predictions_neighbor_${population}.csv"
                 --baseline "ridge_round1=${BASE}/predictions_amrfinder_${population}.csv")
        if [[ "${fs}" != "amr" && -f "${RESULT_ROOT}/predictions_${learner}_${population}.csv" ]]; then
            against+=(--baseline "${learner}_amr=${RESULT_ROOT}/predictions_${learner}_${population}.csv")
        fi
        "${PY}" scripts/compare_mic_models.py \
            --model "${tag}=${RESULT_ROOT}/predictions_${tag}_${population}.csv" "${against[@]}" \
            --cohort "${COHORT}" --output "${RESULT_ROOT}/comparison_${tag}_${population}.json"
      done
    done
done

{
    printf 'job_id=%s\nbase_models_job=%s\ncohort_job=%s\nlearners=%s\nfeature_sets=%s\nembeddings=%s\n' \
        "${SLURM_JOB_ID}" "${AMR_MODELS_JOB}" "${AMR_COHORT_JOB}" "${LEARNERS}" "${FEATURE_SETS}" "${AMR_EMBEDDINGS:-none}"
    printf 'populations=%s\n' "${POPULATIONS}"
    for script in run_mic_models run_censored_regression embedding_features evaluate_mic_predictions compare_mic_models; do
        printf '%s_sha256=%s\n' "${script}" "$(sha256sum "scripts/${script}.py" | awk '{print $1}')"
    done
    printf 'models_environment=%s\n' "$(cat provenance/models-tool-versions.txt 2>/dev/null || echo unrecorded)"
    printf 'study_config_sha256=%s\n' "$(sha256sum config/study.json | awk '{print $1}')"
} > "${RESULT_ROOT}/provenance.txt"
(cd "${RESULT_ROOT}" && find . -maxdepth 1 -type f ! -name artifact.sha256 -printf '%f\n' \
    | LC_ALL=C sort | xargs sha256sum > artifact.sha256)
touch "${RESULT_ROOT}/COMPLETE"
