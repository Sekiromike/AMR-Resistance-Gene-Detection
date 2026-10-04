#!/usr/bin/env bash
# Submit only with an afterok dependency on a complete assembly-audit job.
#SBATCH --job-name=amr-disjoint-compare
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/disjointness-compare-%j.out
#SBATCH --error=logs/disjointness-compare-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
CONTROLLER_ENV="${AMR_CONTROLLER_ENV:-${WORK_ROOT}/conda/envs/amr-ecoli-controller}"
DISJOINTNESS_ENV="${AMR_DISJOINTNESS_ENV:-${WORK_ROOT}/conda/envs/amr-genomic-disjointness}"
REQUESTS="data/genomics/disjointness/requests/sequence_requests.tsv"
CONTRACT="data/genomics/disjointness/requests/disjointness_contract.json"
AUDIT="results/genomics/disjointness/assembly-audit-manifest.json"
RESULT_ROOT="results/genomics/disjointness/comparisons/${SLURM_JOB_ID}"

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
if [[ ! -x "${CONTROLLER_ENV}/bin/python" || ! -x "${DISJOINTNESS_ENV}/bin/skani" ]]; then
    printf 'Verified controller or disjointness environment is unavailable\n' >&2
    exit 2
fi
if [[ -e "${RESULT_ROOT}" ]]; then
    printf 'Refusing to overwrite comparison result: %s\n' "${RESULT_ROOT}" >&2
    exit 2
fi
mkdir -p "${RESULT_ROOT}/tmp"
export TMPDIR="${resolved_project}/${RESULT_ROOT}/tmp"

"${CONTROLLER_ENV}/bin/python" scripts/build_genomic_disjointness.py prepare-comparison \
    --requests "${REQUESTS}" \
    --audit-manifest "${AUDIT}" \
    --assemblies-root results/genomics/disjointness/assemblies \
    --development-list "${RESULT_ROOT}/development-assemblies.txt" \
    --external-list "${RESULT_ROOT}/external-assemblies.txt" \
    --path-map "${RESULT_ROOT}/assembly-path-map.tsv" \
    --completed-external "${RESULT_ROOT}/completed-external.txt" \
    --manifest "${RESULT_ROOT}/comparison-input-manifest.json" \
    --workers "${SLURM_CPUS_PER_TASK}"

"${DISJOINTNESS_ENV}/bin/skani" sketch \
    -l "${RESULT_ROOT}/development-assemblies.txt" \
    -o "${RESULT_ROOT}/development-sketch" \
    -t "${SLURM_CPUS_PER_TASK}" -c 125 -m 1000
"${DISJOINTNESS_ENV}/bin/skani" search \
    -d "${RESULT_ROOT}/development-sketch" \
    --ql "${RESULT_ROOT}/external-assemblies.txt" \
    -o "${RESULT_ROOT}/skani-search.raw.tsv" \
    -t "${SLURM_CPUS_PER_TASK}" -s 95 --both-min-af 50 --keep-refs

"${CONTROLLER_ENV}/bin/python" scripts/build_genomic_disjointness.py normalize-skani \
    --raw "${RESULT_ROOT}/skani-search.raw.tsv" \
    --path-map "${RESULT_ROOT}/assembly-path-map.tsv" \
    --output "${RESULT_ROOT}/comparisons.tsv" \
    --manifest "${RESULT_ROOT}/normalization-manifest.json"
"${CONTROLLER_ENV}/bin/python" scripts/build_genomic_disjointness.py evaluate \
    --requests "${REQUESTS}" \
    --contract "${CONTRACT}" \
    --comparisons "${RESULT_ROOT}/comparisons.tsv" \
    --completed-external "${RESULT_ROOT}/completed-external.txt" \
    --collisions "${RESULT_ROOT}/collisions.tsv" \
    --manifest "${RESULT_ROOT}/evaluation-manifest.json"

{
    printf 'completed_at_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'skani_version=%s\n' "$("${DISJOINTNESS_ENV}/bin/skani" --version)"
    printf 'skani_sketch_parameters=-c 125 -m 1000\n'
    printf 'skani_search_parameters=-s 95 --both-min-af 50 --keep-refs\n'
    sha256sum "${REQUESTS}" "${CONTRACT}" "${AUDIT}"
} > "${RESULT_ROOT}/provenance.txt"
(
    cd "${RESULT_ROOT}"
    sha256sum development-assemblies.txt external-assemblies.txt \
        assembly-path-map.tsv completed-external.txt comparison-input-manifest.json \
        skani-search.raw.tsv comparisons.tsv normalization-manifest.json \
        collisions.tsv evaluation-manifest.json provenance.txt > artifact.sha256
)
touch "${RESULT_ROOT}/COMPLETE"
