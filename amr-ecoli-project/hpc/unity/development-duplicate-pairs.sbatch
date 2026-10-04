#!/usr/bin/env bash
# Extract putative genomic duplicate pairs (amendment 005, rule G) from the
# completed development-versus-development comparison table. Streams the table,
# so memory stays small even at ~106M rows. Submit with
# --dependency=afterok:<development-clusters resume job>.
#SBATCH --job-name=amr-dev-duplicate-pairs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/development-duplicate-pairs-%j.out
#SBATCH --error=logs/development-duplicate-pairs-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_CLUSTERS_JOB:?Set AMR_CLUSTERS_JOB to the COMPLETE development-clusters job}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
CONTROLLER_ENV="${AMR_CONTROLLER_ENV:-${WORK_ROOT}/conda/envs/amr-ecoli-controller}"
CONTRACT="data/genomics/disjointness/requests/disjointness_contract.json"
SOURCE_ROOT="results/genomics/disjointness/development-clusters/${AMR_CLUSTERS_JOB}"
RESULT_ROOT="results/genomics/disjointness/duplicate-pairs/${SLURM_JOB_ID}"

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
if [[ ! -f "${SOURCE_ROOT}/COMPLETE" ]]; then
    printf 'Development clustering is not COMPLETE: %s\n' "${SOURCE_ROOT}" >&2
    exit 2
fi
(cd "${SOURCE_ROOT}" && sha256sum --check --quiet artifact.sha256)
if [[ -e "${RESULT_ROOT}" ]]; then
    printf 'Refusing to overwrite duplicate-pair result: %s\n' "${RESULT_ROOT}" >&2
    exit 2
fi
mkdir -p "${RESULT_ROOT}"

"${CONTROLLER_ENV}/bin/python" scripts/build_genomic_disjointness.py extract-duplicate-pairs \
    --contract "${CONTRACT}" \
    --comparisons "${SOURCE_ROOT}/comparisons.tsv" \
    --pairs "${RESULT_ROOT}/development-duplicate-pairs.tsv" \
    --manifest "${RESULT_ROOT}/duplicate-pairs-manifest.json"

{
    printf 'source_clusters_job=%s\n' "${AMR_CLUSTERS_JOB}"
    printf 'builder_sha256=%s\n' "$(sha256sum scripts/build_genomic_disjointness.py | awk '{print $1}')"
} > "${RESULT_ROOT}/provenance.txt"
(cd "${RESULT_ROOT}" && sha256sum provenance.txt development-duplicate-pairs.tsv \
    duplicate-pairs-manifest.json > artifact.sha256)
touch "${RESULT_ROOT}/COMPLETE"
