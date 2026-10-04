#!/usr/bin/env bash
# External-versus-external genomic clustering.
# Reuses the frozen sketch and path map from a COMPLETE comparison run; supply
# it as AMR_SOURCE_COMPARISON (a job id under results/genomics/disjointness/comparisons).
#SBATCH --job-name=amr-ext-clusters
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/external-clusters-%j.out
#SBATCH --error=logs/external-clusters-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_SOURCE_COMPARISON:?Set AMR_SOURCE_COMPARISON to the COMPLETE comparison job id}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
CONTROLLER_ENV="${AMR_CONTROLLER_ENV:-${WORK_ROOT}/conda/envs/amr-ecoli-controller}"
DISJOINTNESS_ENV="${AMR_DISJOINTNESS_ENV:-${WORK_ROOT}/conda/envs/amr-genomic-disjointness}"
REQUESTS="data/genomics/disjointness/requests/sequence_requests.tsv"
CONTRACT="data/genomics/disjointness/requests/disjointness_contract.json"
SOURCE_ROOT="results/genomics/disjointness/comparisons/${AMR_SOURCE_COMPARISON}"
RESULT_ROOT="results/genomics/disjointness/external-clusters/${SLURM_JOB_ID}"

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

# Fail closed on an incomplete or missing source run rather than clustering
# against a partial sketch.
if [[ ! -f "${SOURCE_ROOT}/COMPLETE" ]]; then
    printf 'Source comparison is not COMPLETE: %s\n' "${SOURCE_ROOT}" >&2
    exit 2
fi
for required in external-assemblies.txt assembly-path-map.tsv artifact.sha256; do
    if [[ ! -e "${SOURCE_ROOT}/${required}" ]]; then
        printf 'Source comparison is missing %s\n' "${required}" >&2
        exit 2
    fi
done
if [[ -e "${RESULT_ROOT}" ]]; then
    printf 'Refusing to overwrite cluster result: %s\n' "${RESULT_ROOT}" >&2
    exit 2
fi
mkdir -p "${RESULT_ROOT}/tmp"
export TMPDIR="${resolved_project}/${RESULT_ROOT}/tmp"

# Reauthenticate the reused inputs against the source run's own checksums.
(cd "${SOURCE_ROOT}" && sha256sum --check --ignore-missing artifact.sha256) \
    > "${RESULT_ROOT}/source-verification.txt"

printf 'source_comparison=%s\n' "${AMR_SOURCE_COMPARISON}" > "${RESULT_ROOT}/provenance.txt"
{
    printf 'requests_sha256=%s\n' "$(sha256sum "${REQUESTS}" | awk '{print $1}')"
    printf 'contract_sha256=%s\n' "$(sha256sum "${CONTRACT}" | awk '{print $1}')"
    printf 'path_map_sha256=%s\n' "$(sha256sum "${SOURCE_ROOT}/assembly-path-map.tsv" | awk '{print $1}')"
    printf 'skani_version=%s\n' "$("${DISJOINTNESS_ENV}/bin/skani" --version)"
    printf 'skani_search_parameters=-s 95 --both-min-af 50 --keep-refs\n'
    printf 'sketch_built_here=external-sketch\n'
    printf 'external_list=%s/external-assemblies.txt\n' "${SOURCE_ROOT}"
} >> "${RESULT_ROOT}/provenance.txt"

# No external sketch exists: 63772662 sketched development only. Build one
# from the same frozen assembly list, with the frozen sketch parameters.
"${DISJOINTNESS_ENV}/bin/skani" sketch \
    -l "${SOURCE_ROOT}/external-assemblies.txt" \
    -o "${RESULT_ROOT}/external-sketch" \
    -t "${SLURM_CPUS_PER_TASK}" -c 125 -m 1000

"${DISJOINTNESS_ENV}/bin/skani" search \
    -d "${RESULT_ROOT}/external-sketch" \
    --ql "${SOURCE_ROOT}/external-assemblies.txt" \
    -o "${RESULT_ROOT}/skani-search.raw.tsv" \
    -t "${SLURM_CPUS_PER_TASK}" -s 95 --both-min-af 50 --keep-refs

"${CONTROLLER_ENV}/bin/python" scripts/build_genomic_disjointness.py normalize-skani \
    --raw "${RESULT_ROOT}/skani-search.raw.tsv" \
    --path-map "${SOURCE_ROOT}/assembly-path-map.tsv" \
    --query-cohort external \
    --reference-cohort external \
    --output "${RESULT_ROOT}/comparisons.tsv" \
    --manifest "${RESULT_ROOT}/normalization-manifest.json"

"${CONTROLLER_ENV}/bin/python" scripts/build_genomic_disjointness.py cluster-development \
    --requests "${REQUESTS}" \
    --contract "${CONTRACT}" \
    --comparisons "${RESULT_ROOT}/comparisons.tsv" \
    --clusters "${RESULT_ROOT}/external-clusters.tsv" \
    --manifest "${RESULT_ROOT}/cluster-manifest.json" \
    --cohort external

(
    cd "${RESULT_ROOT}"
    sha256sum provenance.txt source-verification.txt skani-search.raw.tsv comparisons.tsv \
        normalization-manifest.json external-clusters.tsv cluster-manifest.json \
        > artifact.sha256
)
rm -rf "${RESULT_ROOT}/tmp"
touch "${RESULT_ROOT}/COMPLETE"
