#!/usr/bin/env bash
# Validate the exact skani sketch/search output contract on one frozen pair.
#SBATCH --job-name=amr-skani-canary
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/disjointness-search-canary-%j.out
#SBATCH --error=logs/disjointness-search-canary-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
CONTROLLER_ENV="${AMR_CONTROLLER_ENV:-${WORK_ROOT}/conda/envs/amr-ecoli-controller}"
ENV_PREFIX="${AMR_DISJOINTNESS_ENV:-${WORK_ROOT}/conda/envs/amr-genomic-disjointness}"
DEVELOPMENT_ASSEMBLY="${AMR_CANARY_DEVELOPMENT_ASSEMBLY:-results/genomics/disjointness/assemblies/development/PDT000041778.1/assembly.fna}"
EXTERNAL_ASSEMBLY="${AMR_CANARY_EXTERNAL_ASSEMBLY:-results/genomics/disjointness/assemblies/external/JBABADA-19-0001/assembly.fna}"
DEVELOPMENT_SEQUENCE_ID="${AMR_CANARY_DEVELOPMENT_SEQUENCE_ID:-development::PDT000041778.1}"
EXTERNAL_SEQUENCE_ID="${AMR_CANARY_EXTERNAL_SEQUENCE_ID:-external::JBABADA-19-0001}"
case "${DEVELOPMENT_SEQUENCE_ID}|${EXTERNAL_SEQUENCE_ID}" in
    development::*'|'external::*) ;;
    *)
        printf 'Canary sequence IDs do not preserve cohort prefixes\n' >&2
        exit 2
        ;;
esac
DEVELOPMENT_ISOLATE_ID="${DEVELOPMENT_SEQUENCE_ID#development::}"
EXTERNAL_ISOLATE_ID="${EXTERNAL_SEQUENCE_ID#external::}"
RESULT_ROOT="results/genomics/disjointness/search-canary/${SLURM_JOB_ID}"

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
for assembly in "${DEVELOPMENT_ASSEMBLY}" "${EXTERNAL_ASSEMBLY}"; do
    resolved_assembly="$(realpath -m "${assembly}")"
    case "${resolved_assembly}" in
        "${resolved_project}"/*) ;;
        *)
            printf 'Canary assembly is outside the project: %s\n' "${resolved_assembly}" >&2
            exit 2
            ;;
    esac
    if [[ ! -s "${resolved_assembly}" || ! -f "$(dirname "${resolved_assembly}")/COMPLETE" ]]; then
        printf 'Canary assembly is not complete: %s\n' "${resolved_assembly}" >&2
        exit 2
    fi
done
if [[ ! -x "${ENV_PREFIX}/bin/skani" || ! -x "${CONTROLLER_ENV}/bin/python" || -e "${RESULT_ROOT}" ]]; then
    printf 'Verified environment is unavailable or canary result already exists\n' >&2
    exit 2
fi
mkdir -p "${RESULT_ROOT}"

"${ENV_PREFIX}/bin/skani" sketch \
    -o "${RESULT_ROOT}/development-sketch" \
    -t "${SLURM_CPUS_PER_TASK}" -c 125 -m 1000 \
    "${DEVELOPMENT_ASSEMBLY}"
"${ENV_PREFIX}/bin/skani" search \
    -d "${RESULT_ROOT}/development-sketch" \
    -o "${RESULT_ROOT}/skani-search.raw.tsv" \
    -t "${SLURM_CPUS_PER_TASK}" -s 95 --both-min-af 50 --keep-refs \
    "${EXTERNAL_ASSEMBLY}"
{
    printf 'sequence_id\tcohort\tisolate_id\tassembly_path\n'
    printf '%s\tdevelopment\t%s\t%s\n' \
        "${DEVELOPMENT_SEQUENCE_ID}" "${DEVELOPMENT_ISOLATE_ID}" "${DEVELOPMENT_ASSEMBLY}"
    printf '%s\texternal\t%s\t%s\n' \
        "${EXTERNAL_SEQUENCE_ID}" "${EXTERNAL_ISOLATE_ID}" "${EXTERNAL_ASSEMBLY}"
} > "${RESULT_ROOT}/assembly-path-map.tsv"
"${CONTROLLER_ENV}/bin/python" scripts/build_genomic_disjointness.py normalize-skani \
    --raw "${RESULT_ROOT}/skani-search.raw.tsv" \
    --path-map "${RESULT_ROOT}/assembly-path-map.tsv" \
    --output "${RESULT_ROOT}/comparisons.tsv" \
    --manifest "${RESULT_ROOT}/normalization-manifest.json"
{
    printf 'development_assembly=%s\n' "${DEVELOPMENT_ASSEMBLY}"
    printf 'external_assembly=%s\n' "${EXTERNAL_ASSEMBLY}"
    printf 'skani_version=%s\n' "$("${ENV_PREFIX}/bin/skani" --version)"
    printf 'skani_sketch_parameters=-c 125 -m 1000\n'
    printf 'skani_search_parameters=-s 95 --both-min-af 50 --keep-refs\n'
} > "${RESULT_ROOT}/provenance.txt"
(
    cd "${RESULT_ROOT}"
    sha256sum skani-search.raw.tsv assembly-path-map.tsv comparisons.tsv \
        normalization-manifest.json provenance.txt > artifact.sha256
)
touch "${RESULT_ROOT}/COMPLETE"
