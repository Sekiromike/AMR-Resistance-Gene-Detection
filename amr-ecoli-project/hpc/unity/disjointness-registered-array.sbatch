#!/usr/bin/env bash
# Submit with: sbatch --array=1-100%5 --export=ALL hpc/unity/disjointness-registered-array.sbatch
#SBATCH --job-name=amr-disjoint-ref
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --output=logs/disjointness-ref-%A_%a.out
#SBATCH --error=logs/disjointness-ref-%A_%a.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${SLURM_ARRAY_TASK_ID:?Submit this script as an array}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENV_PREFIX="${AMR_DISJOINTNESS_ENV:-${WORK_ROOT}/conda/envs/amr-genomic-disjointness}"
REQUESTS="${PROJECT_ROOT}/data/genomics/disjointness/requests/sequence_requests.tsv"
ASSEMBLY_ROOT="${PROJECT_ROOT}/results/genomics/disjointness/assemblies/development"
BATCH_SIZE="${AMR_REGISTERED_BATCH_SIZE:-100}"
export PATH="${ENV_PREFIX}/bin:${PATH}"
if [[ ! -x "${ENV_PREFIX}/bin/datasets" ]]; then
    printf 'Verified disjointness environment is unavailable: %s\n' "${ENV_PREFIX}" >&2
    exit 2
fi

cd "${PROJECT_ROOT}"
start=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))
end=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE ))
work="${PROJECT_ROOT}/results/genomics/disjointness/registered-batches/${SLURM_ARRAY_JOB_ID}/${SLURM_ARRAY_TASK_ID}"
mkdir -p "${work}"
awk -F '\t' -v start="${start}" -v end="${end}" \
    'NR>1 && $5=="registered_assembly" {n++; if(n>=start && n<=end) print $6}' \
    "${REQUESTS}" > "${work}/accessions.txt"
if [[ ! -s "${work}/accessions.txt" ]]; then
    printf 'No registered-assembly requests in batch %s\n' "${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
fi

"${ENV_PREFIX}/bin/datasets" download genome accession \
    --inputfile "${work}/accessions.txt" --include genome --dehydrated \
    --filename "${work}/package.zip"
mkdir -p "${work}/package"
unzip -q "${work}/package.zip" -d "${work}/package"
"${ENV_PREFIX}/bin/datasets" rehydrate --directory "${work}/package"

while IFS= read -r accession; do
    isolate="$(awk -F '\t' -v accession="${accession}" \
        'NR>1 && $5=="registered_assembly" && $6==accession {print $3; exit}' "${REQUESTS}")"
    if [[ -z "${isolate}" ]]; then
        printf 'Frozen mapping missing for %s\n' "${accession}" >&2
        exit 3
    fi
    source_dir="${work}/package/ncbi_dataset/data/${accession}"
    fasta="$(find "${source_dir}" -maxdepth 1 -type f -name '*_genomic.fna' -print -quit)"
    if [[ -z "${fasta}" ]]; then
        printf 'Official package missing genomic FASTA for %s\n' "${accession}" >&2
        exit 4
    fi
    result="${ASSEMBLY_ROOT}/${isolate}"
    if [[ -f "${result}/COMPLETE" ]]; then
        (cd "${result}" && sha256sum --check artifact.sha256)
        continue
    fi
    if [[ -e "${result}" ]]; then
        printf 'Incomplete result path requires review before retry: %s\n' "${result}" >&2
        exit 5
    fi
    mkdir -p "${result}"
    cp "${fasta}" "${result}/assembly.fna"
    {
        printf 'completed_at_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'slurm_array_job_id=%s\n' "${SLURM_ARRAY_JOB_ID}"
        printf 'slurm_array_task_id=%s\n' "${SLURM_ARRAY_TASK_ID}"
        printf 'isolate_id=%s\n' "${isolate}"
        printf 'assembly_accession=%s\n' "${accession}"
        printf 'assembly_sha256=%s\n' "$(sha256sum "${result}/assembly.fna" | awk '{print $1}')"
    } > "${result}/provenance.txt"
    (cd "${result}" && sha256sum assembly.fna provenance.txt > artifact.sha256)
    touch "${result}/COMPLETE"
done < "${work}/accessions.txt"
touch "${work}/BATCH_COMPLETE"
