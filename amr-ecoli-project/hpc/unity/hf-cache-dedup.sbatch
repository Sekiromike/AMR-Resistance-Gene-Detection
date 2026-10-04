#!/usr/bin/env bash
# Remove a HuggingFace model copy from the home quota when a byte-identical
# copy already exists in the work allocation.
#
# Fails closed: the home copy is deleted only after the work copy is confirmed
# present, complete, and identical in size. Runs as a job because walking and
# deleting a large tree is login-node work that Unity prohibits.
#SBATCH --job-name=amr-hf-dedup
#SBATCH --partition=cpu-preempt
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:45:00
#SBATCH --output=logs/hf-dedup-%j.out
#SBATCH --error=logs/hf-dedup-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_DEDUP_MODEL:?Set AMR_DEDUP_MODEL, e.g. models--arcinstitute--evo2_7b}"
: "${AMR_DEDUP_WEIGHT:?Set AMR_DEDUP_WEIGHT, the weight filename, e.g. evo2_7b.pt}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
HOME_COPY="${HOME}/.cache/huggingface/hub/${AMR_DEDUP_MODEL}"
WORK_COPY="${AMR_DEDUP_WORK_COPY:-${WORK_ROOT}/evo2_project/cache/hf/hub/${AMR_DEDUP_MODEL}}"

cd "${PROJECT_ROOT}"
mkdir -p logs provenance
REPORT="provenance/hf-dedup-${SLURM_JOB_ID}.txt"

weight_size() {
    find "$1" -name "${AMR_DEDUP_WEIGHT}" -exec stat -Lc %s {} \; 2>/dev/null | head -1
}

{
    printf 'job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'model=%s\n' "${AMR_DEDUP_MODEL}"
    printf 'home_copy=%s\n' "${HOME_COPY}"
    printf 'work_copy=%s\n' "${WORK_COPY}"
} | tee "${REPORT}"

if [[ ! -d "${HOME_COPY}" ]]; then
    printf 'home copy absent, nothing to do\n' | tee -a "${REPORT}"
    exit 0
fi
if [[ ! -d "${WORK_COPY}" ]]; then
    printf 'REFUSING: work copy does not exist\n' | tee -a "${REPORT}" >&2
    exit 2
fi

home_bytes="$(weight_size "${HOME_COPY}")"
work_bytes="$(weight_size "${WORK_COPY}")"
incomplete="$(find "${WORK_COPY}" -name '*.incomplete' 2>/dev/null | wc -l)"

{
    printf 'home_bytes=%s\n' "${home_bytes:-missing}"
    printf 'work_bytes=%s\n' "${work_bytes:-missing}"
    printf 'work_incomplete_blobs=%s\n' "${incomplete}"
} | tee -a "${REPORT}"

if [[ -z "${work_bytes}" || "${home_bytes}" != "${work_bytes}" ]]; then
    printf 'REFUSING: work copy is missing or differs in size\n' | tee -a "${REPORT}" >&2
    exit 2
fi
if [[ "${incomplete}" != "0" ]]; then
    printf 'REFUSING: work copy has incomplete blobs\n' | tee -a "${REPORT}" >&2
    exit 2
fi

du -sh "${HOME_COPY}" 2>/dev/null | sed 's/^/before /' | tee -a "${REPORT}"
rm -rf "${HOME_COPY}"
printf 'deleted %s\n' "${HOME_COPY}" | tee -a "${REPORT}"

{
    printf '\n--- home usage after ---\n'
    du -sh "${HOME}/.cache" "${HOME}/.conda" 2>/dev/null || true
} | tee -a "${REPORT}"

sha256sum "${REPORT}" > "${REPORT}.sha256"
