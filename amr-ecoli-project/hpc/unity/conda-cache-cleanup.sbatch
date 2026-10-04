#!/usr/bin/env bash
# Reclaim home-quota space by removing the Conda package cache only.
#
# The cache is needed at build time and can be recreated. Environments,
# credentials, and the environments registry are NOT caches and are preserved
# explicitly. This runs as a job because walking and deleting a large tree is
# exactly the login-node work Unity prohibits.
#SBATCH --job-name=amr-conda-cache-cleanup
#SBATCH --partition=cpu-preempt
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/conda-cache-cleanup-%j.out
#SBATCH --error=logs/conda-cache-cleanup-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
REPORT_DIR="${PROJECT_ROOT}/provenance"
CACHE="${HOME}/.conda/pkgs"

cd "${PROJECT_ROOT}"
mkdir -p logs "${REPORT_DIR}"
REPORT="${REPORT_DIR}/conda-cache-cleanup-${SLURM_JOB_ID}.txt"

{
    printf 'job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'node=%s\n' "${SLURMD_NODENAME}"
    printf 'home=%s\n' "${HOME}"

    printf '\n--- before ---\n'
    du -sh "${HOME}/.conda" 2>/dev/null || true
    for entry in "${HOME}"/.conda/*; do
        [[ -e "${entry}" ]] || continue
        du -sh "${entry}" 2>/dev/null || true
    done

    printf '\n--- preserved (never deleted) ---\n'
    for keep in envs environments.txt aau_token aau_token_host tos; do
        if [[ -e "${HOME}/.conda/${keep}" ]]; then
            printf 'preserve %s\n' "${keep}"
        fi
    done
    printf '\n--- environments registry ---\n'
    cat "${HOME}/.conda/environments.txt" 2>/dev/null || printf '(none)\n'
} | tee "${REPORT}"

# Refuse to proceed if the cache path is anything other than the expected cache.
if [[ "${CACHE}" != "${HOME}/.conda/pkgs" ]]; then
    printf 'Refusing to delete an unexpected path: %s\n' "${CACHE}" >&2
    exit 2
fi
if [[ -e "${HOME}/.conda/pkgs/envs" ]]; then
    printf 'Refusing to delete: an envs directory exists inside the cache path\n' >&2
    exit 2
fi

if [[ -d "${CACHE}" ]]; then
    rm -rf "${CACHE}"
    printf 'deleted %s\n' "${CACHE}" | tee -a "${REPORT}"
else
    printf 'cache absent, nothing to delete\n' | tee -a "${REPORT}"
fi

{
    printf '\n--- after ---\n'
    du -sh "${HOME}/.conda" 2>/dev/null || true
    printf '\n--- environments still present ---\n'
    ls -1 "${HOME}/.conda/envs" 2>/dev/null || printf '(none)\n'
} | tee -a "${REPORT}"

sha256sum "${REPORT}" > "${REPORT}.sha256"
