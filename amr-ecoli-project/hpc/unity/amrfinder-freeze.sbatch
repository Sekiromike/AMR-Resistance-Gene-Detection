#!/usr/bin/env bash
# Build the AMRFinderPlus environment and discover the current
# database release, so config/study.json can pin it.
#
# quast, mlst and amrfinder cannot share one environment: their perl, zlib and
# python_abi requirements are mutually unsatisfiable. Each has its own spec.
#
# This probes into a scratch location and does NOT write the frozen database.
# Freezing is the job of the genomics_freeze_amrfinder_database rule, which
# verifies that the release it downloads equals the pinned version and fails
# closed on drift.
#
# Runs on the non-preemptible cpu partition because it combines a Conda solve
# with a network download and must not be killed midway.
#SBATCH --job-name=amr-amrfinder-freeze
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/amrfinder-freeze-%j.out
#SBATCH --error=logs/amrfinder-freeze-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENV_PREFIX="${AMR_GENOMICS_ENV:-${WORK_ROOT}/conda/envs/amr-genomics-amrfinder}"
PKG_CACHE="${AMR_CONDA_PKGS:-${WORK_ROOT}/conda/pkgs}"
SPEC="${PROJECT_ROOT}/workflow/envs/genomics-amrfinder.yaml"
STAMP="${ENV_PREFIX}/.amr-amrfinder-spec.sha256"
PROBE="${PROJECT_ROOT}/results/genomics/amrfinder-probe/${SLURM_JOB_ID}"

cd "${PROJECT_ROOT}"
mkdir -p logs provenance "$(dirname "${ENV_PREFIX}")" "${PKG_CACHE}" "${PROBE}"
if ! type module >/dev/null 2>&1; then
    set +u
    source /etc/profile.d/z00-lmod-profile.sh
    source /etc/profile.d/z03-lmod-vars.sh
    source /etc/profile.d/z06-lmod-modulepath.sh
    set -u
fi
module load conda/latest
# Never let the package cache fall back to the 100 GB home quota.
export CONDA_PKGS_DIRS="${PKG_CACHE}"

expected_hash="$(sha256sum "${SPEC}" | awk '{print $1}')"
observed_hash="$(cat "${STAMP}" 2>/dev/null || true)"
if [[ ! -x "${ENV_PREFIX}/bin/amrfinder" ]]; then
    conda env create --prefix "${ENV_PREFIX}" --file "${SPEC}"
elif [[ "${observed_hash}" != "${expected_hash}" ]]; then
    conda env update --prefix "${ENV_PREFIX}" --file "${SPEC}" --prune
fi
printf '%s\n' "${expected_hash}" > "${STAMP}"
conda list --prefix "${ENV_PREFIX}" --explicit > provenance/amrfinder-conda-explicit.txt

export PATH="${ENV_PREFIX}/bin:${PATH}"

# Probe the current release into a scratch directory.
amrfinder_update --force_update --database "${PROBE}/database"
resolved="$(readlink -f "${PROBE}/database/latest")"
version="$(head -n 1 "${resolved}/version.txt" | tr -d '\r\n')"

{
    printf 'job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'node=%s\n' "${SLURMD_NODENAME}"
    printf 'amrfinder_version=%s\n' "$(amrfinder --version 2>&1 | head -1)"
    printf 'amrfinder_database_version=%s\n' "${version}"
    printf 'database_path=%s\n' "${resolved}"
} | tee "provenance/amrfinder-database-probe-${SLURM_JOB_ID}.txt"

amrfinder --database "${resolved}" --database_version \
    > "provenance/amrfinder-database-metadata-${SLURM_JOB_ID}.txt" 2>&1 || true

sha256sum "${SPEC}" provenance/amrfinder-conda-explicit.txt \
    "provenance/amrfinder-database-probe-${SLURM_JOB_ID}.txt" \
    > "provenance/genomics-environment-${SLURM_JOB_ID}.sha256"

printf '\nPIN THIS IN config/study.json genomics.amrfinder_database_version: %s\n' "${version}"
