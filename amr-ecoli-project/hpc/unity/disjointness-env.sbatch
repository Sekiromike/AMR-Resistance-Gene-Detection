#!/usr/bin/env bash
#SBATCH --job-name=amr-disjoint-env
#SBATCH --partition=cpu-preempt
#SBATCH --qos=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/disjointness-env-%j.out
#SBATCH --error=logs/disjointness-env-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in config/unity.local.env}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENV_PREFIX="${AMR_DISJOINTNESS_ENV:-${WORK_ROOT}/conda/envs/amr-genomic-disjointness}"
PKG_CACHE="${AMR_CONDA_PKGS:-${WORK_ROOT}/conda/pkgs}"
SPEC="${PROJECT_ROOT}/workflow/envs/disjointness.yaml"
STAMP="${ENV_PREFIX}/.amr-disjointness-spec.sha256"

cd "${PROJECT_ROOT}"
mkdir -p logs provenance "$(dirname "${ENV_PREFIX}")" "${PKG_CACHE}"
if ! type module >/dev/null 2>&1; then
    set +u
    source /etc/profile.d/z00-lmod-profile.sh
    source /etc/profile.d/z03-lmod-vars.sh
    source /etc/profile.d/z06-lmod-modulepath.sh
    set -u
fi
module load conda/latest
export CONDA_PKGS_DIRS="${PKG_CACHE}"

expected_hash="$(sha256sum "${SPEC}" | awk '{print $1}')"
observed_hash="$(cat "${STAMP}" 2>/dev/null || true)"
if [[ ! -x "${ENV_PREFIX}/bin/skani" ]]; then
    conda env create --prefix "${ENV_PREFIX}" --file "${SPEC}"
elif [[ "${observed_hash}" != "${expected_hash}" ]]; then
    conda env update --prefix "${ENV_PREFIX}" --file "${SPEC}" --prune
fi
printf '%s\n' "${expected_hash}" > "${STAMP}"
conda list --prefix "${ENV_PREFIX}" --explicit \
    > provenance/disjointness-conda-explicit.txt
export PATH="${ENV_PREFIX}/bin:${PATH}"
{
    datasets --version
    prefetch --version
    fasterq-dump --version
    shovill --version
    skani --version
} > provenance/disjointness-tool-versions.txt 2>&1
sha256sum "${SPEC}" provenance/disjointness-conda-explicit.txt \
    provenance/disjointness-tool-versions.txt \
    > provenance/disjointness-environment.sha256
