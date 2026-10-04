#!/usr/bin/env bash
set -euo pipefail

: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT to your private Unity /work allocation path}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENV_PREFIX="${AMR_CONTROLLER_ENV:-${WORK_ROOT}/conda/envs/amr-ecoli-controller}"
PKG_CACHE="${AMR_CONDA_PKGS:-${WORK_ROOT}/conda/pkgs}"
ENV_SPEC="${PROJECT_ROOT}/workflow/envs/controller.yaml"
ENV_SPEC_STAMP="${ENV_PREFIX}/.amr-controller-spec.sha256"

if [[ ! -d "${PROJECT_ROOT}" ]]; then
    printf 'Project directory not found: %s\n' "${PROJECT_ROOT}" >&2
    exit 2
fi

mkdir -p "$(dirname "${ENV_PREFIX}")" "${PKG_CACHE}" "${PROJECT_ROOT}/logs" "${PROJECT_ROOT}/provenance"
if ! type module >/dev/null 2>&1; then
    # Noninteractive SSH does not source Unity's Lmod profile scripts.
    set +u
    source /etc/profile.d/z00-lmod-profile.sh
    source /etc/profile.d/z03-lmod-vars.sh
    source /etc/profile.d/z06-lmod-modulepath.sh
    set -u
fi
module load conda/latest
export CONDA_PKGS_DIRS="${PKG_CACHE}"

expected_spec_hash="$(sha256sum "${ENV_SPEC}" | awk '{print $1}')"
observed_spec_hash="$(cat "${ENV_SPEC_STAMP}" 2>/dev/null || true)"
if [[ ! -x "${ENV_PREFIX}/bin/snakemake" ]]; then
    conda env create --prefix "${ENV_PREFIX}" --file "${ENV_SPEC}"
elif [[ "${observed_spec_hash}" != "${expected_spec_hash}" ]]; then
    conda env update --prefix "${ENV_PREFIX}" --file "${ENV_SPEC}" --prune
fi
printf '%s\n' "${expected_spec_hash}" > "${ENV_SPEC_STAMP}"

conda list --prefix "${ENV_PREFIX}" --explicit > "${PROJECT_ROOT}/provenance/controller-conda-explicit.txt"
"${ENV_PREFIX}/bin/python" --version
"${ENV_PREFIX}/bin/snakemake" --version
printf 'Controller environment: %s\n' "${ENV_PREFIX}"
