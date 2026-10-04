#!/usr/bin/env bash
# Build the modelling environment (workflow/envs/models.yaml) on a compute
# node, never the login node; package cache stays in the allocation. Runs on
# the non-preemptible cpu partition so the solve cannot be killed midway.
#SBATCH --job-name=amr-models-env
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/models-env-%j.out
#SBATCH --error=logs/models-env-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENV_PREFIX="${AMR_MODELS_ENV:-${WORK_ROOT}/conda/envs/amr-models}"
PKG_CACHE="${AMR_CONDA_PKGS:-${WORK_ROOT}/conda/pkgs}"
SPEC="${PROJECT_ROOT}/workflow/envs/models.yaml"
STAMP="${ENV_PREFIX}/.amr-models-spec.sha256"

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
if [[ ! -x "${ENV_PREFIX}/bin/python" ]]; then
    conda env create --prefix "${ENV_PREFIX}" --file "${SPEC}"
elif [[ "$(cat "${STAMP}" 2>/dev/null || true)" != "${expected_hash}" ]]; then
    conda env update --prefix "${ENV_PREFIX}" --file "${SPEC}" --prune
fi
printf '%s\n' "${expected_hash}" > "${STAMP}"
conda list --prefix "${ENV_PREFIX}" --explicit > provenance/models-conda-explicit.txt
"${ENV_PREFIX}/bin/python" -c 'import xgboost, sklearn, scipy, numpy, pandas; print("xgboost", xgboost.__version__, "sklearn", sklearn.__version__, "scipy", scipy.__version__, "numpy", numpy.__version__, "pandas", pandas.__version__)' \
    > provenance/models-tool-versions.txt
cat provenance/models-tool-versions.txt
sha256sum "${SPEC}" provenance/models-conda-explicit.txt provenance/models-tool-versions.txt \
    > provenance/models-environment.sha256
