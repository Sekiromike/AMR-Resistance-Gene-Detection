#!/usr/bin/env bash
set -euo pipefail

: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT to your private Unity /work allocation path}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENV_PREFIX="${AMR_CONTROLLER_ENV:-${WORK_ROOT}/conda/envs/amr-ecoli-controller}"

cd "${PROJECT_ROOT}"
mkdir -p logs provenance
if ! type module >/dev/null 2>&1; then
    set +u
    source /etc/profile.d/z00-lmod-profile.sh
    source /etc/profile.d/z03-lmod-vars.sh
    source /etc/profile.d/z06-lmod-modulepath.sh
    set -u
fi
module load conda/latest
export CONDA_PKGS_DIRS="${WORK_ROOT}/conda/pkgs"
exec "${ENV_PREFIX}/bin/snakemake" \
    --snakefile workflow/Snakefile \
    --configfile config/study.json \
    --workflow-profile profiles/unity \
    --conda-prefix "${WORK_ROOT}/conda/snakemake" \
    "$@"
