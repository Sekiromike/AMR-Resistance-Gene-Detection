#!/usr/bin/env bash
#SBATCH --job-name=amr-ecoli-smoke
#SBATCH --partition=cpu-preempt
#SBATCH --qos=short
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/unity-smoke-%j.out

set -euo pipefail

: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT before submitting this job}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENV_PREFIX="${AMR_CONTROLLER_ENV:-${WORK_ROOT}/conda/envs/amr-ecoli-controller}"

cd "${PROJECT_ROOT}"
if ! type module >/dev/null 2>&1; then
    set +u
    source /etc/profile.d/z00-lmod-profile.sh
    source /etc/profile.d/z03-lmod-vars.sh
    source /etc/profile.d/z06-lmod-modulepath.sh
    set -u
fi
module load conda/latest
export CONDA_PKGS_DIRS="${WORK_ROOT}/conda/pkgs"
"${ENV_PREFIX}/bin/python" -m compileall -q scripts tests
"${ENV_PREFIX}/bin/python" -m unittest discover -s tests -v
"${ENV_PREFIX}/bin/snakemake" \
    --snakefile workflow/Snakefile \
    --configfile workflow/testdata/config.json \
    --executor local \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --software-deployment-method conda \
    --conda-prefix "${WORK_ROOT}/conda/snakemake" \
    --forceall \
    --printshellcmds \
    --rerun-incomplete
