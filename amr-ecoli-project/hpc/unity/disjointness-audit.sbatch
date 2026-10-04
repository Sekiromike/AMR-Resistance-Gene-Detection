#!/usr/bin/env bash
#SBATCH --job-name=amr-disjoint-audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/disjointness-audit-%j.out
#SBATCH --error=logs/disjointness-audit-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENV_PREFIX="${AMR_CONTROLLER_ENV:-${WORK_ROOT}/conda/envs/amr-ecoli-controller}"

cd "${PROJECT_ROOT}"
"${ENV_PREFIX}/bin/python" scripts/build_genomic_disjointness.py audit-assemblies \
    --requests data/genomics/disjointness/requests/sequence_requests.tsv \
    --assemblies-root results/genomics/disjointness/assemblies \
    --manifest results/genomics/disjointness/assembly-audit-manifest.json \
    --workers "${SLURM_CPUS_PER_TASK}"
