#!/usr/bin/env bash
#SBATCH --job-name=amr-gpu-smoke
#SBATCH --partition=gpu-preempt
#SBATCH --time=00:20:00
#SBATCH --gpus=1
#SBATCH --constraint=sm_70&vram16
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/gpu-smoke-%j.out

set -euo pipefail

: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT before submitting this job}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"

cd "${PROJECT_ROOT}"
mkdir -p logs provenance

{
    printf 'job_id\t%s\n' "${SLURM_JOB_ID}"
    printf 'node\t%s\n' "${SLURMD_NODENAME}"
    printf 'allocated_gpus\t%s\n' "${SLURM_GPUS_ON_NODE:-unknown}"
    printf 'cuda_visible_devices\t%s\n' "${CUDA_VISIBLE_DEVICES:-unknown}"
    nvidia-smi --query-gpu=name,uuid,compute_cap,memory.total,driver_version \
        --format=csv,noheader
} | tee "provenance/gpu-smoke-${SLURM_JOB_ID}.txt"

nvidia-smi
