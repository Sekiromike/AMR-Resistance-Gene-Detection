#!/usr/bin/env bash
# Build the GPU embedding environment. This MUST run as a job on a GPU node,
# not on a login node: Unity prohibits heavy work on login nodes, and conda's
# __cuda virtual package is only detected where a GPU is present, so a CUDA
# build of pytorch cannot be solved on a login node at all.
#SBATCH --job-name=amr-embed-env
#SBATCH --partition=gpu-preempt
#SBATCH --gpus=1
#SBATCH --constraint=vram16
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --output=logs/embeddings-env-%j.out
#SBATCH --error=logs/embeddings-env-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in config/unity.local.env}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENV_PREFIX="${AMR_EMBEDDING_ENV:-${WORK_ROOT}/conda/envs/amr-genome-embeddings}"
PKG_CACHE="${AMR_CONDA_PKGS:-${WORK_ROOT}/conda/pkgs}"
SPEC="${PROJECT_ROOT}/workflow/envs/embeddings.yaml"
STAMP="${ENV_PREFIX}/.amr-embeddings-spec.sha256"

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
# Keep the package cache in the work allocation. The default ~/.conda/pkgs
# consumes the 100 GB home quota and is only needed at build time.
export CONDA_PKGS_DIRS="${PKG_CACHE}"

expected_hash="$(sha256sum "${SPEC}" | awk '{print $1}')"
observed_hash="$(cat "${STAMP}" 2>/dev/null || true)"
if [[ ! -x "${ENV_PREFIX}/bin/python" ]]; then
    conda env create --prefix "${ENV_PREFIX}" --file "${SPEC}"
elif [[ "${observed_hash}" != "${expected_hash}" ]]; then
    conda env update --prefix "${ENV_PREFIX}" --file "${SPEC}" --prune
fi
printf '%s\n' "${expected_hash}" > "${STAMP}"
conda list --prefix "${ENV_PREFIX}" --explicit \
    > provenance/embeddings-conda-explicit.txt

# Record the device the environment was solved against, because the CUDA build
# selected by conda depends on the driver present at solve time.
{
    printf 'node=%s\n' "${SLURMD_NODENAME}"
    nvidia-smi --query-gpu=name,compute_cap,memory.total,driver_version \
        --format=csv,noheader | sed 's/^/gpu=/'
    "${ENV_PREFIX}/bin/python" --version
    "${ENV_PREFIX}/bin/python" -c 'import torch; print("torch", torch.__version__)'
    "${ENV_PREFIX}/bin/python" -c 'import torch; print("cuda_available", torch.cuda.is_available())'
    "${ENV_PREFIX}/bin/python" -c 'import torch; print("cuda_version", torch.version.cuda)'
    "${ENV_PREFIX}/bin/python" -c 'import transformers; print("transformers", transformers.__version__)'
} > provenance/embeddings-tool-versions.txt 2>&1

sha256sum "${SPEC}" provenance/embeddings-conda-explicit.txt \
    provenance/embeddings-tool-versions.txt \
    > provenance/embeddings-environment.sha256
