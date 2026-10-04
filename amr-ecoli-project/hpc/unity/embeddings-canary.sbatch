#!/usr/bin/env bash
# Bounded foundation-model embedding canary. Measures peak VRAM, throughput and
# reverse-complement consistency on a handful of assemblies. No production array
# may be submitted before this passes. See docs/FOUNDATION_MODEL_PROTOCOL.md.
#SBATCH --job-name=amr-embed-canary
#SBATCH --partition=gpu-preempt
#SBATCH --time=01:30:00
#SBATCH --gpus=1
#SBATCH --constraint=vram143
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=logs/embeddings-canary-%j.out
#SBATCH --error=logs/embeddings-canary-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
EMBEDDING_ENV="${AMR_EMBEDDING_ENV:-${HOME}/.conda/envs/evo2}"
ASSEMBLIES="results/genomics/disjointness/assemblies"
RESULT_ROOT="results/genomics/embeddings/canary/${SLURM_JOB_ID}"
CANARY_COUNT="${AMR_CANARY_COUNT:-3}"

resolved_work="$(realpath -m "${WORK_ROOT}")"
resolved_project="$(realpath -m "${PROJECT_ROOT}")"
case "${resolved_project}" in
    "${resolved_work}"/*) ;;
    *)
        printf 'Project root is outside AMR_WORK_ROOT: %s\n' "${resolved_project}" >&2
        exit 2
        ;;
esac
cd "${resolved_project}"
if [[ ! -x "${EMBEDDING_ENV}/bin/python" ]]; then
    printf 'Embedding environment is unavailable: %s\n' "${EMBEDDING_ENV}" >&2
    exit 2
fi
if [[ ! -d "${ASSEMBLIES}" ]]; then
    printf 'Assembly root is unavailable: %s\n' "${ASSEMBLIES}" >&2
    exit 2
fi
if [[ -e "${RESULT_ROOT}" ]]; then
    printf 'Refusing to overwrite canary result: %s\n' "${RESULT_ROOT}" >&2
    exit 2
fi
mkdir -p "${RESULT_ROOT}"
export TMPDIR="${resolved_project}/${RESULT_ROOT}"
# Keep the model cache inside the allocation, never in the 100 GB home directory.
export HF_HOME="${WORK_ROOT}/huggingface"
export HF_HUB_OFFLINE="${AMR_HF_OFFLINE:-0}"
# transformer_engine._load_nvrtc globs CUDA_HOME first and otherwise falls back
# to `ldconfig -p | grep libnvrtc`, which raises on Unity because CUDA arrives
# through modules that set LD_LIBRARY_PATH without populating the ldconfig
# cache. Point CUDA_HOME at the libnvrtc bundled in the environment.
export CUDA_HOME="${AMR_CUDA_HOME:-${EMBEDDING_ENV}/lib/python3.12/site-packages/nvidia}"
if [[ ! -e "${CUDA_HOME}/cuda_nvrtc/lib/libnvrtc.so.12" ]]; then
    printf 'libnvrtc not found under CUDA_HOME: %s\n' "${CUDA_HOME}" >&2
    exit 2
fi

{
    printf 'job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'node=%s\n' "${SLURMD_NODENAME}"
    printf 'canary_count=%s\n' "${CANARY_COUNT}"
    printf 'allocated_gpus=%s\n' "${SLURM_GPUS_ON_NODE:-unknown}"
    nvidia-smi --query-gpu=name,compute_cap,memory.total,driver_version \
        --format=csv,noheader | sed 's/^/gpu=/'
    printf 'python=%s\n' "$("${EMBEDDING_ENV}/bin/python" --version 2>&1)"
    printf 'torch=%s\n' "$("${EMBEDDING_ENV}/bin/python" -c 'import torch;print(torch.__version__)')"
    printf 'cuda_available=%s\n' "$("${EMBEDDING_ENV}/bin/python" -c 'import torch;print(torch.cuda.is_available())')"
    printf 'cuda_home=%s\n' "${CUDA_HOME}"
    printf 'evo2=%s\n' "$("${EMBEDDING_ENV}/bin/python" -c 'import evo2,importlib.metadata as m;print(m.version("evo2"))' 2>/dev/null || echo unavailable)"
    # flash-attn is newer than transformer_engine declares support for; record
    # the observed version with every run rather than assuming it is benign.
    printf 'flash_attn=%s\n' "$("${EMBEDDING_ENV}/bin/python" -c 'import importlib.metadata as m;print(m.version("flash_attn"))' 2>/dev/null || echo unavailable)"
    printf 'model=%s\n' "${AMR_MODEL_ID:-unset}"
    printf 'layer=%s\n' "${AMR_MODEL_LAYER:-unset}"
} | tee "${RESULT_ROOT}/provenance.txt"

"${EMBEDDING_ENV}/bin/python" scripts/extract_genome_embeddings.py plan \
    --assemblies-root "${ASSEMBLIES}" \
    --output "${RESULT_ROOT}/embedding-plan.json" \
    --model-id "${AMR_MODEL_ID:?Set AMR_MODEL_ID}" \
    --model-revision "${AMR_MODEL_REVISION:?Set AMR_MODEL_REVISION to a pinned commit}" \
    --layer "${AMR_MODEL_LAYER:-blocks.6.mlp.l3}" \
    --window "${AMR_WINDOW:-8192}" \
    --minimum-window "${AMR_MINIMUM_WINDOW:-1024}"

(
    cd "${RESULT_ROOT}"
    sha256sum provenance.txt embedding-plan.json > artifact.sha256
)
touch "${RESULT_ROOT}/COMPLETE"
