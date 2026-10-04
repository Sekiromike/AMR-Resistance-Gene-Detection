#!/usr/bin/env bash
# Evo 2 embedding of deduplicated resistance units (protocol revision 0.2).
# H200 NVL GPUs are only in gpu-preempt, so tasks stay under 2 hours and are
# resumable: shards are written atomically and finished shards are skipped.
# Array task i handles shards with index % AMR_SHARD_COUNT == i.
# Canary: AMR_CANARY=1 embeds 20 sequences of one shard to measure VRAM,
# throughput and reverse-complement consistency before any production run.
#SBATCH --job-name=amr-embed-units
#SBATCH --partition=gpu-preempt
#SBATCH --gpus=1
#SBATCH --constraint=h200_nvl
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/embed-units-%A_%a.out
#SBATCH --error=logs/embed-units-%A_%a.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_COHORT_JOB:?Set AMR_COHORT_JOB to the cohort whose loci were collected}"
: "${AMR_SHARD_COUNT:?Set AMR_SHARD_COUNT to the number of array tasks}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
EVO2_ENV="${AMR_EVO2_ENV:-${WORK_ROOT}/.conda/envs/evo2_20b_py311}"
MODEL="${AMR_EVO2_MODEL:-evo2_20b}"
LOCI="results/genomics/loci/${AMR_COHORT_JOB}"
OUT="results/genomics/embeddings/loci/${AMR_COHORT_JOB}/${MODEL}"
[[ "${AMR_CANARY:-0}" == "1" ]] && OUT="${OUT}-canary"

resolved_work="$(realpath -m "${WORK_ROOT}")"
resolved_project="$(realpath -m "${PROJECT_ROOT}")"
case "${resolved_project}" in
    "${resolved_work}"/*) ;;
    *) printf 'Project root is outside AMR_WORK_ROOT: %s\n' "${resolved_project}" >&2; exit 2 ;;
esac
cd "${resolved_project}"
[[ -f "${LOCI}/COMPLETE" ]] || { printf 'Loci not collected: %s\n' "${LOCI}" >&2; exit 2; }
(cd "${LOCI}" && sha256sum --check --quiet artifact.sha256)
mkdir -p "${OUT}" tmp
# Read the existing Evo 2 weights offline; nothing is written to that cache.
export HF_HOME="${AMR_HF_HOME:-${WORK_ROOT}/evo2_project/cache/hf}"
export HF_HUB_OFFLINE=1
export TMPDIR="${resolved_project}/tmp"
# Batch shells do not load Lmod; source it as the other Unity jobs do.
if ! type module >/dev/null 2>&1; then
    set +u
    source /etc/profile.d/z00-lmod-profile.sh
    source /etc/profile.d/z03-lmod-vars.sh
    source /etc/profile.d/z06-lmod-modulepath.sh
    set -u
fi
module purge
module load conda/latest
module load cuda/12.8
set +u; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "${EVO2_ENV}"; set -u
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader

args=(--unique-fasta "${LOCI}/unique_units.fna" --output-root "${OUT}" --model "${MODEL}"
      --shard-index "${SLURM_ARRAY_TASK_ID:-0}" --shard-count "${AMR_SHARD_COUNT}")
if [[ "${AMR_CANARY:-0}" == "1" ]]; then
    args+=(--max-shards 1 --canary-sequences 20)
fi
python -u scripts/embed_resistance_units.py "${args[@]}"
