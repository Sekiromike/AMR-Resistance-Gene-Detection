#!/usr/bin/env bash
# docs/LEAVE_VARIANT_OUT_PLAN.md: one array task per target allele; each fits
# 2 learners x 4 feature sets with the target's carriers held out.
# Development data only.
#SBATCH --job-name=amr-leave-variant
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-5
#SBATCH --output=logs/leave-variant-%A_%a.out
#SBATCH --error=logs/leave-variant-%A_%a.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${SLURM_ARRAY_JOB_ID:?Submit this script as a Slurm array}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PY="${WORK_ROOT}/conda/envs/amr-models/bin/python"
COHORT_JOB="${AMR_COHORT_JOB:-64948707}"
BASE="results/models/${AMR_MODELS_JOB:-64960056}"
FAMILY="results/posthoc/family/${AMR_FAMILY_JOB:-65022960}/amr_family_features.tsv"
TARGETS=(blaCTX-M-15 blaCTX-M-27 blaCTX-M-14 blaCTX-M-55 blaCMY-2 blaOXA-1)
TARGET="${TARGETS[${SLURM_ARRAY_TASK_ID}]}"
OUT="results/leave-variant-out/${SLURM_ARRAY_JOB_ID}"

cd "${PROJECT_ROOT}"
for root in "results/cohort/${COHORT_JOB}" "${BASE}" "$(dirname "${FAMILY}")"; do
    (cd "${root}" && sha256sum --check --quiet artifact.sha256)
done
mkdir -p "${OUT}"
"${PY}" scripts/leave_variant_out.py run --target "${TARGET}" \
    --cohort "results/cohort/${COHORT_JOB}/research_cohort.csv" --splits "${BASE}/isolate_splits.csv" \
    --features "${BASE}/amr_features.tsv" --family-features "${FAMILY}" \
    --embeddings "results/genomics/embeddings/loci/${COHORT_JOB}/evo2_20b" \
    --unit-mapping "results/genomics/loci/${COHORT_JOB}/unit_mapping.tsv" \
    --output "${OUT}/predictions_${TARGET}.csv"
