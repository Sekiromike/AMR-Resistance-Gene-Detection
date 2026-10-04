#!/usr/bin/env bash
# Extract resistance-relevant units (panel loci + acquired determinants) from
# every eligible genome of a COMPLETE cohort build; AMR_BATCH_SIZE genomes per
# task. Submit e.g. --array=0-66 for 13,399 genomes at 200 per task.
#SBATCH --job-name=amr-loci
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/loci-%A_%a.out
#SBATCH --error=logs/loci-%A_%a.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_COHORT_JOB:?Set AMR_COHORT_JOB to the COMPLETE cohort-build job}"
: "${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENVS="${WORK_ROOT}/conda/envs"
PY="${ENVS}/amr-ecoli-controller/bin/python"
BLASTN="${ENVS}/amr-genomics-mlst/bin/blastn"
BATCH="${AMR_BATCH_SIZE:-200}"
PANEL="data/reference/genomics/panel-reference/${AMR_PANEL_REFERENCE:-GCF_000005845.2}"
OUT="results/genomics/loci/${AMR_COHORT_JOB}/units"

cd "${PROJECT_ROOT}"
[[ -f "${PANEL}/COMPLETE" ]] || { printf 'Panel reference not built\n' >&2; exit 2; }
[[ -f "results/cohort/${AMR_COHORT_JOB}/COMPLETE" ]] || { printf 'Cohort not COMPLETE\n' >&2; exit 2; }
mkdir -p "${OUT}"
task_out="${OUT}/units-$(printf '%04d' "${SLURM_ARRAY_TASK_ID}").tsv"
[[ -s "${task_out}" ]] && exit 0
PATH="${ENVS}/amr-genomics-mlst/bin:${PATH}" "${PY}" scripts/extract_resistance_loci.py extract \
    --genome-manifest "results/cohort/${AMR_COHORT_JOB}/genome_manifest.tsv" \
    --panel "${PANEL}/panel.fna" --blastn "${BLASTN}" \
    --start $(( SLURM_ARRAY_TASK_ID * BATCH )) --count "${BATCH}" --output "${task_out}.partial"
mv "${task_out}.partial" "${task_out}"
