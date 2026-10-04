#!/usr/bin/env bash
# Deduplicate extracted units across genomes for embedding; checks every task
# output exists and every sequence hash.
#SBATCH --job-name=amr-loci-collect
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/loci-collect-%j.out
#SBATCH --error=logs/loci-collect-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_COHORT_JOB:?Set AMR_COHORT_JOB to the COMPLETE cohort-build job}"
: "${AMR_TASKS:?Set AMR_TASKS to the number of loci-array tasks}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PY="${WORK_ROOT}/conda/envs/amr-ecoli-controller/bin/python"
ROOT="results/genomics/loci/${AMR_COHORT_JOB}"

cd "${PROJECT_ROOT}"
tables=()
for (( i = 0; i < AMR_TASKS; i++ )); do
    t="${ROOT}/units/units-$(printf '%04d' "${i}").tsv"
    [[ -s "${t}" ]] || { printf 'Missing task output: %s\n' "${t}" >&2; exit 2; }
    tables+=("${t}")
done
"${PY}" scripts/extract_resistance_loci.py collect --units "${tables[@]}" \
    --unique-fasta "${ROOT}/unique_units.fna" --mapping "${ROOT}/unit_mapping.tsv" \
    --manifest "${ROOT}/loci_manifest.json"
(cd "${ROOT}" && sha256sum unique_units.fna unit_mapping.tsv loci_manifest.json units/*.tsv > artifact.sha256)
touch "${ROOT}/COMPLETE"
