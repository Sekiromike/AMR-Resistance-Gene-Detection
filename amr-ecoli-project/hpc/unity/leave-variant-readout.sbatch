#!/usr/bin/env bash
# Pre-stated read-out of docs/LEAVE_VARIANT_OUT_PLAN.md over all targets.
#SBATCH --job-name=amr-leave-readout
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/leave-variant-readout-%j.out
#SBATCH --error=logs/leave-variant-readout-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_LEAVE_JOB:?Set AMR_LEAVE_JOB to the leave-variant-out array job}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PY="${WORK_ROOT}/conda/envs/amr-models/bin/python"
OUT="results/leave-variant-out/${AMR_LEAVE_JOB}"

cd "${PROJECT_ROOT}"
n="$(ls "${OUT}"/predictions_*.csv | wc -l)"
[[ "${n}" == "6" ]] || { printf 'Expected 6 target outputs, found %s\n' "${n}" >&2; exit 2; }
"${PY}" scripts/leave_variant_out.py readout --predictions "${OUT}"/predictions_*.csv --output "${OUT}/readout.json"
(cd "${OUT}" && sha256sum predictions_*.csv predictions_*.json readout.json > artifact.sha256)
touch "${OUT}/COMPLETE"
