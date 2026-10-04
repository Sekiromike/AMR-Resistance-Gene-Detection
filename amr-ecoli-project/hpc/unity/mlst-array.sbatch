#!/usr/bin/env bash
# Re-type every verified assembly with the MLST scheme frozen in
# config/study.json genomics.mlst (Achtman, ecoli_achtman_4).
#
# Genome-analysis array 64927454 ran "--scheme ecoli", which in mlst >= 2.23 is
# the Pasteur scheme. Its outputs are left untouched as a record; the Achtman
# calls go to a separate tree that finalize_genome_analyses.py reads with
# --mlst-root, verifying each genome's artifact.sha256 and the scheme/loci.
#
# Each task handles AMR_BATCH_SIZE rows; genomes with COMPLETE are skipped.
# Submit with e.g.  --array=0-269%60  for 13,493 rows at batch size 50.
#SBATCH --job-name=amr-mlst
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=logs/mlst-%A_%a.out
#SBATCH --error=logs/mlst-%A_%a.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENVS="${WORK_ROOT}/conda/envs"
REQUESTS="${AMR_REQUEST_LIST:-data/genomics/disjointness/requests/sequence_requests.tsv}"
BATCH="${AMR_BATCH_SIZE:-50}"
ASSEMBLIES="results/genomics/disjointness/assemblies"

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

PY="${ENVS}/amr-ecoli-controller/bin/python"
MLST_ENV="${ENVS}/amr-genomics-mlst"
[[ -x "${MLST_ENV}/bin/mlst" ]] || { printf 'mlst unavailable\n' >&2; exit 2; }
read -r SCHEME EXPECTED_VERSION < <("${PY}" -c '
import json
m = json.load(open("config/study.json"))["genomics"]["mlst"]
print(m["scheme"], m["tool_version"])')
# mlst is "#!/usr/bin/env perl": its environment must lead PATH.
export PATH="${MLST_ENV}/bin:${PATH}"
version="$(mlst --version 2>&1 | tail -1)"
[[ "${version}" == "mlst ${EXPECTED_VERSION}" ]] || { printf 'mlst version %s != %s\n' "${version}" "${EXPECTED_VERSION}" >&2; exit 2; }
mlst --list 2>/dev/null | tr ' ' '\n' | grep -qx "${SCHEME}" || { printf 'Scheme %s not installed\n' "${SCHEME}" >&2; exit 2; }
OUT_ROOT="results/genomics/mlst/${SCHEME}"

start=$(( SLURM_ARRAY_TASK_ID * BATCH ))
mapfile -t rows < <("${PY}" - "${REQUESTS}" "${start}" "${BATCH}" <<'PY'
import csv, sys
path, start, batch = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
for row in rows[start:start + batch]:
    print(f"{row['cohort']}\t{row['isolate_id']}")
PY
)

failures=0
for entry in "${rows[@]}"; do
    cohort="${entry%%$'\t'*}"
    isolate="${entry#*$'\t'}"
    source="${ASSEMBLIES}/${cohort}/${isolate}"
    out="${OUT_ROOT}/${cohort}/${isolate}"
    [[ -f "${out}/COMPLETE" ]] && continue
    if [[ -e "${out}" ]]; then
        mv "${out}" "${out}.partial-${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}"
    fi
    if [[ ! -f "${source}/COMPLETE" ]] || ! (cd "${source}" && grep ' assembly.fna$' artifact.sha256 | sha256sum --check --quiet); then
        printf 'Assembly missing or failed checksum: %s\n' "${source}" >&2
        failures=$(( failures + 1 ))
        continue
    fi
    mkdir -p "${out}"
    set +e
    mlst --scheme "${SCHEME}" --full --threads "${SLURM_CPUS_PER_TASK:-2}" "${source}/assembly.fna" \
        > "${out}/mlst.tsv" 2> "${out}/command.stderr"
    printf '%s\n' "$?" > "${out}/exit_code.txt"
    set -e
    {
        printf 'cohort=%s\n' "${cohort}"
        printf 'isolate_id=%s\n' "${isolate}"
        printf 'assembly_sha256=%s\n' "$(sha256sum "${source}/assembly.fna" | awk '{print $1}')"
        printf 'mlst_scheme=%s\n' "${SCHEME}"
        printf 'mlst_version=%s\n' "${version}"
        printf 'job=%s_%s\n' "${SLURM_JOB_ID}" "${SLURM_ARRAY_TASK_ID}"
    } > "${out}/provenance.txt"
    (cd "${out}" && sha256sum command.stderr exit_code.txt mlst.tsv provenance.txt > artifact.sha256)
    touch "${out}/COMPLETE"
done

if (( failures > 0 )); then
    printf '%d genome(s) in this batch could not be typed\n' "${failures}" >&2
    exit 3
fi
