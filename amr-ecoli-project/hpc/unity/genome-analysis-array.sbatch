#!/usr/bin/env bash
# Per-genome QC (QUAST), typing (MLST), AMR determinants (AMRFinderPlus) and
# species confirmation (skani vs the frozen E. coli type-strain reference), run
# on the SAME checksum-verified assemblies used for genomic disjointness.
#
# Each array task processes AMR_BATCH_SIZE consecutive rows of the request
# table. A genome with a COMPLETE marker is skipped, so resubmitting the array
# only fills gaps. Tool commands match workflow/rules/genomics.smk.
#
# Submit with e.g.  --array=0-1349%60  for 13,493 rows at batch size 10.
# Sized from canary 64926291: 6 genomes in 3m16s, peak RSS 0.5 GB.
#SBATCH --job-name=amr-genome-analysis
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/genome-analysis-%A_%a.out
#SBATCH --error=logs/genome-analysis-%A_%a.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${SLURM_ARRAY_TASK_ID:?Submit this script as a Slurm array}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENVS="${WORK_ROOT}/conda/envs"
REQUESTS="${AMR_REQUEST_LIST:-data/genomics/disjointness/requests/sequence_requests.tsv}"
BATCH="${AMR_BATCH_SIZE:-10}"
ASSEMBLIES="results/genomics/disjointness/assemblies"
OUT_ROOT="results/genomics/analysis"
REFERENCE_ACCESSION="${AMR_SPECIES_REFERENCE:-GCF_003697165.2}"
REFERENCE="data/reference/genomics/species-reference/${REFERENCE_ACCESSION}/reference.fna"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

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
QUAST="${ENVS}/amr-genomics-qc/bin/quast.py"
MLST="${ENVS}/amr-genomics-mlst/bin/mlst"
AMRFINDER="${ENVS}/amr-genomics-amrfinder/bin/amrfinder"
SKANI="${ENVS}/amr-genomic-disjointness/bin/skani"
for tool in "${PY}" "${QUAST}" "${MLST}" "${AMRFINDER}" "${SKANI}"; do
    [[ -x "${tool}" ]] || { printf 'Tool unavailable: %s\n' "${tool}" >&2; exit 2; }
done
pinned="$("${PY}" -c 'import json;print(json.load(open("config/study.json"))["genomics"]["amrfinder_database_version"])')"
DATABASE="data/reference/genomics/amrfinderplus/${pinned}"
# Achtman scheme from config; in mlst >= 2.23 plain "ecoli" is Pasteur.
MLST_SCHEME="$("${PY}" -c 'import json;print(json.load(open("config/study.json"))["genomics"]["mlst"]["scheme"])')"
[[ -f "${DATABASE}/FROZEN.sha256" ]] || { printf 'AMRFinderPlus database not frozen: %s\n' "${DATABASE}" >&2; exit 2; }
[[ -s "${REFERENCE}" ]] || { printf 'Species reference missing: %s\n' "${REFERENCE}" >&2; exit 2; }
reference_sha="$(sha256sum "${REFERENCE}" | awk '{print $1}')"

# Header-aware extraction of this task's rows.
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
    if [[ -f "${out}/COMPLETE" ]]; then
        continue
    fi
    if [[ -e "${out}" ]]; then
        mv "${out}" "${out}.partial-${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}"
    fi
    if [[ ! -f "${source}/COMPLETE" ]] || ! (cd "${source}" && grep ' assembly.fna$' artifact.sha256 | sha256sum --check --quiet); then
        printf 'Assembly missing or failed checksum: %s\n' "${source}" >&2
        failures=$(( failures + 1 ))
        continue
    fi
    fasta="${source}/assembly.fna"
    mkdir -p "${out}"/{quast,mlst,amrfinder,species}

    set +e
    PATH="${ENVS}/amr-genomics-qc/bin:${PATH}" "${QUAST}" --threads "${THREADS}" --min-contig 0 --output-dir "${out}/quast" "${fasta}" \
        > "${out}/quast/command.stdout" 2> "${out}/quast/command.stderr"
    printf '%s\n' "$?" > "${out}/quast/exit_code.txt"

    # mlst is '#!/usr/bin/env perl': its environment must lead PATH.
    PATH="${ENVS}/amr-genomics-mlst/bin:${PATH}" "${MLST}" --scheme "${MLST_SCHEME}" --full "${fasta}" \
        > "${out}/mlst/mlst.tsv" 2> "${out}/mlst/command.stderr"
    printf '%s\n' "$?" > "${out}/mlst/exit_code.txt"

    PATH="${ENVS}/amr-genomics-amrfinder/bin:${PATH}" "${AMRFINDER}" --nucleotide "${fasta}" --organism Escherichia --plus \
        --name "${isolate}" --threads "${THREADS}" --database "${DATABASE}" \
        --output "${out}/amrfinder/amrfinder.tsv" \
        > "${out}/amrfinder/command.stdout" 2> "${out}/amrfinder/command.stderr"
    printf '%s\n' "$?" > "${out}/amrfinder/exit_code.txt"

    "${SKANI}" dist -q "${fasta}" -r "${REFERENCE}" -o "${out}/species/skani.tsv" \
        > "${out}/species/command.stdout" 2> "${out}/species/command.stderr"
    printf '%s\n' "$?" > "${out}/species/exit_code.txt"
    set -e

    {
        printf 'cohort=%s\n' "${cohort}"
        printf 'isolate_id=%s\n' "${isolate}"
        printf 'assembly_sha256=%s\n' "$(sha256sum "${fasta}" | awk '{print $1}')"
        printf 'amrfinder_database=%s\n' "${pinned}"
        printf 'species_reference=%s\n' "${REFERENCE_ACCESSION}"
        printf 'species_reference_sha256=%s\n' "${reference_sha}"
        printf 'job=%s_%s\n' "${SLURM_JOB_ID}" "${SLURM_ARRAY_TASK_ID}"
    } > "${out}/provenance.txt"
    (
        cd "${out}"
        find . -type f ! -name artifact.sha256 -print0 | LC_ALL=C sort -z | xargs -0 sha256sum > artifact.sha256
    )
    touch "${out}/COMPLETE"
done

if (( failures > 0 )); then
    printf '%d genome(s) in this batch could not be analysed\n' "${failures}" >&2
    exit 3
fi
