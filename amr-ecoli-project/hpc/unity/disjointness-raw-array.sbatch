#!/usr/bin/env bash
# Submit with: sbatch --array=1-3507%40 --export=ALL hpc/unity/disjointness-raw-array.sbatch
#SBATCH --job-name=amr-disjoint-raw
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/disjointness-raw-%A_%a.out
#SBATCH --error=logs/disjointness-raw-%A_%a.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${SLURM_ARRAY_TASK_ID:?Submit this script as an array}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENV_PREFIX="${AMR_DISJOINTNESS_ENV:-${WORK_ROOT}/conda/envs/amr-genomic-disjointness}"
REQUESTS="${PROJECT_ROOT}/data/genomics/disjointness/requests/sequence_requests.tsv"
ASSEMBLY_ROOT="${PROJECT_ROOT}/results/genomics/disjointness/assemblies"
export PATH="${ENV_PREFIX}/bin:${PATH}"
if [[ ! -x "${ENV_PREFIX}/bin/shovill" || ! -x "${ENV_PREFIX}/bin/skani" ]]; then
    printf 'Verified disjointness environment is unavailable: %s\n' "${ENV_PREFIX}" >&2
    exit 2
fi

cd "${PROJECT_ROOT}"
row="$(awk -F '\t' -v target="${SLURM_ARRAY_TASK_ID}" \
    'NR>1 && $5=="raw_reads" {n++; if(n==target) {print $2 "\t" $3 "\t" $6 "\t" $7; exit}}' \
    "${REQUESTS}")"
IFS=$'\t' read -r cohort isolate run expected_md5 <<< "${row}"
expected_md5="${expected_md5//$'\r'/}"
if [[ -z "${cohort}" || -z "${isolate}" || -z "${run}" ]]; then
    printf 'No raw-read request at array index %s\n' "${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
fi

result="${ASSEMBLY_ROOT}/${cohort}/${isolate}"
if [[ -f "${result}/COMPLETE" ]]; then
    (cd "${result}" && sha256sum --check artifact.sha256)
    exit 0
fi
if [[ -e "${result}" ]]; then
    printf 'Incomplete result path requires review before retry: %s\n' "${result}" >&2
    exit 2
fi
work="${result}.partial-${SLURM_JOB_ID}"
if [[ -e "${work}" ]]; then
    printf 'Task-private partial path already exists: %s\n' "${work}" >&2
    exit 2
fi
mkdir -p "${work}/sra" "${work}/reads" "${work}/tmp"
export TMPDIR="${work}/tmp"
export VDB_CONFIG="${work}/vdb-config"

{
    printf 'started_at_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'slurm_array_task_id=%s\n' "${SLURM_ARRAY_TASK_ID}"
    printf 'cohort=%s\n' "${cohort}"
    printf 'isolate_id=%s\n' "${isolate}"
    printf 'run_accession=%s\n' "${run}"
    printf 'expected_sra_md5=%s\n' "${expected_md5}"
} > "${work}/provenance.txt"

prefetch_attempt_limit=3
prefetch_attempt=1
while true; do
    attempt_dir="${work}/sra/attempt-${prefetch_attempt}"
    mkdir -p "${attempt_dir}"
    if "${ENV_PREFIX}/bin/prefetch" "${run}" \
        --output-directory "${attempt_dir}" --max-size u; then
        sra_file="${attempt_dir}/${run}/${run}.sra"
        break
    else
        prefetch_exit=$?
    fi
    printf 'prefetch_attempt_%s_exit=%s\nprefetch_attempt_%s_failed_at_utc=%s\n' \
        "${prefetch_attempt}" "${prefetch_exit}" "${prefetch_attempt}" \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "${work}/provenance.txt"
    if (( prefetch_attempt >= prefetch_attempt_limit )); then
        printf 'Prefetch failed after %s attempts for %s\n' \
            "${prefetch_attempt_limit}" "${run}" >&2
        exit "${prefetch_exit}"
    fi
    retry_delay=$((prefetch_attempt * 30))
    printf 'Prefetch attempt %s failed for %s; retrying in %s seconds\n' \
        "${prefetch_attempt}" "${run}" "${retry_delay}" >&2
    sleep "${retry_delay}"
    prefetch_attempt=$((prefetch_attempt + 1))
done
printf 'prefetch_attempts=%s\n' "${prefetch_attempt}" >> "${work}/provenance.txt"
"${ENV_PREFIX}/bin/vdb-validate" "${sra_file}" > "${work}/vdb-validate.txt" 2>&1
observed_md5="$(md5sum "${sra_file}" | awk '{print $1}')"
observed_sha256="$(sha256sum "${sra_file}" | awk '{print $1}')"
printf 'observed_sra_md5=%s\nobserved_sra_sha256=%s\n' \
    "${observed_md5}" "${observed_sha256}" >> "${work}/provenance.txt"
if [[ -n "${expected_md5}" && "${observed_md5}" != "${expected_md5}" ]]; then
    printf 'SRA MD5 mismatch for %s\n' "${run}" >&2
    exit 4
fi

"${ENV_PREFIX}/bin/fasterq-dump" "${sra_file}" --split-files \
    --threads "${SLURM_CPUS_PER_TASK}" --outdir "${work}/reads" --temp "${work}/tmp"
r1="${work}/reads/${run}_1.fastq"
r2="${work}/reads/${run}_2.fastq"
if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    printf 'Expected paired FASTQ files were not produced for %s\n' "${run}" >&2
    exit 5
fi
{
    printf '%s  %s_1.fastq\n' "$(sha256sum "${r1}" | awk '{print $1}')" "${run}"
    printf '%s  %s_2.fastq\n' "$(sha256sum "${r2}" | awk '{print $1}')" "${run}"
} > "${work}/read-files.sha256"

"${ENV_PREFIX}/bin/shovill" --R1 "${r1}" --R2 "${r2}" \
    --outdir "${work}/shovill" --cpus "${SLURM_CPUS_PER_TASK}" \
    --ram 14 --depth 150 --minlen 200
mkdir -p "${result}"
cp "${work}/shovill/contigs.fa" "${result}/assembly.fna"
cp "${work}/provenance.txt" "${work}/vdb-validate.txt" \
    "${work}/read-files.sha256" "${result}/"
cp -r "${work}/shovill" "${result}/shovill"
{
    printf 'completed_at_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'assembly_sha256=%s\n' "$(sha256sum "${result}/assembly.fna" | awk '{print $1}')"
} >> "${result}/provenance.txt"
(cd "${result}" && sha256sum assembly.fna provenance.txt vdb-validate.txt \
    read-files.sha256 > artifact.sha256)
touch "${result}/COMPLETE"
resolved_work="$(realpath -m "${work}")"
case "${resolved_work}" in
    "${ASSEMBLY_ROOT}"/*.partial-"${SLURM_JOB_ID}")
        rm -rf -- "${resolved_work}"
        ;;
    *)
        printf 'Refusing cleanup outside the task-private assembly root: %s\n' \
            "${resolved_work}" >&2
        exit 6
        ;;
esac
