#!/usr/bin/env bash
#SBATCH --job-name=amr-disjoint-pilot
#SBATCH --partition=cpu-preempt
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=logs/disjointness-pilot-%j.out
#SBATCH --error=logs/disjointness-pilot-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in config/unity.local.env}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENV_PREFIX="${AMR_DISJOINTNESS_ENV:-${WORK_ROOT}/conda/envs/amr-genomic-disjointness}"
REQUESTS="${PROJECT_ROOT}/data/genomics/disjointness/requests/sequence_requests.tsv"
PILOT_ROOT="${PROJECT_ROOT}/results/genomics/disjointness/pilot/${SLURM_JOB_ID}"

cd "${PROJECT_ROOT}"
mkdir -p logs "${PILOT_ROOT}" "${PILOT_ROOT}/tmp" "${PILOT_ROOT}/sra"
export TMPDIR="${PILOT_ROOT}/tmp"
export VDB_CONFIG="${PILOT_ROOT}/vdb-config"

if [[ ! -x "${ENV_PREFIX}/bin/skani" ]]; then
    printf 'Disjointness environment missing: %s\n' "${ENV_PREFIX}" >&2
    exit 2
fi
export PATH="${ENV_PREFIX}/bin:${PATH}"
if [[ ! -s "${REQUESTS}" ]]; then
    printf 'Frozen request table missing: %s\n' "${REQUESTS}" >&2
    exit 2
fi

dev_accession="$(awk -F '\t' '$2=="development" && $5=="registered_assembly" {print $6; exit}' "${REQUESTS}")"
external_row="$(awk -F '\t' '$2=="external" && $5=="raw_reads" {print $1 "\t" $6 "\t" $7; exit}' "${REQUESTS}")"
IFS=$'\t' read -r external_id external_run expected_md5 <<< "${external_row}"
expected_md5="${expected_md5//$'\r'/}"
if [[ -z "${dev_accession}" || -z "${external_run}" ]]; then
    printf 'Pilot accessions could not be selected from the frozen requests\n' >&2
    exit 2
fi

{
    printf 'started_at_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'development_accession=%s\n' "${dev_accession}"
    printf 'external_sequence_id=%s\n' "${external_id}"
    printf 'external_run=%s\n' "${external_run}"
    printf 'expected_sra_md5=%s\n' "${expected_md5}"
    sha256sum "${REQUESTS}" workflow/envs/disjointness.yaml
    sacct -j "${SLURM_JOB_ID}" --format=JobID,State,ExitCode,Elapsed,AllocCPUS,ReqMem \
        -n -P
} > "${PILOT_ROOT}/provenance.txt"

"${ENV_PREFIX}/bin/datasets" download genome accession "${dev_accession}" \
    --include genome --filename "${PILOT_ROOT}/development.zip"
unzip -q "${PILOT_ROOT}/development.zip" -d "${PILOT_ROOT}/development-package"
dev_fasta="$(find "${PILOT_ROOT}/development-package/ncbi_dataset/data/${dev_accession}" \
    -maxdepth 1 -type f -name '*_genomic.fna' -print -quit)"
if [[ -z "${dev_fasta}" ]]; then
    printf 'Development FASTA was not present in the official package\n' >&2
    exit 3
fi

"${ENV_PREFIX}/bin/prefetch" "${external_run}" \
    --output-directory "${PILOT_ROOT}/sra" --max-size u
sra_file="${PILOT_ROOT}/sra/${external_run}/${external_run}.sra"
"${ENV_PREFIX}/bin/vdb-validate" "${sra_file}" \
    > "${PILOT_ROOT}/vdb-validate.txt" 2>&1
observed_md5="$(md5sum "${sra_file}" | awk '{print $1}')"
printf 'observed_sra_md5=%s\n' "${observed_md5}" >> "${PILOT_ROOT}/provenance.txt"
if [[ "${observed_md5}" != "${expected_md5}" ]]; then
    printf 'SRA MD5 mismatch for %s: expected %s observed %s\n' \
        "${external_run}" "${expected_md5}" "${observed_md5}" >&2
    exit 4
fi

mkdir -p "${PILOT_ROOT}/reads"
"${ENV_PREFIX}/bin/fasterq-dump" "${sra_file}" --split-files \
    --threads "${SLURM_CPUS_PER_TASK}" --outdir "${PILOT_ROOT}/reads" \
    --temp "${PILOT_ROOT}/tmp"
r1="${PILOT_ROOT}/reads/${external_run}_1.fastq"
r2="${PILOT_ROOT}/reads/${external_run}_2.fastq"
if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    printf 'Expected paired FASTQ files were not produced\n' >&2
    exit 5
fi

"${ENV_PREFIX}/bin/shovill" --R1 "${r1}" --R2 "${r2}" \
    --outdir "${PILOT_ROOT}/external-assembly" \
    --cpus "${SLURM_CPUS_PER_TASK}" --ram 20 --depth 150 --minlen 200
external_fasta="${PILOT_ROOT}/external-assembly/contigs.fa"
"${ENV_PREFIX}/bin/skani" dist "${dev_fasta}" "${external_fasta}" \
    -o "${PILOT_ROOT}/skani.tsv"

{
    printf '%s  development.fna\n' "$(sha256sum "${dev_fasta}" | awk '{print $1}')"
    printf '%s  %s.sra\n' "$(sha256sum "${sra_file}" | awk '{print $1}')" "${external_run}"
    printf '%s  %s_1.fastq\n' "$(sha256sum "${r1}" | awk '{print $1}')" "${external_run}"
    printf '%s  %s_2.fastq\n' "$(sha256sum "${r2}" | awk '{print $1}')" "${external_run}"
    printf '%s  external-assembly.fna\n' "$(sha256sum "${external_fasta}" | awk '{print $1}')"
    printf '%s  skani.tsv\n' "$(sha256sum "${PILOT_ROOT}/skani.tsv" | awk '{print $1}')"
} > "${PILOT_ROOT}/artifact.sha256"
printf 'completed_at_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    >> "${PILOT_ROOT}/provenance.txt"
touch "${PILOT_ROOT}/PILOT_COMPLETE"
