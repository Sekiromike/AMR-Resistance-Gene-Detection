#!/usr/bin/env bash
# Resume a development clustering run whose skani search completed but whose
# in-memory normalization was OOM-killed (job 64633773: 26h search finished,
# normalize-skani killed at 96 GB on a ~106M-row, 44.6 GB table).
#
# The search is NOT rerun. The saved raw table is normalized and clustered by
# the streaming implementation.
#
# Verification gate: before touching the development run, the streaming code is
# rerun on two completed runs produced by the earlier in-memory code, and must
# reproduce their recorded output SHA-256 values byte for byte. Any mismatch
# stops the job before the development result is written.
#SBATCH --job-name=amr-dev-clusters-resume
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/development-clusters-resume-%j.out
#SBATCH --error=logs/development-clusters-resume-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_RESUME_JOB:?Set AMR_RESUME_JOB to the development-clusters job whose search completed}"
: "${AMR_SOURCE_COMPARISON:?Set AMR_SOURCE_COMPARISON to the COMPLETE dev-vs-external comparison job}"
: "${AMR_VERIFY_EXTERNAL_CLUSTERS:?Set AMR_VERIFY_EXTERNAL_CLUSTERS to the COMPLETE external-clusters job}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
CONTROLLER_ENV="${AMR_CONTROLLER_ENV:-${WORK_ROOT}/conda/envs/amr-ecoli-controller}"
PY="${CONTROLLER_ENV}/bin/python"
REQUESTS="data/genomics/disjointness/requests/sequence_requests.tsv"
CONTRACT="data/genomics/disjointness/requests/disjointness_contract.json"
SOURCE_ROOT="results/genomics/disjointness/comparisons/${AMR_SOURCE_COMPARISON}"
EXT_ROOT="results/genomics/disjointness/external-clusters/${AMR_VERIFY_EXTERNAL_CLUSTERS}"
RESUME_ROOT="results/genomics/disjointness/development-clusters/${AMR_RESUME_JOB}"
VERIFY_ROOT="results/genomics/disjointness/streaming-verification/${SLURM_JOB_ID}"

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
if [[ ! -x "${PY}" ]]; then
    printf 'Verified controller environment is unavailable\n' >&2
    exit 2
fi

for complete in "${SOURCE_ROOT}/COMPLETE" "${EXT_ROOT}/COMPLETE"; do
    if [[ ! -f "${complete}" ]]; then
        printf 'Verification source is not COMPLETE: %s\n' "${complete}" >&2
        exit 2
    fi
done
for required in skani-search.raw.tsv provenance.txt source-verification.txt; do
    if [[ ! -s "${RESUME_ROOT}/${required}" ]]; then
        printf 'Resume run is missing %s\n' "${required}" >&2
        exit 2
    fi
done
for forbidden in COMPLETE comparisons.tsv development-clusters.tsv; do
    if [[ -e "${RESUME_ROOT}/${forbidden}" ]]; then
        printf 'Refusing to overwrite existing %s in %s\n' "${forbidden}" "${RESUME_ROOT}" >&2
        exit 2
    fi
done
if [[ -e "${VERIFY_ROOT}" ]]; then
    printf 'Refusing to overwrite verification result: %s\n' "${VERIFY_ROOT}" >&2
    exit 2
fi
mkdir -p "${VERIFY_ROOT}"
export TMPDIR="${resolved_project}/${VERIFY_ROOT}"

manifest_value() {
    "${PY}" -c 'import json,sys; d=json.load(open(sys.argv[1]));
for k in sys.argv[2].split("."): d=d[k]
print(d)' "$1" "$2"
}

require_equal() {
    local label="$1" expected="$2" observed="$3"
    printf '%s expected=%s observed=%s\n' "${label}" "${expected}" "${observed}" \
        | tee -a "${VERIFY_ROOT}/verification.txt"
    if [[ "${expected}" != "${observed}" ]]; then
        printf 'VERIFICATION FAILED: %s\n' "${label}" | tee -a "${VERIFY_ROOT}/verification.txt" >&2
        exit 3
    fi
}

# ---- Gate 1: dev-vs-external normalization reproduces the recorded output ----
"${PY}" scripts/build_genomic_disjointness.py normalize-skani \
    --raw "${SOURCE_ROOT}/skani-search.raw.tsv" \
    --path-map "${SOURCE_ROOT}/assembly-path-map.tsv" \
    --query-cohort external --reference-cohort development \
    --output "${VERIFY_ROOT}/dev-vs-external.comparisons.tsv" \
    --manifest "${VERIFY_ROOT}/dev-vs-external.normalization-manifest.json" > /dev/null
require_equal "gate1_dev_vs_external_comparisons_sha256" \
    "$(manifest_value "${SOURCE_ROOT}/normalization-manifest.json" outputs.normalized_comparisons.sha256)" \
    "$(manifest_value "${VERIFY_ROOT}/dev-vs-external.normalization-manifest.json" outputs.normalized_comparisons.sha256)"

# ---- Gate 2: external-vs-external normalization and clustering reproduce ----
"${PY}" scripts/build_genomic_disjointness.py normalize-skani \
    --raw "${EXT_ROOT}/skani-search.raw.tsv" \
    --path-map "${SOURCE_ROOT}/assembly-path-map.tsv" \
    --query-cohort external --reference-cohort external \
    --output "${VERIFY_ROOT}/external.comparisons.tsv" \
    --manifest "${VERIFY_ROOT}/external.normalization-manifest.json" > /dev/null
require_equal "gate2_external_comparisons_sha256" \
    "$(manifest_value "${EXT_ROOT}/normalization-manifest.json" outputs.normalized_comparisons.sha256)" \
    "$(manifest_value "${VERIFY_ROOT}/external.normalization-manifest.json" outputs.normalized_comparisons.sha256)"

"${PY}" scripts/build_genomic_disjointness.py cluster-development \
    --requests "${REQUESTS}" --contract "${CONTRACT}" \
    --comparisons "${VERIFY_ROOT}/external.comparisons.tsv" \
    --clusters "${VERIFY_ROOT}/external-clusters.tsv" \
    --manifest "${VERIFY_ROOT}/external.cluster-manifest.json" \
    --cohort external > /dev/null
require_equal "gate2_external_clusters_sha256" \
    "$(manifest_value "${EXT_ROOT}/cluster-manifest.json" outputs.clusters.sha256)" \
    "$(manifest_value "${VERIFY_ROOT}/external.cluster-manifest.json" outputs.clusters.sha256)"

printf 'verification=PASSED\n' | tee -a "${VERIFY_ROOT}/verification.txt"
# The verified tables are byte-identical to retained originals; keep manifests
# and the verification record, drop the duplicate bulk tables.
rm -f "${VERIFY_ROOT}/dev-vs-external.comparisons.tsv" "${VERIFY_ROOT}/external.comparisons.tsv"
(cd "${VERIFY_ROOT}" && sha256sum verification.txt ./*.json external-clusters.tsv > artifact.sha256)

# ---- Resume the development run ----
{
    printf 'resumed_by_job=%s\n' "${SLURM_JOB_ID}"
    printf 'resume_reason=in-memory normalize-skani OOM-killed; search output reused, not rerun\n'
    printf 'normalizer_implementation=streaming\n'
    printf 'streaming_verification=%s\n' "${VERIFY_ROOT}"
    printf 'builder_sha256=%s\n' "$(sha256sum scripts/build_genomic_disjointness.py | awk '{print $1}')"
} >> "${RESUME_ROOT}/provenance.txt"

"${PY}" scripts/build_genomic_disjointness.py normalize-skani \
    --raw "${RESUME_ROOT}/skani-search.raw.tsv" \
    --path-map "${SOURCE_ROOT}/assembly-path-map.tsv" \
    --query-cohort development --reference-cohort development \
    --output "${RESUME_ROOT}/comparisons.tsv" \
    --manifest "${RESUME_ROOT}/normalization-manifest.json" > /dev/null

# Completeness: every development genome must appear as a query and as a
# reference (each matches itself), or the saved search is suspect.
expected_n="$(grep -c . "${SOURCE_ROOT}/development-assemblies.txt")"
require_equal "resume_distinct_queries" "${expected_n}" \
    "$(manifest_value "${RESUME_ROOT}/normalization-manifest.json" counts.distinct_queries)"
require_equal "resume_distinct_references" "${expected_n}" \
    "$(manifest_value "${RESUME_ROOT}/normalization-manifest.json" counts.distinct_references)"

"${PY}" scripts/build_genomic_disjointness.py cluster-development \
    --requests "${REQUESTS}" --contract "${CONTRACT}" \
    --comparisons "${RESUME_ROOT}/comparisons.tsv" \
    --clusters "${RESUME_ROOT}/development-clusters.tsv" \
    --manifest "${RESUME_ROOT}/cluster-manifest.json" \
    --cohort development

(
    cd "${RESUME_ROOT}"
    sha256sum provenance.txt source-verification.txt skani-search.raw.tsv comparisons.tsv \
        normalization-manifest.json development-clusters.tsv cluster-manifest.json \
        > artifact.sha256
)
rm -rf "${RESUME_ROOT}/tmp"
touch "${RESUME_ROOT}/COMPLETE"
