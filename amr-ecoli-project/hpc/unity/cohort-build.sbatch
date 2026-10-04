#!/usr/bin/env bash
# Build the MIC research cohort from the completed genome-analysis array:
#   1. finalize per-genome outputs into the genome manifest;
#   2. apply the amendment 004 MLST adequacy criteria (stops unless ST is kept);
#   3. assemble isolate metadata (amendment 005 rules I and G);
#   4. construct the cohort from the NCBI and JARBS normalized MICs;
#   5. validate it under the mic_regression profile.
# Every source is checksum-verified first and every stage fails closed. The
# external set stays sealed: logs carry counts only, never MIC values,
# distributions or resistance prevalence. Submit with
# --dependency=afterany:<genome-analysis array>; completeness is checked here.
#SBATCH --job-name=amr-cohort-build
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/cohort-build-%j.out
#SBATCH --error=logs/cohort-build-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_DEV_CLUSTERS_JOB:?Set AMR_DEV_CLUSTERS_JOB to the COMPLETE development-clusters job}"
: "${AMR_EXT_CLUSTERS_JOB:?Set AMR_EXT_CLUSTERS_JOB to the COMPLETE external-clusters job}"
: "${AMR_PAIRS_JOB:?Set AMR_PAIRS_JOB to the COMPLETE duplicate-pairs job}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PY="${AMR_CONTROLLER_ENV:-${WORK_ROOT}/conda/envs/amr-ecoli-controller}/bin/python"
REQUESTS="data/genomics/disjointness/requests/sequence_requests.tsv"
DEV_CLUSTERS_ROOT="results/genomics/disjointness/development-clusters/${AMR_DEV_CLUSTERS_JOB}"
EXT_CLUSTERS_ROOT="results/genomics/disjointness/external-clusters/${AMR_EXT_CLUSTERS_JOB}"
PAIRS_ROOT="results/genomics/disjointness/duplicate-pairs/${AMR_PAIRS_JOB}"
MEMBERSHIP="data/genomics/disjointness/membership/external_membership.tsv"
SOURCE_ATTRIBUTES="data/interim/source_attributes.tsv"
NCBI_AST="data/interim/ncbi_ast_normalized.csv"
JARBS_AST="data/interim/jarbs_ast_normalized.csv"
RESULT_ROOT="results/cohort/${SLURM_JOB_ID}"
MLST_SCHEME="$("${PY}" -c 'import json;print(json.load(open("config/study.json"))["genomics"]["mlst"]["scheme"])')"
MLST_ROOT="${AMR_MLST_ROOT:-results/genomics/mlst/${MLST_SCHEME}}"

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

# ---- sources: complete and unchanged ----
for root in "${DEV_CLUSTERS_ROOT}" "${EXT_CLUSTERS_ROOT}" "${PAIRS_ROOT}"; do
    [[ -f "${root}/COMPLETE" ]] || { printf 'Source is not COMPLETE: %s\n' "${root}" >&2; exit 2; }
    (cd "${root}" && sha256sum --check --quiet artifact.sha256)
done
for file in "${MEMBERSHIP}" "${SOURCE_ATTRIBUTES}" "${NCBI_AST}" "${JARBS_AST}"; do
    [[ -s "${file}" ]] || { printf 'Input missing: %s\n' "${file}" >&2; exit 2; }
done
requested="$(tail -n +2 "${REQUESTS}" | wc -l)"
count_complete() {
    "${PY}" - "${REQUESTS}" "$1" <<'PY'
import csv, sys
from pathlib import Path
rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8"), delimiter="\t"))
root = Path(sys.argv[2])
print(sum((root / r["cohort"] / r["isolate_id"] / "COMPLETE").is_file() for r in rows))
PY
}
completed="$(count_complete results/genomics/analysis)"
typed="$(count_complete "${MLST_ROOT}")"
if [[ "${typed}" != "${requested}" ]]; then
    printf 'MLST re-typing incomplete: %s of %s genomes COMPLETE under %s\n' \
        "${typed}" "${requested}" "${MLST_ROOT}" >&2
    exit 2
fi
if [[ "${completed}" != "${requested}" ]]; then
    printf 'Genome analysis incomplete: %s of %s genomes COMPLETE; resubmit the array to fill gaps\n' \
        "${completed}" "${requested}" >&2
    exit 2
fi
if [[ -e "${RESULT_ROOT}" ]]; then
    printf 'Refusing to overwrite cohort result: %s\n' "${RESULT_ROOT}" >&2
    exit 2
fi
mkdir -p "${RESULT_ROOT}"

# ---- 1. genome manifest ----
"${PY}" scripts/finalize_genome_analyses.py \
    --requests "${REQUESTS}" \
    --mlst-root "${MLST_ROOT}" \
    --output "${RESULT_ROOT}/genome_manifest.tsv" \
    --manifest "${RESULT_ROOT}/genome_manifest.json"

# ---- 2. lineage adequacy (amendment 004); exit 4 means escalate or reopen ----
set +e
"${PY}" scripts/assess_lineage_adequacy.py \
    --genome-manifest "${RESULT_ROOT}/genome_manifest.tsv" \
    --source-attributes "${SOURCE_ATTRIBUTES}" \
    --output "${RESULT_ROOT}/lineage_adequacy.json"
adequacy_exit=$?
set -e
if (( adequacy_exit != 0 )); then
    printf 'Lineage adequacy did not keep ST (exit %s); see %s. Stopping before any split.\n' \
        "${adequacy_exit}" "${RESULT_ROOT}/lineage_adequacy.json" >&2
    (cd "${RESULT_ROOT}" && sha256sum genome_manifest.tsv genome_manifest.json lineage_adequacy.json > artifact.sha256)
    exit "${adequacy_exit}"
fi

# ---- 3. isolate metadata (amendment 005) ----
"${PY}" scripts/build_isolate_metadata.py \
    --requests "${REQUESTS}" \
    --clusters "${DEV_CLUSTERS_ROOT}/development-clusters.tsv" \
    --clusters "${EXT_CLUSTERS_ROOT}/external-clusters.tsv" \
    --membership "${MEMBERSHIP}" \
    --genome-manifest "${RESULT_ROOT}/genome_manifest.tsv" \
    --source-attributes "${SOURCE_ATTRIBUTES}" \
    --duplicate-pairs "${PAIRS_ROOT}/development-duplicate-pairs.tsv" \
    --output "${RESULT_ROOT}/isolate_metadata.csv" \
    --exclusions "${RESULT_ROOT}/isolate_metadata_exclusions.csv" \
    --cross-lineage-ledger "${RESULT_ROOT}/cross_lineage_ledger.csv" \
    --manifest "${RESULT_ROOT}/isolate_metadata_manifest.json"

# ---- 4. cohort (both MIC sources) ----
"${PY}" scripts/construct_research_cohort.py \
    --interpreted-ast "${NCBI_AST}" \
    --interpreted-ast "${JARBS_AST}" \
    --isolate-metadata "${RESULT_ROOT}/isolate_metadata.csv" \
    --genome-manifest "${RESULT_ROOT}/genome_manifest.tsv" \
    --output "${RESULT_ROOT}/research_cohort.csv" \
    --exclusions "${RESULT_ROOT}/cohort_exclusions.csv" \
    --manifest "${RESULT_ROOT}/cohort_manifest.json" \
    --endpoint-profile mic_regression

# ---- 5. validation ----
set +e
"${PY}" scripts/validate_research_cohort.py \
    --cohort "${RESULT_ROOT}/research_cohort.csv" \
    --report "${RESULT_ROOT}/cohort_validation.json" \
    --endpoint-profile mic_regression
validation_exit=$?
set -e

{
    printf 'job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'development_clusters_job=%s\n' "${AMR_DEV_CLUSTERS_JOB}"
    printf 'external_clusters_job=%s\n' "${AMR_EXT_CLUSTERS_JOB}"
    printf 'duplicate_pairs_job=%s\n' "${AMR_PAIRS_JOB}"
    printf 'mlst_root=%s\n' "${MLST_ROOT}"
    printf 'validation_exit=%s\n' "${validation_exit}"
    for script in finalize_genome_analyses assess_lineage_adequacy build_isolate_metadata \
                  construct_research_cohort validate_research_cohort; do
        printf '%s_sha256=%s\n' "${script}" "$(sha256sum "scripts/${script}.py" | awk '{print $1}')"
    done
    printf 'study_config_sha256=%s\n' "$(sha256sum config/study.json | awk '{print $1}')"
} > "${RESULT_ROOT}/provenance.txt"
(cd "${RESULT_ROOT}" && find . -maxdepth 1 -type f ! -name artifact.sha256 -printf '%f\n' \
    | LC_ALL=C sort | xargs sha256sum > artifact.sha256)
if (( validation_exit != 0 )); then
    printf 'Cohort validation FAILED; see %s\n' "${RESULT_ROOT}/cohort_validation.json" >&2
    exit "${validation_exit}"
fi
touch "${RESULT_ROOT}/COMPLETE"
