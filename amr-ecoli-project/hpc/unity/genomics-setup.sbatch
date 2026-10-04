#!/usr/bin/env bash
# One-time genomics setup for analysis of the verified disjointness assemblies:
#   1. build the QUAST and MLST environments (the AMRFinderPlus one exists);
#   2. freeze the pinned AMRFinderPlus database from the probe download of job
#      64689133, verifying its version file against config before use;
#   3. download and freeze the E. coli species reference, verifying it is
#      type material before use.
# Runs on the non-preemptible cpu partition: Conda solves and network
# downloads must not be killed midway.
#SBATCH --job-name=amr-genomics-setup
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/genomics-setup-%j.out
#SBATCH --error=logs/genomics-setup-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_AMRFINDER_PROBE:?Set AMR_AMRFINDER_PROBE to the probe job id that downloaded the pinned release}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
PKG_CACHE="${AMR_CONDA_PKGS:-${WORK_ROOT}/conda/pkgs}"
ENVS="${WORK_ROOT}/conda/envs"
DISJOINTNESS_ENV="${AMR_DISJOINTNESS_ENV:-${ENVS}/amr-genomic-disjointness}"
CONTROLLER_PY="${ENVS}/amr-ecoli-controller/bin/python"
REFERENCE_ACCESSION="${AMR_SPECIES_REFERENCE:-GCF_003697165.2}"

cd "${PROJECT_ROOT}"
mkdir -p logs provenance "${PKG_CACHE}"
if ! type module >/dev/null 2>&1; then
    set +u
    source /etc/profile.d/z00-lmod-profile.sh
    source /etc/profile.d/z03-lmod-vars.sh
    source /etc/profile.d/z06-lmod-modulepath.sh
    set -u
fi
module load conda/latest
export CONDA_PKGS_DIRS="${PKG_CACHE}"

pinned="$("${CONTROLLER_PY}" -c 'import json;print(json.load(open("config/study.json"))["genomics"]["amrfinder_database_version"])')"

# ---- 1. environments ----
for name in qc mlst; do
    spec="workflow/envs/genomics-${name}.yaml"
    prefix="${ENVS}/amr-genomics-${name}"
    stamp="${prefix}/.amr-genomics-${name}-spec.sha256"
    expected="$(sha256sum "${spec}" | awk '{print $1}')"
    if [[ ! -d "${prefix}/conda-meta" ]]; then
        conda env create --prefix "${prefix}" --file "${spec}"
    elif [[ "$(cat "${stamp}" 2>/dev/null || true)" != "${expected}" ]]; then
        conda env update --prefix "${prefix}" --file "${spec}" --prune
    fi
    printf '%s\n' "${expected}" > "${stamp}"
    conda list --prefix "${prefix}" --explicit > "provenance/genomics-${name}-conda-explicit.txt"
done
# Tools run with their own environment first on PATH: mlst uses
# '#!/usr/bin/env perl' and otherwise picks up the system perl.
in_env() { local env="$1"; shift; PATH="${env}/bin:${PATH}" "$@"; }
quast_version="$(in_env "${ENVS}/amr-genomics-qc" quast.py --version 2>&1 | tail -1)"
mlst_version="$(in_env "${ENVS}/amr-genomics-mlst" mlst --version 2>&1 | tail -1)"
in_env "${ENVS}/amr-genomics-mlst" mlst --check
[[ "${quast_version}" == *"v5.3.0"* ]] || { printf 'QUAST unusable: %s\n' "${quast_version}" >&2; exit 2; }
[[ "${mlst_version}" == "mlst 2.35.0" ]] || { printf 'mlst unusable: %s\n' "${mlst_version}" >&2; exit 2; }
{
    printf 'quast\t%s\n' "${quast_version}"
    printf 'mlst\t%s\n' "${mlst_version}"
    printf 'amrfinder\t%s\n' "$("${ENVS}/amr-genomics-amrfinder/bin/amrfinder" --version 2>&1 | tail -1)"
    printf 'skani\t%s\n' "$("${DISJOINTNESS_ENV}/bin/skani" --version 2>&1 | tail -1)"
} > provenance/genomics-tool-versions.tsv

# ---- 2. AMRFinderPlus database, frozen from the verified probe download ----
probe="results/genomics/amrfinder-probe/${AMR_AMRFINDER_PROBE}/database/${pinned}"
frozen="data/reference/genomics/amrfinderplus/${pinned}"
if [[ ! -f "${probe}/version.txt" ]]; then
    printf 'Probe database not found: %s\n' "${probe}" >&2
    exit 2
fi
observed="$(head -n 1 "${probe}/version.txt" | tr -d '\r\n')"
if [[ "${observed}" != "${pinned}" ]]; then
    printf 'Probe release %s does not match pinned %s\n' "${observed}" "${pinned}" >&2
    exit 2
fi
if [[ ! -d "${frozen}" ]]; then
    mkdir -p "$(dirname "${frozen}")"
    cp -a "${probe}" "${frozen}"
fi
(cd "${frozen}" && find . -type f ! -name 'FROZEN.sha256' -print0 | LC_ALL=C sort -z | xargs -0 sha256sum > FROZEN.sha256)
"${ENVS}/amr-genomics-amrfinder/bin/amrfinder" --database "${frozen}" --database_version \
    > provenance/amrfinder-database-frozen.txt 2>&1

# ---- 3. species reference: must be E. coli type material ----
reference_root="data/reference/genomics/species-reference/${REFERENCE_ACCESSION}"
if [[ ! -f "${reference_root}/COMPLETE" ]]; then
    rm -rf "${reference_root}"
    mkdir -p "${reference_root}"
    "${DISJOINTNESS_ENV}/bin/datasets" summary genome accession "${REFERENCE_ACCESSION}" \
        > "${reference_root}/summary.json"
    # NCBI names type strains with their culture-collection designations
    # (e.g. 'Escherichia coli DSM 30083 = ...') under a strain-level taxid, so
    # species is checked through the taxonomy classification, not the name.
    taxid="$("${CONTROLLER_PY}" -c 'import json,sys;print(json.load(open(sys.argv[1]))["reports"][0]["organism"]["tax_id"])' "${reference_root}/summary.json")"
    "${DISJOINTNESS_ENV}/bin/datasets" summary taxonomy taxon "${taxid}" > "${reference_root}/taxonomy.json"
    "${CONTROLLER_PY}" - "${reference_root}/summary.json" "${reference_root}/taxonomy.json" <<'PY'
import json, sys
report = json.load(open(sys.argv[1]))["reports"][0]
taxonomy = json.load(open(sys.argv[2]))["reports"][0]["taxonomy"]
species = taxonomy.get("classification", {}).get("species", {})
type_material = report.get("type_material", {}) or report.get("assembly_info", {}).get("type_material", {})
if species.get("id") != 562:
    raise SystemExit(f"Reference species is {species!r}, not Escherichia coli (taxid 562)")
if not str(type_material.get("type_label", "")).startswith("TYPE_MATERIAL"):
    raise SystemExit(f"Reference assembly is not flagged as type material: {type_material!r}")
print("reference_organism=" + report["organism"]["organism_name"])
print("reference_species_taxid=562")
print("reference_type_material=" + json.dumps(type_material, sort_keys=True))
PY
    "${DISJOINTNESS_ENV}/bin/datasets" download genome accession "${REFERENCE_ACCESSION}" \
        --include genome --filename "${reference_root}/package.zip"
    (cd "${reference_root}" && unzip -q package.zip)
    fna="$(find "${reference_root}/ncbi_dataset" -name '*.fna' | head -1)"
    cp "${fna}" "${reference_root}/reference.fna"
    (cd "${reference_root}" && sha256sum package.zip reference.fna summary.json taxonomy.json > artifact.sha256)
    touch "${reference_root}/COMPLETE"
fi
(cd "${reference_root}" && sha256sum --check --quiet artifact.sha256)

{
    printf 'job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'amrfinder_database=%s\n' "${frozen}"
    printf 'amrfinder_database_frozen_from_probe_job=%s\n' "${AMR_AMRFINDER_PROBE}"
    printf 'species_reference=%s\n' "${REFERENCE_ACCESSION}"
    printf 'species_reference_sha256=%s\n' "$(sha256sum "${reference_root}/reference.fna" | awk '{print $1}')"
} | tee provenance/genomics-setup.txt
