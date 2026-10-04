#!/usr/bin/env bash
# Build the chromosomal resistance-locus panel from E. coli K-12 MG1655
# (FOUNDATION_MODEL_PROTOCOL.md revision 0.2): download genome + GFF3, verify
# the organism, extract each panel gene with 300 bp upstream.
#SBATCH --job-name=amr-loci-reference
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/loci-reference-%j.out
#SBATCH --error=logs/loci-reference-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENVS="${WORK_ROOT}/conda/envs"
PY="${ENVS}/amr-ecoli-controller/bin/python"
DATASETS="${ENVS}/amr-genomic-disjointness/bin/datasets"
ACCESSION="${AMR_PANEL_REFERENCE:-GCF_000005845.2}"
ROOT="data/reference/genomics/panel-reference/${ACCESSION}"

cd "${PROJECT_ROOT}"
if [[ ! -f "${ROOT}/COMPLETE" ]]; then
    rm -rf "${ROOT}"; mkdir -p "${ROOT}"
    "${DATASETS}" summary genome accession "${ACCESSION}" > "${ROOT}/summary.json"
    "${PY}" - "${ROOT}/summary.json" <<'PY'
import json, sys
report = json.load(open(sys.argv[1]))["reports"][0]
name = report["organism"]["organism_name"]
strain = report["organism"].get("infraspecific_names", {}).get("strain", "")
if not name.startswith("Escherichia coli") or "MG1655" not in (name + strain):
    raise SystemExit(f"Panel reference is {name!r} strain {strain!r}, not E. coli K-12 MG1655")
print("panel_reference_organism=" + name, strain)
PY
    "${DATASETS}" download genome accession "${ACCESSION}" --include genome,gff3 --filename "${ROOT}/package.zip"
    (cd "${ROOT}" && unzip -q package.zip)
    cp "$(find "${ROOT}/ncbi_dataset" -name '*.fna' | head -1)" "${ROOT}/genome.fna"
    cp "$(find "${ROOT}/ncbi_dataset" -name '*.gff' | head -1)" "${ROOT}/genomic.gff"
    "${PY}" scripts/extract_resistance_loci.py reference --genome "${ROOT}/genome.fna" \
        --gff "${ROOT}/genomic.gff" --output "${ROOT}/panel.fna" | tee "${ROOT}/panel-lengths.json"
    (cd "${ROOT}" && sha256sum package.zip genome.fna genomic.gff panel.fna summary.json > artifact.sha256)
    touch "${ROOT}/COMPLETE"
fi
(cd "${ROOT}" && sha256sum --check --quiet artifact.sha256)
cat "${ROOT}/panel-lengths.json"
