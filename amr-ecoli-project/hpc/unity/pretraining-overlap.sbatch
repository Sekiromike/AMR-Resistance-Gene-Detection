#!/usr/bin/env bash
# Evo 2 pretraining-overlap audit: search every unique acquired resistance unit
# against the OpenGenome2 IMG/PR plasmid set Evo 2 was trained on.
#SBATCH --job-name=amr-overlap-audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/pretraining-overlap-%j.out
#SBATCH --error=logs/pretraining-overlap-%j.err

set -euo pipefail
: "${AMR_WORK_ROOT:?Set AMR_WORK_ROOT in the submission environment}"
: "${AMR_COHORT_JOB:?Set AMR_COHORT_JOB to the cohort whose loci were collected}"
: "${SLURM_JOB_ID:?Submit this script through Slurm}"
WORK_ROOT="${AMR_WORK_ROOT%/}"
PROJECT_ROOT="${AMR_PROJECT_ROOT:-${WORK_ROOT}/amr-ecoli-project}"
ENVS="${WORK_ROOT}/conda/envs"
PY="${ENVS}/amr-ecoli-controller/bin/python"
BLAST_BIN="${ENVS}/amr-genomics-mlst/bin"
URL="https://huggingface.co/datasets/arcinstitute/opengenome2/resolve/main/fasta/plasmids_phage/imgpr.fasta.gz"
REF="data/reference/opengenome2/imgpr"
LOCI="results/genomics/loci/${AMR_COHORT_JOB}"
OUT="results/genomics/pretraining-overlap/${SLURM_JOB_ID}"

cd "${PROJECT_ROOT}"
[[ -f "${LOCI}/COMPLETE" ]] || { printf 'Loci not collected\n' >&2; exit 2; }
(cd "${LOCI}" && sha256sum --check --quiet artifact.sha256)
[[ -e "${OUT}" ]] && { printf 'Refusing to overwrite: %s\n' "${OUT}" >&2; exit 2; }
mkdir -p "${OUT}" "${REF}"
if [[ ! -f "${REF}/COMPLETE" ]]; then
    curl -sSL --retry 3 -o "${REF}/imgpr.fasta.gz" "${URL}"
    gunzip -c "${REF}/imgpr.fasta.gz" > "${REF}/imgpr.fasta"
    "${BLAST_BIN}/makeblastdb" -in "${REF}/imgpr.fasta" -dbtype nucl -out "${REF}/imgpr" > "${REF}/makeblastdb.log"
    rm "${REF}/imgpr.fasta"
    (cd "${REF}" && sha256sum imgpr.fasta.gz > artifact.sha256)
    touch "${REF}/COMPLETE"
fi
(cd "${REF}" && sha256sum --check --quiet artifact.sha256)
imgpr_sha="$(awk '{print $1}' "${REF}/artifact.sha256")"

"${PY}" - "${LOCI}/unique_units.fna" "${LOCI}/unit_mapping.tsv" "${OUT}/acquired_units.fna" <<'PY'
import csv, sys
fasta, mapping, out = sys.argv[1:4]
acquired = {r["sequence_sha256"] for r in csv.DictReader(open(mapping), delimiter="\t")
            if r["unit"].startswith("acquired:") and r["sequence_sha256"]}
keep = False
with open(fasta) as src, open(out, "w") as dst:
    for line in src:
        if line.startswith(">"):
            keep = line[1:].strip() in acquired
        if keep:
            dst.write(line)
print(f"acquired unique units: {len(acquired)}")
PY
"${BLAST_BIN}/blastn" -task megablast -query "${OUT}/acquired_units.fna" -db "${REF}/imgpr" \
    -perc_identity 90 -evalue 1e-20 -max_target_seqs 5 -max_hsps 1 -num_threads "${SLURM_CPUS_PER_TASK:-8}" \
    -outfmt "6 qseqid sseqid pident length qlen qcovhsp bitscore" > "${OUT}/blast.tsv"
representative="no"
grep -q "003697165" data/genomics/disjointness/requests/sequence_requests.tsv && representative="yes"
"${PY}" scripts/audit_pretraining_overlap.py --blast "${OUT}/blast.tsv" --mapping "${LOCI}/unit_mapping.tsv" \
    --imgpr-sha256 "${imgpr_sha}" --representative-in-cohort "${representative}" --output "${OUT}/overlap.json"
(cd "${OUT}" && sha256sum acquired_units.fna blast.tsv overlap.json > artifact.sha256)
touch "${OUT}/COMPLETE"
