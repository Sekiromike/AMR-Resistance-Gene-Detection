"""Pretraining-overlap audit for Evo 2 (docs/FOUNDATION_MODEL_PROTOCOL.md).

Evo 2 was trained on OpenGenome2: one representative genome per GTDB species
(releases 214.1 and 220), IMG/PR plasmids, IMG/VR phages and metagenomes. It
saw no phenotypes, so overlap cannot leak MIC labels; the audit reports how
far the sequences embedded here were in-distribution for the model.

- Chromosomes: E. coli is represented in GTDB by a single species
  representative, checked against the cohort's genome accessions.
- Plasmid-borne determinants: every unique acquired unit is searched against
  IMG/PR with blastn; each unit is classed near_identical (>= 99% identity
  over >= 95% of the unit), close (>= 90% over >= 80%) or not_found.

Genotype only. Run from the repository root on the blastn tabular output:
    python scripts/audit_pretraining_overlap.py --blast ... --mapping ... --unique-fasta ... --output ...
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

BLAST_FIELDS = ["qseqid", "sseqid", "pident", "length", "qlen", "qcovhsp", "bitscore"]
CLASSES = ("near_identical", "close", "not_found")


def classify(pident: float, coverage: float) -> str:
    if pident >= 99.0 and coverage >= 95.0:
        return "near_identical"
    if pident >= 90.0 and coverage >= 80.0:
        return "close"
    return "not_found"


def best_hits(blast_tsv: Path) -> dict[str, str]:
    best: dict[str, tuple[float, str]] = {}
    with blast_tsv.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = dict(zip(BLAST_FIELDS, line.rstrip("\n").split("\t")))
            label = classify(float(row["pident"]), float(row["qcovhsp"]))
            rank = CLASSES.index(label)
            current = best.get(row["qseqid"])
            if current is None or rank < CLASSES.index(current[1]):
                best[row["qseqid"]] = (float(row["bitscore"]), label)
    return {query: label for query, (_, label) in best.items()}


def summarize(best: dict[str, str], mapping: Path) -> dict[str, Any]:
    unit_classes: Counter[str] = Counter()
    row_classes: Counter[str] = Counter()
    by_gene: dict[str, Counter[str]] = defaultdict(Counter)
    seen_units: set[str] = set()
    with mapping.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row["unit"].startswith("acquired:") or not row["sequence_sha256"]:
                continue
            label = best.get(row["sequence_sha256"], "not_found")
            row_classes[label] += 1
            gene = row["unit"].split(":", 1)[1].split("#")[0]
            by_gene[gene][label] += 1
            if row["sequence_sha256"] not in seen_units:
                seen_units.add(row["sequence_sha256"])
                unit_classes[label] += 1
    common = sorted(by_gene, key=lambda g: -sum(by_gene[g].values()))[:15]
    return {
        "unique_acquired_units": {c: unit_classes[c] for c in CLASSES},
        "acquired_unit_rows": {c: row_classes[c] for c in CLASSES},
        "most_common_genes": {g: {c: by_gene[g][c] for c in CLASSES} for g in common},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--blast", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--imgpr-sha256", required=True)
    parser.add_argument("--representative-in-cohort", choices=("yes", "no"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summary = summarize(best_hits(args.blast), args.mapping)
    payload = {
        "schema_version": "1.0.0",
        "operation": "audit_pretraining_overlap",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "training_data": "OpenGenome2 (GTDB 214.1/220 species representatives, IMG/PR, IMG/VR, metagenomes)",
        "chromosomes": {"escherichia_coli_representative": "GCF_003697165.2",
                        "representative_in_cohort": args.representative_in_cohort},
        "plasmids": {"database": "OpenGenome2 fasta/plasmids_phage/imgpr.fasta.gz",
                     "sha256": args.imgpr_sha256, **summary,
                     "thresholds": {"near_identical": "identity >= 99, coverage >= 95",
                                    "close": "identity >= 90, coverage >= 80"}},
        "scientific_boundary": {"phenotypes_read": False,
                                "interpretation": "overlap cannot leak MIC labels; it measures familiarity"},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary["unique_acquired_units"]), json.dumps(summary["acquired_unit_rows"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
