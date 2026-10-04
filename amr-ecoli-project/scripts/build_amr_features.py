"""Extract known AMR determinants per genome from frozen AMRFinderPlus output.

Rules (fixed before any model is fitted; docs/EVALUATION_DESIGN_AMENDMENT.md):

- Only AMRFinderPlus elements of Type "AMR" (acquired genes and point
  mutations, core and plus scope). VIRULENCE and STRESS elements are excluded:
  they are not resistance mechanisms and act mainly as lineage markers.
- An element is present when its Method is a complete call (EXACT, ALLELE,
  BLAST, POINT, HMM) or PARTIAL_CONTIG_END (split by assembly). PARTIAL and
  INTERNAL_STOP calls are absent: truncated or disrupted.
- Features are the element symbols (for example blaCTX-M-15, gyrA_S83L).
  Frequency filtering happens inside training folds, never here.

Each amrfinder.tsv is checked against the genome's artifact.sha256 before use.
Genotype only; no phenotype is read.

Run from the repository root:
    python scripts/build_amr_features.py --genome-manifest ... --output ... --manifest ...
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

COMPLETE_METHOD_PREFIXES = ("EXACT", "ALLELE", "BLAST", "POINT", "HMM", "PARTIAL_CONTIG_END")
OUTPUT_FIELDS = ["isolate_id", "cohort", "element_symbol", "subtype", "class", "subclass", "method"]


class FeatureError(ValueError):
    """Raised when AMRFinderPlus output is missing, altered or malformed."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def recorded_hash(directory: Path, relative: str) -> str:
    for line in (directory / "artifact.sha256").read_text(encoding="utf-8").splitlines():
        digest, _, name = line.partition("  ")
        if name.strip().removeprefix("./") == relative:
            return digest.strip().lower()
    raise FeatureError(f"{relative} is not recorded in {directory}/artifact.sha256")


def is_present(method: str) -> bool:
    return method.strip().upper().startswith(COMPLETE_METHOD_PREFIXES)


def parse_amrfinder(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        needed = {"Element symbol", "Type", "Subtype", "Class", "Subclass", "Method"}
        if reader.fieldnames is None or not needed <= set(reader.fieldnames):
            raise FeatureError(f"{path} lacks AMRFinderPlus 4 columns {sorted(needed)}")
        return list(reader)


def extract(genome_manifest: Path, analysis_root: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    with genome_manifest.open("r", encoding="utf-8", newline="") as handle:
        genomes = [row for row in csv.DictReader(handle, delimiter="\t") if row["eligible_for_modeling"] == "true"]
    rows: list[dict[str, str]] = []
    counts: Counter[str] = Counter()
    for genome in genomes:
        directory = analysis_root / genome["cohort"] / genome["isolate_id"]
        report = directory / "amrfinder" / "amrfinder.tsv"
        if not report.is_file():
            raise FeatureError(f"Missing AMRFinderPlus report: {report}")
        if sha256_file(report) != recorded_hash(directory, "amrfinder/amrfinder.tsv"):
            raise FeatureError(f"AMRFinderPlus report changed since analysis: {report}")
        seen: set[str] = set()
        for hit in parse_amrfinder(report):
            if hit["Type"].strip() != "AMR":
                counts["excluded_type_" + hit["Type"].strip().lower()] += 1
                continue
            if not is_present(hit["Method"]):
                counts["excluded_method_" + hit["Method"].strip().upper()] += 1
                continue
            symbol = hit["Element symbol"].strip()
            if symbol in seen:
                counts["duplicate_calls_collapsed"] += 1
                continue
            seen.add(symbol)
            rows.append({
                "isolate_id": genome["isolate_id"], "cohort": genome["cohort"], "element_symbol": symbol,
                "subtype": hit["Subtype"].strip(), "class": hit["Class"].strip(),
                "subclass": hit["Subclass"].strip(), "method": hit["Method"].strip(),
            })
        counts["genomes"] += 1
        counts["genomes_without_amr_elements"] += not seen
    rows.sort(key=lambda row: (row["cohort"], row["isolate_id"], row["element_symbol"]))
    summary = {
        "genomes": counts.pop("genomes", 0),
        "element_calls": len(rows),
        "distinct_elements": len({row["element_symbol"] for row in rows}),
        "calls_by_subtype": dict(Counter(row["subtype"] for row in rows)),
        "other_counts": dict(sorted(counts.items())),
    }
    return rows, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genome-manifest", type=Path, required=True)
    parser.add_argument("--analysis-root", type=Path, default=Path("results/genomics/analysis"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    try:
        rows, summary = extract(args.genome_manifest, args.analysis_root)
    except FeatureError as exc:
        raise SystemExit(f"AMR feature extraction failed: {exc}") from exc
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "schema_version": "1.0.0",
        "operation": "build_amr_features",
        "generated_at_utc": utc_now(),
        "rules": {"types": ["AMR"], "present_methods": list(COMPLETE_METHOD_PREFIXES),
                  "absent_methods": ["PARTIAL", "INTERNAL_STOP"]},
        "inputs": {"genome_manifest": {"path": str(args.genome_manifest), "sha256": sha256_file(args.genome_manifest)}},
        "summary": summary,
        "outputs": {"features": {"path": str(args.output), "sha256": sha256_file(args.output)}},
        "scientific_boundary": {"phenotypes_read": False},
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
