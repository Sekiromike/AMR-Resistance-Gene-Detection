"""Apply the MLST adequacy criteria of amendment 004 before any split is made.

The criteria were written and adopted before any cohort-level MLST output was
read (docs/GROUPING_VARIABLE_AMENDMENT.md). They are applied once, to the
development genomes that are eligible for modelling (species confirmed and QC
PASS), separately for the primary (human clinical) and sensitivity (all
sources) populations. No phenotype is read. External genomes are never used.

Per population:

1. typing coverage = resolved / eligible. Below the minimum, the whole cohort
   escalates to cgMLST (checks 2 and 3 are then re-run on cgMLST).
2. largest ST as a fraction of resolved isolates. Above the maximum, the
   amendment is reopened.
3. inverse Simpson index over STs. Below the minimum, the amendment is
   reopened.

Across populations the most severe outcome wins: escalate_cgmlst, then
reopen_amendment_004, then keep_st. Anything but keep_st exits non-zero, so a
downstream split cannot silently proceed.

Run from the repository root:
    python scripts/assess_lineage_adequacy.py --genome-manifest ... \\
        --source-attributes ... --output ...
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

PRIMARY = "human_clinical"
SENSITIVITY = "all_sources"
OUTCOME_ORDER = ("keep_st", "reopen_amendment_004", "escalate_cgmlst")


class LineageAdequacyError(ValueError):
    """Raised when inputs are inconsistent with the adequacy contract."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_criteria(study_config: Path) -> dict[str, Any]:
    config = json.loads(study_config.read_text(encoding="utf-8"))
    try:
        criteria = config["splitting"]["lineage_adequacy"]
    except KeyError as exc:
        raise LineageAdequacyError("config has no splitting.lineage_adequacy block") from exc
    for key in ("min_typing_coverage", "max_largest_lineage_fraction", "min_inverse_simpson"):
        if not isinstance(criteria.get(key), (int, float)):
            raise LineageAdequacyError(f"lineage_adequacy.{key} must be numeric")
    return {**criteria, "folds": config["splitting"].get("folds")}


def measure(lineages: list[str], eligible: int, criteria: dict[str, Any]) -> dict[str, Any]:
    """Return metrics and the outcome for one population.

    ``lineages`` holds one entry per eligible isolate; blank means unresolved.
    """
    resolved = [value for value in lineages if value]
    counts = Counter(resolved)
    metrics: dict[str, Any] = {
        "eligible_isolates": eligible,
        "resolved_isolates": len(resolved),
        "unresolved_isolates": eligible - len(resolved),
        "typing_coverage": len(resolved) / eligible if eligible else 0.0,
        "sequence_types": len(counts),
    }
    if resolved:
        # Integer arithmetic: summed float fractions fall just short of 0.5.
        cumulative, covering = 0, 0
        for _, count in counts.most_common():
            cumulative += count
            covering += 1
            if 2 * cumulative >= len(resolved):
                break
        metrics.update(
            largest_lineage=counts.most_common(1)[0][0],
            largest_lineage_fraction=counts.most_common(1)[0][1] / len(resolved),
            inverse_simpson=len(resolved) ** 2 / sum(n * n for n in counts.values()),  # exact n^2 / sum(n_i^2)
            lineages_covering_half=covering,
            singleton_lineages=sum(1 for count in counts.values() if count == 1),
            top_lineages=[{"lineage_group": st, "isolates": n} for st, n in counts.most_common(10)],
        )
    else:
        metrics.update(largest_lineage="", largest_lineage_fraction=1.0, inverse_simpson=0.0,
                       lineages_covering_half=0, singleton_lineages=0, top_lineages=[])

    checks = {
        "typing_coverage": metrics["typing_coverage"] >= criteria["min_typing_coverage"],
        "concentration": metrics["largest_lineage_fraction"] <= criteria["max_largest_lineage_fraction"],
        "effective_diversity": metrics["inverse_simpson"] >= criteria["min_inverse_simpson"],
    }
    if not checks["typing_coverage"]:
        outcome = "escalate_cgmlst"
    elif not (checks["concentration"] and checks["effective_diversity"]):
        outcome = "reopen_amendment_004"
    else:
        outcome = "keep_st"
    metrics["fold_share"] = 1.0 / criteria["folds"] if criteria.get("folds") else None
    metrics["largest_lineage_exceeds_fold_share"] = (
        metrics["fold_share"] is not None and metrics["largest_lineage_fraction"] > metrics["fold_share"]
    )
    return {"metrics": metrics, "checks": checks, "outcome": outcome}


def assess(
    genome_manifest: Path,
    source_attributes: Path,
    criteria: dict[str, Any],
) -> dict[str, Any]:
    genomes = _read_tsv(genome_manifest)
    attributes = _read_tsv(source_attributes)
    for label, rows, needed in (
        ("genome manifest", genomes, {"isolate_id", "cohort", "lineage_group", "eligible_for_modeling"}),
        ("source attributes", attributes, {"isolate_id", "cohort", "intended_use_population"}),
    ):
        missing = needed.difference(rows[0] if rows else {})
        if missing:
            raise LineageAdequacyError(f"{label} lacks columns {sorted(missing)}")

    population = {
        row["isolate_id"]: row["intended_use_population"]
        for row in attributes
        if row["cohort"] == "development"
    }
    development = [row for row in genomes if row["cohort"] == "development"]
    eligible = [row for row in development if row["eligible_for_modeling"] == "true"]
    unknown = sorted(row["isolate_id"] for row in eligible if row["isolate_id"] not in population)
    if unknown:
        raise LineageAdequacyError(
            f"{len(unknown)} eligible development genomes have no source attributes, e.g. {unknown[:3]}"
        )

    subsets = {
        PRIMARY: [row for row in eligible if population[row["isolate_id"]] == PRIMARY],
        SENSITIVITY: eligible,
    }
    results = {
        name: measure([row["lineage_group"].strip() for row in rows], len(rows), criteria)
        for name, rows in subsets.items()
    }
    outcome = max((result["outcome"] for result in results.values()), key=OUTCOME_ORDER.index)
    return {
        "development_genomes": len(development),
        "eligible_development_genomes": len(eligible),
        "populations": results,
        "outcome": outcome,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genome-manifest", type=Path, required=True)
    parser.add_argument("--source-attributes", type=Path, required=True)
    parser.add_argument("--study-config", type=Path, default=Path("config/study.json"))
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        criteria = load_criteria(args.study_config)
        result = assess(args.genome_manifest, args.source_attributes, criteria)
    except LineageAdequacyError as exc:
        raise SystemExit(f"Lineage adequacy assessment failed: {exc}") from exc
    payload = {
        "schema_version": "1.0.0",
        "operation": "assess_lineage_adequacy",
        "generated_at_utc": utc_now(),
        "amendment": "docs/GROUPING_VARIABLE_AMENDMENT.md",
        "criteria": criteria,
        "inputs": {
            "genome_manifest": {"path": str(args.genome_manifest), "sha256": sha256_file(args.genome_manifest)},
            "source_attributes": {"path": str(args.source_attributes), "sha256": sha256_file(args.source_attributes)},
            "study_config": {"path": str(args.study_config), "sha256": sha256_file(args.study_config)},
        },
        **result,
        "scientific_boundary": {"phenotypes_read": False, "external_genomes_used": False},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"outcome": payload["outcome"], "populations": {
        name: {"outcome": r["outcome"], "checks": r["checks"]} for name, r in payload["populations"].items()
    }}, indent=2))
    return 0 if payload["outcome"] == "keep_st" else 4


if __name__ == "__main__":
    raise SystemExit(main())
