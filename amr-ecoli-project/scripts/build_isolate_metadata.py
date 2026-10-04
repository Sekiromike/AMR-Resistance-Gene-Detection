"""Assemble the frozen isolate-metadata table the cohort builder consumes.

This joins sequence-derived facts only: genomic clusters from the
development-versus-development comparison, evaluation split from the frozen
external membership, and species/QC status from the executed genome manifest.

It never invents a value. Every required field must arrive from a named input,
and an isolate missing any of them is written to the exclusion ledger with a
reason code rather than given a placeholder. Absent AST site, clinical
indication, and surveillance network stay absent; they belong to the blocked
categorical endpoint (docs/ENDPOINT_AMENDMENT.md).

Run from the repository root:
    python scripts/build_isolate_metadata.py --requests ... --output ...
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

OUTPUT_FIELDS = [
    "isolate_id",
    "biosample_accession",
    "assembly_accession",
    "species_method",
    "species_result",
    "genome_qc_status",
    "specimen_source",
    "deduplication_group",
    "lineage_group",
    "genomic_cluster",
    "evaluation_split",
    "intended_use_population",
    "genome_source",
    "genome_source_accession",
    "genome_sha256",
]
POPULATIONS = frozenset({"human_clinical", "non_human_or_environmental", "undetermined"})

EXCLUSION_FIELDS = ["isolate_id", "cohort", "reason", "detail"]
LEDGER_FIELDS = ["isolate_a", "lineage_a", "isolate_b", "lineage_b", "ani_percent", "action"]


class MetadataError(ValueError):
    """Raised when a required input is missing or inconsistent."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_table(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def _require(rows: list[dict[str, str]], columns: tuple[str, ...], label: str) -> None:
    if not rows:
        raise MetadataError(f"{label} is empty")
    missing = sorted(set(columns).difference(rows[0]))
    if missing:
        raise MetadataError(f"{label} is missing required columns: {missing}")


def _blank(value: str | None) -> bool:
    return value is None or not str(value).strip()


class _UnionFind:
    """Size-tracking union-find whose roots are always the smallest member."""

    def __init__(self, items: list[str]) -> None:
        self.parent = {item: item for item in items}
        self.size = {item: 1 for item in items}

    def find(self, item: str) -> str:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != root:
            self.parent[item], item = root, self.parent[item]
        return root

    def union(self, left: str, right: str) -> str:
        a, b = self.find(left), self.find(right)
        if a == b:
            return a
        keep, drop = (a, b) if a < b else (b, a)
        self.parent[drop] = keep
        self.size[keep] += self.size[drop]
        return keep


def build_isolate_metadata(
    requests_path: Path,
    clusters_path: Path | list[Path],
    membership_path: Path,
    genome_manifest_path: Path,
    source_attributes_path: Path,
    duplicate_pairs_path: Path,
    output_path: Path,
    exclusions_path: Path,
    cross_lineage_ledger_path: Path,
    manifest_path: Path,
    *,
    merge_cap_fraction: float = 0.10,
    reopen_cross_lineage_fraction: float = 0.01,
) -> dict[str, Any]:
    """Assemble isolate metadata, applying amendment 005 rules I and G."""
    requests = _read_table(requests_path)
    _require(
        requests,
        ("sequence_id", "cohort", "isolate_id", "biosample_accession", "sequence_accession"),
        "Sequence requests",
    )
    # Development and external clusterings are separate files with distinct
    # cluster-ID prefixes; a sequence may appear in only one of them.
    cluster_paths = list(clusters_path) if isinstance(clusters_path, (list, tuple)) else [clusters_path]
    cluster_by_sequence: dict[str, str] = {}
    for path in cluster_paths:
        clusters = _read_table(path)
        _require(clusters, ("sequence_id", "genomic_cluster"), f"Genomic clusters ({path.name})")
        for row in clusters:
            if row["sequence_id"] in cluster_by_sequence:
                raise MetadataError(f"Sequence {row['sequence_id']} has more than one cluster assignment.")
            cluster_by_sequence[row["sequence_id"]] = row["genomic_cluster"]

    membership = _read_table(membership_path)
    _require(membership, ("sequence_id", "membership"), "External membership")
    membership_by_sequence = {row["sequence_id"]: row for row in membership}

    genomes = _read_table(genome_manifest_path)
    _require(
        genomes,
        ("isolate_id", "assembly_accession", "species_method", "species_result",
         "genome_qc_status", "lineage_group", "eligible_for_modeling", "genome_source",
         "genome_source_accession", "assembly_sha256", "exclusion_reasons"),
        "Genome manifest",
    )
    genome_by_isolate = {row["isolate_id"]: row for row in genomes}

    attributes = _read_table(source_attributes_path)
    _require(attributes, ("isolate_id", "specimen_source", "intended_use_population"), "Source attributes")
    attribute_by_isolate = {row["isolate_id"]: row for row in attributes}
    if len(attribute_by_isolate) != len(attributes):
        raise MetadataError("Source attributes contain duplicate isolate identifiers.")
    unknown = sorted({row["intended_use_population"] for row in attributes} - POPULATIONS)
    if unknown:
        raise MetadataError(f"Unknown intended_use_population values: {unknown}")

    rows: list[dict[str, str]] = []
    excluded: list[dict[str, str]] = []
    request_by_isolate: dict[str, dict[str, str]] = {}
    sequence_to_isolate = {row["sequence_id"]: row["isolate_id"] for row in requests}

    def drop(isolate_id: str, cohort: str, reason: str, detail: str) -> None:
        excluded.append({"isolate_id": isolate_id, "cohort": cohort, "reason": reason, "detail": detail})

    for request in requests:
        sequence_id, isolate_id, cohort = request["sequence_id"], request["isolate_id"], request["cohort"]
        genome = genome_by_isolate.get(isolate_id)
        if genome is None:
            drop(isolate_id, cohort, "genome_analysis_absent", "No executed genome-manifest row exists.")
            continue
        missing = [f for f in ("species_method", "species_result", "genome_qc_status") if _blank(genome.get(f))]
        if missing:
            drop(isolate_id, cohort, "genome_qc_incomplete", f"Genome manifest field is blank: {missing[0]}")
            continue
        if str(genome.get("eligible_for_modeling", "")).strip().lower() != "true":
            drop(isolate_id, cohort, "genome_not_eligible",
                 str(genome.get("exclusion_reasons", "")).strip() or "Genome manifest marks it ineligible.")
            continue
        if cohort == "development":
            evaluation_split = "development"
        elif cohort == "external":
            assignment = membership_by_sequence.get(sequence_id)
            if assignment is None:
                drop(isolate_id, cohort, "external_membership_absent", "No frozen membership row exists.")
                continue
            if assignment["membership"] != "external_locked":
                drop(isolate_id, cohort, "external_excluded_by_membership",
                     f"Membership is {assignment['membership']} under the frozen rule.")
                continue
            evaluation_split = "external"
        else:
            drop(isolate_id, cohort, "unknown_cohort", f"Unrecognised cohort: {cohort}")
            continue
        genomic_cluster = cluster_by_sequence.get(sequence_id)
        if genomic_cluster is None:
            drop(isolate_id, cohort, "genomic_cluster_absent",
                 f"No {cohort} cluster assignment exists; run {cohort} clustering.")
            continue
        attribute = attribute_by_isolate.get(isolate_id, {})
        specimen_source = attribute.get("specimen_source", "")
        population = attribute.get("intended_use_population", "")
        lineage_group = genome.get("lineage_group", "")
        # specimen_source may be absent in source: it stays blank (never
        # inferred) and the isolate is kept; only specimen-type subgroup
        # analyses exclude it (decision owner, 2026-09-27).
        absent = [
            f
            for f, v in (
                ("lineage_group", lineage_group),
                ("intended_use_population", population),
            )
            if _blank(v)
        ]
        if absent:
            drop(isolate_id, cohort, "required_field_absent", f"Absent and not substitutable: {absent[0]}")
            continue
        request_by_isolate[isolate_id] = request
        rows.append(
            {
                "isolate_id": isolate_id,
                "biosample_accession": request["biosample_accession"],
                "assembly_accession": genome["assembly_accession"],
                "species_method": genome["species_method"],
                "species_result": genome["species_result"],
                "genome_qc_status": genome["genome_qc_status"],
                "genome_source": genome["genome_source"],
                "genome_source_accession": genome["genome_source_accession"],
                "genome_sha256": genome["assembly_sha256"],
                "specimen_source": specimen_source,
                "deduplication_group": "",
                "lineage_group": lineage_group,
                "genomic_cluster": genomic_cluster,
                "evaluation_split": evaluation_split,
                "intended_use_population": population,
            }
        )

    isolates = sorted(row["isolate_id"] for row in rows)
    row_by_isolate = {row["isolate_id"]: row for row in rows}

    # ---- Rule I: identifier-level duplicates are one physical isolate ----
    dedup = _UnionFind(isolates)
    owner: dict[str, str] = {}
    identifier_merges = 0
    for isolate_id in isolates:
        request = request_by_isolate[isolate_id]
        values = {
            request["biosample_accession"],
            request["sequence_accession"],
            row_by_isolate[isolate_id]["assembly_accession"],
        }
        for value in sorted(v.strip() for v in values if not _blank(v)):
            if value in owner:
                if dedup.find(owner[value]) != dedup.find(isolate_id):
                    identifier_merges += 1
                dedup.union(owner[value], isolate_id)
            else:
                owner[value] = isolate_id

    # ---- Rule G: genomic duplicates between distinct development isolates ----
    development = [i for i in isolates if row_by_isolate[i]["evaluation_split"] == "development"]
    development_set = set(development)
    lineage_size: dict[str, int] = {}
    for isolate_id in development:
        lineage = row_by_isolate[isolate_id]["lineage_group"]
        lineage_size[lineage] = lineage_size.get(lineage, 0) + 1
    cv = _UnionFind(development)
    first_of_lineage: dict[str, str] = {}
    for isolate_id in development:
        lineage = row_by_isolate[isolate_id]["lineage_group"]
        if lineage in first_of_lineage:
            cv.union(first_of_lineage[lineage], isolate_id)
        else:
            first_of_lineage[lineage] = isolate_id
    # Identifier duplicates (rule I) must share a fold too. A root outside
    # development means a dedup group spans the split; that is rejected below.
    for isolate_id in development:
        root = dedup.find(isolate_id)
        if root in development_set:
            cv.union(root, isolate_id)

    best: dict[tuple[str, str], float] = {}
    duplicate_rows = _read_table(duplicate_pairs_path)
    for pair in duplicate_rows:
        a = sequence_to_isolate.get(pair.get("development_sequence_id", "").strip())
        b = sequence_to_isolate.get(pair.get("external_sequence_id", "").strip())
        if a is None or b is None:
            raise MetadataError("Duplicate-pair table references an unknown sequence ID.")
        if a == b or a not in development_set or b not in development_set:
            continue
        key = (a, b) if a < b else (b, a)
        ani = float(pair["ani_percent"])
        best[key] = max(best.get(key, ani), ani)
    ordered = sorted(best.items(), key=lambda item: (-item[1], item[0][0], item[0][1]))

    cap = int(merge_cap_fraction * len(development))
    removed: set[str] = set()
    ledger: list[dict[str, str]] = []
    within = cross = merged = already = excluded_g = skipped = 0
    for (a, b), ani in ordered:
        if a in removed or b in removed:
            skipped += 1
            continue
        lineage_a = row_by_isolate[a]["lineage_group"]
        lineage_b = row_by_isolate[b]["lineage_group"]
        if lineage_a == lineage_b:
            within += 1
            continue
        cross += 1
        root_a, root_b = cv.find(a), cv.find(b)
        if root_a == root_b:
            already += 1
            action = "already_same_cv_group"
        elif cv.size[root_a] + cv.size[root_b] <= cap:
            cv.union(a, b)
            dedup.union(a, b)
            merged += 1
            action = "merged_lineages"
        else:
            size_a, size_b = lineage_size[lineage_a], lineage_size[lineage_b]
            if size_a != size_b:
                loser = a if size_a < size_b else b
            else:
                loser = max(a, b)
            removed.add(loser)
            excluded_g += 1
            action = f"excluded:{loser}"
            drop(loser, "development", "cross_lineage_genomic_duplicate",
                 f"Duplicate of {b if loser == a else a} in lineage "
                 f"{lineage_b if loser == a else lineage_a}; merge would exceed the {cap}-isolate cap.")
        ledger.append({
            "isolate_a": a, "lineage_a": lineage_a, "isolate_b": b, "lineage_b": lineage_b,
            "ani_percent": f"{ani:.6f}", "action": action,
        })

    total_pairs = within + cross
    if total_pairs and cross / total_pairs > reopen_cross_lineage_fraction:
        raise MetadataError(
            f"{cross} of {total_pairs} genomic duplicate pairs cross lineages "
            f"(>{reopen_cross_lineage_fraction:.0%}); the MLST lineage definition is suspect and "
            "amendment 004 must be reopened rather than patched."
        )

    rows = [row for row in rows if row["isolate_id"] not in removed]
    for row in rows:
        row["deduplication_group"] = f"dedup-{dedup.find(row['isolate_id'])}"
    rows.sort(key=lambda row: row["isolate_id"])
    excluded.sort(key=lambda row: (row["reason"], row["isolate_id"]))

    for field, label in (("genomic_cluster", "Genomic clusters"), ("deduplication_group", "Deduplication groups")):
        spans: dict[str, set[str]] = {}
        for row in rows:
            spans.setdefault(row[field], set()).add(row["evaluation_split"])
        crossing = sorted(k for k, v in spans.items() if len(v) > 1)
        if crossing:
            raise MetadataError(f"{label} cross the development/external boundary: {crossing[:20]}")

    remaining_dev = [r["isolate_id"] for r in rows if r["evaluation_split"] == "development"]
    group_sizes: dict[str, int] = {}
    for isolate_id in remaining_dev:
        root = cv.find(isolate_id)
        group_sizes[root] = group_sizes.get(root, 0) + 1
    largest = max(group_sizes.values()) if group_sizes else 0

    for path, fields, data in (
        (output_path, OUTPUT_FIELDS, rows),
        (exclusions_path, EXCLUSION_FIELDS, excluded),
        (cross_lineage_ledger_path, LEDGER_FIELDS, ledger),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(data)

    payload = {
        "schema_version": "1.1.0",
        "operation": "build_isolate_metadata",
        "amendment": "docs/DEDUPLICATION_POLICY.md",
        "generated_at_utc": utc_now(),
        "inputs": {
            name: {"path": str(path), "sha256": sha256_file(path)}
            for name, path in (
                ("requests", requests_path),
                *((f"clusters_{index}", path) for index, path in enumerate(cluster_paths)),
                ("membership", membership_path), ("genome_manifest", genome_manifest_path),
                ("source_attributes", source_attributes_path), ("duplicate_pairs", duplicate_pairs_path),
            )
        },
        "counts": {
            "requested_isolates": len(requests),
            "metadata_rows": len(rows),
            "excluded": len(excluded),
            "development": len(remaining_dev),
            "development_by_population": dict(Counter(
                row["intended_use_population"] for row in rows if row["evaluation_split"] == "development"
            )),
            "external": sum(1 for row in rows if row["evaluation_split"] == "external"),
            "specimen_source_absent": sum(1 for row in rows if _blank(row["specimen_source"])),
            "genomic_clusters": len({row["genomic_cluster"] for row in rows}),
            "lineages_development": len(lineage_size),
        },
        "rule_i": {"identifier_merges": identifier_merges},
        "rule_g": {
            "genomic_duplicate_pairs": total_pairs,
            "within_lineage": within,
            "cross_lineage": cross,
            "already_same_cv_group": already,
            "merged": merged,
            "excluded": excluded_g,
            "skipped_after_exclusion": skipped,
            "merge_cap_isolates": cap,
            "cross_lineage_fraction": round(cross / total_pairs, 6) if total_pairs else 0.0,
            "largest_cv_group": largest,
            "largest_cv_group_fraction": round(largest / len(remaining_dev), 6) if remaining_dev else 0.0,
        },
        "exclusion_counts": dict(Counter(row["reason"] for row in excluded)),
        "outputs": {
            "metadata": {"path": str(output_path), "sha256": sha256_file(output_path)},
            "exclusions": {"path": str(exclusions_path), "sha256": sha256_file(exclusions_path)},
            "cross_lineage_ledger": {"path": str(cross_lineage_ledger_path),
                                     "sha256": sha256_file(cross_lineage_ledger_path)},
        },
        "scientific_boundary": {
            "phenotypes_read": False,
            "absent_fields_not_substituted": ["site", "clinical_indication", "surveillance_network"],
            "patient_and_outbreak_identifiers": "absent from public sources; not proxied",
        },
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--requests", type=Path, required=True)
    parser.add_argument("--clusters", type=Path, action="append", required=True,
                        help="Cluster table; repeat for development and external.")
    parser.add_argument("--membership", type=Path, required=True)
    parser.add_argument("--genome-manifest", type=Path, required=True)
    parser.add_argument("--source-attributes", type=Path, required=True)
    parser.add_argument("--duplicate-pairs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--exclusions", type=Path, required=True)
    parser.add_argument("--cross-lineage-ledger", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--merge-cap-fraction", type=float, default=0.10)
    parser.add_argument("--reopen-cross-lineage-fraction", type=float, default=0.01)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        payload = build_isolate_metadata(
            args.requests, args.clusters, args.membership, args.genome_manifest,
            args.source_attributes, args.duplicate_pairs, args.output, args.exclusions,
            args.cross_lineage_ledger, args.manifest,
            merge_cap_fraction=args.merge_cap_fraction,
            reopen_cross_lineage_fraction=args.reopen_cross_lineage_fraction,
        )
    except MetadataError as exc:
        raise SystemExit(f"Isolate metadata construction failed: {exc}") from exc
    print(json.dumps({"counts": payload["counts"], "rule_g": payload["rule_g"]}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
