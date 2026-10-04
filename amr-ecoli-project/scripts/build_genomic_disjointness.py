"""Freeze sequence requests and evaluate development/external genome disjointness.

This module never reads external phenotypes.  It consumes only accession-level
sequence manifests, then applies a predeclared ANI/aligned-fraction contract to
normalized skani comparisons.  External membership is therefore unavailable
to model or threshold selection.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np


REQUEST_FIELDS = [
    "sequence_id",
    "cohort",
    "isolate_id",
    "biosample_accession",
    "source_type",
    "sequence_accession",
    "expected_read_md5",
]
COLLISION_FIELDS = [
    "component_id",
    "development_sequence_id",
    "external_sequence_id",
    "ani_percent",
    "aligned_fraction_development",
    "aligned_fraction_external",
    "minimum_reciprocal_aligned_fraction",
    "collision_class",
    "resolution_status",
]
PATH_MAP_FIELDS = ["sequence_id", "cohort", "isolate_id", "assembly_path"]
NORMALIZED_COMPARISON_FIELDS = [
    "development_sequence_id",
    "external_sequence_id",
    "ani_percent",
    "aligned_fraction_development",
    "aligned_fraction_external",
]
SKANI_RAW_FIELDS = [
    "Ref_file",
    "Query_file",
    "ANI",
    "Align_fraction_ref",
    "Align_fraction_query",
]


class DisjointnessError(ValueError):
    """Raised when the comparison contract cannot be verified."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact_entries(path: Path, result_dir: Path) -> tuple[dict[str, str], list[str]]:
    """Parse a GNU sha256sum file without allowing paths outside one result."""
    entries: dict[str, str] = {}
    issues: list[str] = []
    if not path.is_file():
        return entries, ["MISSING_ARTIFACT_MANIFEST"]
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        parts = raw.split(maxsplit=1)
        if len(parts) != 2 or len(parts[0]) != 64 or any(
            char not in "0123456789abcdefABCDEF" for char in parts[0]
        ):
            issues.append(f"INVALID_ARTIFACT_MANIFEST_LINE:{line_number}")
            continue
        relative = parts[1].lstrip("* ")
        candidate = Path(relative)
        if not relative or candidate.is_absolute() or ".." in candidate.parts:
            issues.append(f"UNSAFE_ARTIFACT_PATH:{line_number}")
            continue
        resolved = (result_dir / candidate).resolve()
        try:
            resolved.relative_to(result_dir.resolve())
        except ValueError:
            issues.append(f"UNSAFE_ARTIFACT_PATH:{line_number}")
            continue
        if relative in entries:
            issues.append(f"DUPLICATE_ARTIFACT_ENTRY:{relative}")
            continue
        entries[relative] = parts[0].lower()
    return entries, issues


def _audit_assembly_result(row: Mapping[str, str], assemblies_root: Path) -> dict[str, Any]:
    sequence_id = row["sequence_id"]
    cohort = row["cohort"]
    isolate_id = row["isolate_id"]
    source_type = row["source_type"]
    result_dir = assemblies_root / cohort / isolate_id
    issues: list[str] = []
    if not (result_dir / "COMPLETE").is_file():
        issues.append("MISSING_COMPLETE")
        return {
            "sequence_id": sequence_id,
            "cohort": cohort,
            "isolate_id": isolate_id,
            "source_type": source_type,
            "issues": issues,
        }

    entries, manifest_issues = _artifact_entries(result_dir / "artifact.sha256", result_dir)
    issues.extend(manifest_issues)
    required = {"assembly.fna", "provenance.txt"}
    if source_type == "raw_reads":
        required.update({"read-files.sha256", "vdb-validate.txt"})
    for relative in sorted(required.difference(entries)):
        issues.append(f"REQUIRED_ARTIFACT_NOT_DECLARED:{relative}")
    for relative, expected in sorted(entries.items()):
        artifact = result_dir / relative
        if not artifact.is_file():
            issues.append(f"MISSING_ARTIFACT:{relative}")
        elif sha256_file(artifact) != expected:
            issues.append(f"CHECKSUM_MISMATCH:{relative}")
    provenance: dict[str, str] = {}
    provenance_path = result_dir / "provenance.txt"
    if provenance_path.is_file():
        for raw in provenance_path.read_text(encoding="utf-8").splitlines():
            if "=" in raw:
                key, value = raw.split("=", 1)
                provenance[key] = value
        expected_fields = {
            "isolate_id": isolate_id,
            "assembly_sha256": sha256_file(result_dir / "assembly.fna")
            if (result_dir / "assembly.fna").is_file()
            else "",
        }
        if source_type == "raw_reads":
            expected_fields.update(
                {
                    "cohort": cohort,
                    "run_accession": row["sequence_accession"],
                    "expected_sra_md5": row["expected_read_md5"],
                }
            )
            if row["expected_read_md5"]:
                expected_fields["observed_sra_md5"] = row["expected_read_md5"]
        else:
            expected_fields["assembly_accession"] = row["sequence_accession"]
        for field, expected in expected_fields.items():
            if provenance.get(field) != expected:
                issues.append(f"PROVENANCE_MISMATCH:{field}")
    return {
        "sequence_id": sequence_id,
        "cohort": cohort,
        "isolate_id": isolate_id,
        "source_type": source_type,
        "issues": issues,
    }


def audit_assemblies(
    requests_path: Path,
    assemblies_root: Path,
    manifest_path: Path,
    *,
    workers: int,
) -> dict[str, Any]:
    """Reconcile every frozen request with completion and checksum evidence."""
    if workers < 1:
        raise DisjointnessError("Audit workers must be at least one")
    requests = _read_tsv(requests_path)
    _require_columns(requests, set(REQUEST_FIELDS), "sequence requests")
    expected_keys = {(row["cohort"], row["isolate_id"]) for row in requests}
    if len(expected_keys) != len(requests):
        raise DisjointnessError("Sequence requests contain duplicate cohort/isolate keys")

    with ThreadPoolExecutor(max_workers=workers) as executor:
        rows = list(executor.map(lambda row: _audit_assembly_result(row, assemblies_root), requests))
    failures = [row for row in rows if row["issues"]]
    completed = [row for row in rows if "MISSING_COMPLETE" not in row["issues"]]
    verified = [row for row in rows if not row["issues"]]

    extra_complete: list[str] = []
    for cohort in ("development", "external"):
        cohort_root = assemblies_root / cohort
        if not cohort_root.is_dir():
            continue
        for marker in cohort_root.glob("*/COMPLETE"):
            key = (cohort, marker.parent.name)
            if key not in expected_keys:
                extra_complete.append(f"{cohort}::{marker.parent.name}")
    partials = sorted(
        str(path.relative_to(assemblies_root)).replace("\\", "/")
        for cohort in ("development", "external")
        for path in (assemblies_root / cohort).glob("*.partial-*")
        if path.is_dir()
    )
    complete_and_valid = not failures and not extra_complete
    if complete_and_valid and partials:
        status = "complete_and_verified_with_retained_partials"
    elif complete_and_valid:
        status = "complete_and_verified"
    else:
        status = "incomplete_or_invalid"
    payload = {
        "schema_version": "1.0.0",
        "operation": "audit_genomic_disjointness_assemblies",
        "generated_at_utc": utc_now(),
        "scientific_boundary": "sequence provenance and integrity only; no phenotype or model evidence",
        "inputs": {
            "requests": {"path": str(requests_path), "sha256": sha256_file(requests_path)},
            "assemblies_root": str(assemblies_root),
        },
        "counts": {
            "requests": len(requests),
            "complete_markers": len(completed),
            "checksum_verified": len(verified),
            "failed_or_incomplete": len(failures),
            "extra_complete": len(extra_complete),
            "retained_partial_directories": len(partials),
            "development_verified": sum(
                row["cohort"] == "development" and not row["issues"] for row in rows
            ),
            "external_verified": sum(
                row["cohort"] == "external" and not row["issues"] for row in rows
            ),
        },
        "status": status,
        "failures": failures,
        "extra_complete": sorted(extra_complete),
        "retained_partial_directories": partials,
    }
    _write_json(manifest_path, payload)
    return payload


def prepare_comparison_inputs(
    requests_path: Path,
    audit_manifest_path: Path,
    assemblies_root: Path,
    development_list_path: Path,
    external_list_path: Path,
    path_map_path: Path,
    completed_external_path: Path,
    manifest_path: Path,
    *,
    workers: int = 1,
) -> dict[str, Any]:
    """Build complete comparison lists only after the assembly audit passes."""
    if workers < 1:
        raise DisjointnessError("Comparison preparation workers must be at least one")
    requests = _read_tsv(requests_path)
    _require_columns(requests, set(REQUEST_FIELDS), "sequence requests")
    audit = json.loads(audit_manifest_path.read_text(encoding="utf-8"))
    accepted_status = {
        "complete_and_verified",
        "complete_and_verified_with_retained_partials",
    }
    if audit.get("status") not in accepted_status:
        raise DisjointnessError("Assembly audit is not complete and verified")
    counts = audit.get("counts", {})
    if (
        counts.get("requests") != len(requests)
        or counts.get("checksum_verified") != len(requests)
        or counts.get("failed_or_incomplete") != 0
        or counts.get("extra_complete") != 0
    ):
        raise DisjointnessError("Assembly audit counts do not prove a complete frozen population")
    expected_request_hash = audit.get("inputs", {}).get("requests", {}).get("sha256")
    if expected_request_hash != sha256_file(requests_path):
        raise DisjointnessError("Assembly audit does not authenticate the frozen request table")
    audited_root = audit.get("inputs", {}).get("assemblies_root")
    if not audited_root or Path(audited_root).resolve() != assemblies_root.resolve():
        raise DisjointnessError("Assembly audit does not authenticate the comparison assembly root")

    with ThreadPoolExecutor(max_workers=workers) as executor:
        reauthenticated = list(
            executor.map(lambda row: _audit_assembly_result(row, assemblies_root), requests)
        )
    integrity_failures = [row for row in reauthenticated if row["issues"]]
    if integrity_failures:
        first = integrity_failures[0]
        raise DisjointnessError(
            "Assembly integrity changed after audit: "
            f"{first['sequence_id']} issues={first['issues']} "
            f"(failures={len(integrity_failures)})"
        )

    mapping: list[dict[str, str]] = []
    development_paths: list[str] = []
    external_paths: list[str] = []
    completed_external: list[str] = []
    for row in requests:
        cohort = row["cohort"]
        if cohort not in {"development", "external"}:
            raise DisjointnessError(f"Invalid request cohort: {cohort}")
        result_dir = assemblies_root / cohort / row["isolate_id"]
        assembly = result_dir / "assembly.fna"
        if not (result_dir / "COMPLETE").is_file() or not assembly.is_file():
            raise DisjointnessError(f"Verified assembly is unavailable for {row['sequence_id']}")
        assembly_path = assembly.as_posix()
        mapping.append(
            {
                "sequence_id": row["sequence_id"],
                "cohort": cohort,
                "isolate_id": row["isolate_id"],
                "assembly_path": assembly_path,
            }
        )
        if cohort == "development":
            development_paths.append(assembly_path)
        else:
            external_paths.append(assembly_path)
            completed_external.append(row["sequence_id"])
    if not development_paths or not external_paths:
        raise DisjointnessError("Both development and external assemblies are required")
    if len({row["assembly_path"] for row in mapping}) != len(mapping):
        raise DisjointnessError("Assembly paths are not unique across frozen requests")

    development_list_path.parent.mkdir(parents=True, exist_ok=True)
    development_list_path.write_text(
        "".join(f"{value}\n" for value in development_paths), encoding="utf-8"
    )
    external_list_path.parent.mkdir(parents=True, exist_ok=True)
    external_list_path.write_text(
        "".join(f"{value}\n" for value in external_paths), encoding="utf-8"
    )
    completed_external_path.parent.mkdir(parents=True, exist_ok=True)
    completed_external_path.write_text(
        "".join(f"{value}\n" for value in sorted(completed_external)), encoding="utf-8"
    )
    _write_tsv(path_map_path, mapping, PATH_MAP_FIELDS)
    payload = {
        "schema_version": "1.0.0",
        "operation": "prepare_genomic_disjointness_comparison",
        "generated_at_utc": utc_now(),
        "scientific_boundary": "sequence identity only; external phenotypes remain sealed",
        "inputs": {
            "requests": {"path": str(requests_path), "sha256": sha256_file(requests_path)},
            "assembly_audit": {
                "path": str(audit_manifest_path),
                "sha256": sha256_file(audit_manifest_path),
            },
        },
        "counts": {
            "development_assemblies": len(development_paths),
            "external_assemblies": len(external_paths),
            "total_assemblies": len(mapping),
            "reauthenticated_assemblies": len(reauthenticated),
        },
        "outputs": {
            "development_list": {
                "path": str(development_list_path),
                "sha256": sha256_file(development_list_path),
            },
            "external_list": {
                "path": str(external_list_path),
                "sha256": sha256_file(external_list_path),
            },
            "path_map": {"path": str(path_map_path), "sha256": sha256_file(path_map_path)},
            "completed_external": {
                "path": str(completed_external_path),
                "sha256": sha256_file(completed_external_path),
            },
        },
    }
    _write_json(manifest_path, payload)
    return payload


def normalize_skani_search(
    raw_path: Path,
    path_map_path: Path,
    output_path: Path,
    manifest_path: Path,
    *,
    query_cohort: str = "external",
    reference_cohort: str = "development",
    chunk_rows: int = 2_000_000,
) -> dict[str, Any]:
    """Map skani search paths to frozen sequence IDs and normalize AF to fractions.

    `reference_cohort` and `query_cohort` name the cohorts on each side of the
    search. The output column names keep their cross-cohort spelling in every
    mode; a within-cohort table therefore carries same-cohort IDs in both
    columns, which `cluster_development` verifies before grouping.

    Streams the raw table. A development-versus-development search produces
    ~106 million rows (44.6 GB), which an in-memory list of dicts cannot hold
    in 96 GB. Each row is reduced to an integer pair key plus three float64
    values; because pair keys are built from ranks in the sorted sequence-ID
    list, sorting the integer keys reproduces the original sort on
    (reference_id, query_id) strings exactly, and the output is byte-identical
    to the earlier in-memory implementation.
    """
    if query_cohort not in {"external", "development"}:
        raise DisjointnessError(f"Unsupported query cohort: {query_cohort}")
    if reference_cohort not in {"external", "development"}:
        raise DisjointnessError(f"Unsupported reference cohort: {reference_cohort}")
    if chunk_rows < 1:
        raise DisjointnessError("chunk_rows must be positive")

    mapping_rows = _read_tsv(path_map_path)
    _require_columns(mapping_rows, set(PATH_MAP_FIELDS), "assembly path map")
    path_map = {row["assembly_path"]: row for row in mapping_rows}
    if len(path_map) != len(mapping_rows):
        raise DisjointnessError("Assembly path map contains duplicate paths")
    sequence_ids = sorted({row["sequence_id"] for row in mapping_rows})
    n_ids = len(sequence_ids)
    rank = {sequence_id: index for index, sequence_id in enumerate(sequence_ids)}
    path_index = {
        path: (row["cohort"], rank[row["sequence_id"]]) for path, row in path_map.items()
    }

    key_chunks: list[np.ndarray] = []
    ani_chunks: list[np.ndarray] = []
    afr_chunks: list[np.ndarray] = []
    afq_chunks: list[np.ndarray] = []
    keys: list[int] = []
    anis: list[float] = []
    afrs: list[float] = []
    afqs: list[float] = []
    references_seen: set[int] = set()
    queries_seen: set[int] = set()
    # ~10^4 distinct paths recur across ~10^8 rows; resolve each spelling once.
    resolved: dict[str, tuple[str, int] | None] = {}

    def flush() -> None:
        if keys:
            key_chunks.append(np.asarray(keys, dtype=np.int64))
            ani_chunks.append(np.asarray(anis, dtype=np.float64))
            afr_chunks.append(np.asarray(afrs, dtype=np.float64))
            afq_chunks.append(np.asarray(afqs, dtype=np.float64))
            keys.clear()
            anis.clear()
            afrs.clear()
            afqs.clear()

    def as_float(value: str, field: str) -> float:
        try:
            return float(value)
        except (TypeError, ValueError) as error:
            raise DisjointnessError(f"Invalid {field} in comparison row") from error

    with raw_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None) or []
        missing_raw_fields = sorted(set(SKANI_RAW_FIELDS).difference(header))
        if missing_raw_fields:
            raise DisjointnessError(
                f"skani search output is missing required columns: {missing_raw_fields}"
            )
        i_ref = header.index("Ref_file")
        i_query = header.index("Query_file")
        i_ani = header.index("ANI")
        i_afr = header.index("Align_fraction_ref")
        i_afq = header.index("Align_fraction_query")
        width = max(i_ref, i_query, i_ani, i_afr, i_afq) + 1
        for row in reader:
            if not row:
                continue
            if len(row) < width:
                raise DisjointnessError("skani output row is truncated")
            raw_reference, raw_query = row[i_ref], row[i_query]
            if raw_reference not in resolved:
                resolved[raw_reference] = path_index.get(Path(raw_reference).as_posix())
            if raw_query not in resolved:
                resolved[raw_query] = path_index.get(Path(raw_query).as_posix())
            reference = resolved[raw_reference]
            query = resolved[raw_query]
            if reference is None or query is None:
                raise DisjointnessError("skani output contains an unmapped assembly path")
            if reference[0] != reference_cohort or query[0] != query_cohort:
                raise DisjointnessError(
                    f"skani output is not a {reference_cohort}-reference/"
                    f"{query_cohort}-query pair"
                )
            ani = as_float(row[i_ani], "ANI")
            af_reference = as_float(row[i_afr], "Align_fraction_ref")
            af_query = as_float(row[i_afq], "Align_fraction_query")
            if not (0 <= ani <= 100 and 0 <= af_reference <= 100 and 0 <= af_query <= 100):
                raise DisjointnessError("skani ANI or aligned fraction is outside its valid range")
            references_seen.add(reference[1])
            queries_seen.add(query[1])
            keys.append(reference[1] * n_ids + query[1])
            anis.append(ani)
            afrs.append(af_reference)
            afqs.append(af_query)
            if len(keys) >= chunk_rows:
                flush()
    flush()

    if key_chunks:
        key = np.concatenate(key_chunks)
        ani_all = np.concatenate(ani_chunks)
        afr_all = np.concatenate(afr_chunks)
        afq_all = np.concatenate(afq_chunks)
    else:
        key = np.empty(0, dtype=np.int64)
        ani_all = afr_all = afq_all = np.empty(0, dtype=np.float64)
    del key_chunks, ani_chunks, afr_chunks, afq_chunks

    order = np.argsort(key, kind="stable")
    sorted_key = key[order]
    if sorted_key.size > 1:
        repeats = np.flatnonzero(sorted_key[1:] == sorted_key[:-1])
        if repeats.size:
            duplicate = int(sorted_key[repeats[0]])
            pair = (sequence_ids[duplicate // n_ids], sequence_ids[duplicate % n_ids])
            raise DisjointnessError(f"Duplicate skani comparison pair: {pair}")
    del sorted_key

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(NORMALIZED_COMPARISON_FIELDS)
        for begin in range(0, order.size, chunk_rows):
            index = order[begin : begin + chunk_rows]
            writer.writerows(
                [
                    sequence_ids[k // n_ids],
                    sequence_ids[k % n_ids],
                    f"{a:.6f}",
                    f"{r / 100:.6f}",
                    f"{q / 100:.6f}",
                ]
                for k, a, r, q in zip(
                    key[index].tolist(),
                    ani_all[index].tolist(),
                    afr_all[index].tolist(),
                    afq_all[index].tolist(),
                )
            )

    payload = {
        "schema_version": "1.0.0",
        "operation": "normalize_skani_search",
        "generated_at_utc": utc_now(),
        "query_cohort": query_cohort,
        "reference_cohort": reference_cohort,
        "implementation": "streaming",
        "scientific_boundary": "sequence identity only; external phenotypes remain sealed",
        "inputs": {
            "raw_skani": {"path": str(raw_path), "sha256": sha256_file(raw_path)},
            "path_map": {"path": str(path_map_path), "sha256": sha256_file(path_map_path)},
        },
        "counts": {
            "candidate_comparisons": int(order.size),
            "distinct_references": len(references_seen),
            "distinct_queries": len(queries_seen),
        },
        "outputs": {
            "normalized_comparisons": {
                "path": str(output_path),
                "sha256": sha256_file(output_path),
            }
        },
    }
    _write_json(manifest_path, payload)
    return payload


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, rows: Iterable[Mapping[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _require_columns(rows: list[dict[str, str]], required: set[str], label: str) -> None:
    if not rows:
        raise DisjointnessError(f"{label} is empty")
    missing = sorted(required.difference(rows[0]))
    if missing:
        raise DisjointnessError(f"{label} is missing required columns: {missing}")


def _one_run(raw: str, isolate_id: str) -> str:
    values = sorted({value.strip() for value in raw.replace(";", ",").split(",") if value.strip()})
    if len(values) != 1:
        raise DisjointnessError(
            f"Development isolate {isolate_id} without an assembly has {len(values)} SRA runs"
        )
    return values[0]


def build_sequence_requests(
    development_enrichment: Path,
    external_sequence_manifest: Path,
    *,
    expected_development: int,
    expected_external: int,
) -> list[dict[str, str]]:
    """Build one unambiguous sequence request per isolate without phenotype access."""
    development = _read_csv(development_enrichment)
    external = _read_csv(external_sequence_manifest)
    _require_columns(
        development,
        {
            "target_acc",
            "source_biosample_accession",
            "source_assembly_accession",
            "source_species_taxonomy_id",
            "sra_run_accessions_raw",
        },
        "development enrichment",
    )
    _require_columns(
        external,
        {
            "isolate_id",
            "supplement_species",
            "biosample_accession",
            "primary_short_read_run",
            "primary_short_read_md5",
            "sequence_link_status",
        },
        "external sequence manifest",
    )
    if len(development) != expected_development:
        raise DisjointnessError(
            f"Expected {expected_development} development isolates, found {len(development)}"
        )
    external_ecoli = [row for row in external if row["supplement_species"].strip() == "Escherichia coli"]
    if len(external_ecoli) != expected_external:
        raise DisjointnessError(
            f"Expected {expected_external} external E. coli isolates, found {len(external_ecoli)}"
        )

    requests: list[dict[str, str]] = []
    for row in development:
        isolate = row["target_acc"].strip()
        if row["source_species_taxonomy_id"].strip() != "562":
            raise DisjointnessError(f"Development isolate {isolate} is not taxonomy ID 562")
        assembly = row["source_assembly_accession"].strip()
        source_type = "registered_assembly" if assembly else "raw_reads"
        accession = assembly or _one_run(row["sra_run_accessions_raw"], isolate)
        requests.append(
            {
                "sequence_id": f"development::{isolate}",
                "cohort": "development",
                "isolate_id": isolate,
                "biosample_accession": row["source_biosample_accession"].strip(),
                "source_type": source_type,
                "sequence_accession": accession,
                "expected_read_md5": "",
            }
        )

    for row in external_ecoli:
        isolate = row["isolate_id"].strip()
        if row["sequence_link_status"].strip() != "linked_one_short_read_run":
            raise DisjointnessError(f"External isolate {isolate} lacks one designated short-read run")
        accession = row["primary_short_read_run"].strip()
        checksum = row["primary_short_read_md5"].strip().lower()
        if not accession or len(checksum) != 32 or any(ch not in "0123456789abcdef" for ch in checksum):
            raise DisjointnessError(f"External isolate {isolate} has invalid run identity metadata")
        requests.append(
            {
                "sequence_id": f"external::{isolate}",
                "cohort": "external",
                "isolate_id": isolate,
                "biosample_accession": row["biosample_accession"].strip(),
                "source_type": "raw_reads",
                "sequence_accession": accession,
                "expected_read_md5": checksum,
            }
        )

    for field in ("sequence_id", "biosample_accession", "sequence_accession"):
        seen: dict[str, str] = {}
        for row in requests:
            value = row[field]
            if not value:
                raise DisjointnessError(f"Missing {field} for {row['sequence_id']}")
            if value in seen:
                raise DisjointnessError(f"Duplicate {field} {value}: {seen[value]} and {row['sequence_id']}")
            seen[value] = row["sequence_id"]
    return sorted(requests, key=lambda row: row["sequence_id"])


def write_plan(
    development_enrichment: Path,
    external_sequence_manifest: Path,
    requests_path: Path,
    assembly_accessions_path: Path,
    raw_runs_path: Path,
    contract_path: Path,
    *,
    expected_development: int,
    expected_external: int,
    search_min_ani: float,
    search_min_af: float,
    near_ani: float,
    near_af: float,
    duplicate_ani: float,
    duplicate_af: float,
    skani_version: str,
) -> dict[str, Any]:
    if not (0 <= search_min_ani <= near_ani <= duplicate_ani <= 100):
        raise DisjointnessError("ANI thresholds must satisfy search <= near <= duplicate <= 100")
    if not (0 <= search_min_af <= near_af <= duplicate_af <= 1):
        raise DisjointnessError("Aligned-fraction thresholds must satisfy search <= near <= duplicate <= 1")
    requests = build_sequence_requests(
        development_enrichment,
        external_sequence_manifest,
        expected_development=expected_development,
        expected_external=expected_external,
    )
    _write_tsv(requests_path, requests, REQUEST_FIELDS)
    assemblies = [row["sequence_accession"] for row in requests if row["source_type"] == "registered_assembly"]
    raw_runs = [row["sequence_accession"] for row in requests if row["source_type"] == "raw_reads"]
    assembly_accessions_path.parent.mkdir(parents=True, exist_ok=True)
    assembly_accessions_path.write_text("".join(f"{value}\n" for value in assemblies), encoding="ascii")
    raw_runs_path.parent.mkdir(parents=True, exist_ok=True)
    raw_runs_path.write_text("".join(f"{value}\n" for value in raw_runs), encoding="ascii")
    payload = {
        "schema_version": "1.0.0",
        "operation": "freeze_genomic_disjointness_plan",
        "generated_at_utc": utc_now(),
        "scientific_boundary": {
            "external_phenotypes_read": False,
            "external_membership_available_to_model_selection": False,
            "cross_cohort_collision_policy": "unresolved collisions keep the external cohort unlocked",
            "threshold_interpretation": "conservative leakage-control rule, not a universal biological species definition",
        },
        "inputs": {
            "development_enrichment": {"path": str(development_enrichment), "sha256": sha256_file(development_enrichment)},
            "external_sequence_manifest": {"path": str(external_sequence_manifest), "sha256": sha256_file(external_sequence_manifest)},
        },
        "counts": {
            "development_isolates": sum(row["cohort"] == "development" for row in requests),
            "external_isolates": sum(row["cohort"] == "external" for row in requests),
            "registered_assemblies": len(assemblies),
            "raw_read_runs": len(raw_runs),
        },
        "comparison": {
            "tool": "skani",
            "version": skani_version,
            "learned_ani_regression": True,
            "compression_factor": 125,
            "marker_compression_factor": 1000,
            "search_min_ani_percent": search_min_ani,
            "search_min_reciprocal_aligned_fraction": search_min_af,
            "near_neighbor_min_ani_percent": near_ani,
            "near_neighbor_min_reciprocal_aligned_fraction": near_af,
            "putative_duplicate_min_ani_percent": duplicate_ani,
            "putative_duplicate_min_reciprocal_aligned_fraction": duplicate_af,
            "aligned_fraction_rule": "minimum of development and external aligned fractions",
            "reported_value_rule": "compare the values emitted by skani 0.3.1; ANI and AF are reported to two decimal places",
            "sensitivity_thresholds_reported": [99.5, 99.9, 99.99],
        },
        "outputs": {
            "requests": {"path": str(requests_path), "sha256": sha256_file(requests_path)},
            "assembly_accessions": {"path": str(assembly_accessions_path), "sha256": sha256_file(assembly_accessions_path)},
            "raw_runs": {"path": str(raw_runs_path), "sha256": sha256_file(raw_runs_path)},
        },
    }
    _write_json(contract_path, payload)
    return payload


def _as_float(row: Mapping[str, str], field: str) -> float:
    try:
        return float(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise DisjointnessError(f"Invalid {field} in comparison row") from error


def classify_comparisons(
    comparisons: list[dict[str, str]],
    requests: list[dict[str, str]],
    *,
    near_ani: float,
    near_af: float,
    duplicate_ani: float,
    duplicate_af: float,
) -> list[dict[str, Any]]:
    request_by_id = {row["sequence_id"]: row for row in requests}
    edges: list[dict[str, Any]] = []
    for row in comparisons:
        dev_id = row.get("development_sequence_id", "").strip()
        ext_id = row.get("external_sequence_id", "").strip()
        if request_by_id.get(dev_id, {}).get("cohort") != "development":
            raise DisjointnessError(f"Unknown development sequence ID: {dev_id}")
        if request_by_id.get(ext_id, {}).get("cohort") != "external":
            raise DisjointnessError(f"Unknown external sequence ID: {ext_id}")
        ani = _as_float(row, "ani_percent")
        af_dev = _as_float(row, "aligned_fraction_development")
        af_ext = _as_float(row, "aligned_fraction_external")
        if ani < 0 or ani > 100 or not (0 <= af_dev <= 1 and 0 <= af_ext <= 1):
            raise DisjointnessError("ANI or aligned fraction is outside its valid range")
        minimum_af = min(af_dev, af_ext)
        collision_class = ""
        if ani >= duplicate_ani and minimum_af >= duplicate_af:
            collision_class = "PUTATIVE_GENOMIC_DUPLICATE"
        elif ani >= near_ani and minimum_af >= near_af:
            collision_class = "NEAR_NEIGHBOR"
        if collision_class:
            edges.append(
                {
                    "development_sequence_id": dev_id,
                    "external_sequence_id": ext_id,
                    "ani_percent": f"{ani:.6f}",
                    "aligned_fraction_development": f"{af_dev:.6f}",
                    "aligned_fraction_external": f"{af_ext:.6f}",
                    "minimum_reciprocal_aligned_fraction": f"{minimum_af:.6f}",
                    "collision_class": collision_class,
                    "resolution_status": "UNRESOLVED",
                }
            )

    # Components matter when one external genome collides with multiple development genomes.
    parent: dict[str, str] = {}

    def find(value: str) -> str:
        parent.setdefault(value, value)
        if parent[value] != value:
            parent[value] = find(parent[value])
        return parent[value]

    def union(left: str, right: str) -> None:
        a, b = find(left), find(right)
        if a != b:
            parent[max(a, b)] = min(a, b)

    for edge in edges:
        union(edge["development_sequence_id"], edge["external_sequence_id"])
    roots = sorted({find(value) for value in parent})
    component_ids = {root: f"cross-cohort-{index:06d}" for index, root in enumerate(roots, start=1)}
    for edge in edges:
        edge["component_id"] = component_ids[find(edge["development_sequence_id"])]
    return sorted(edges, key=lambda row: (row["component_id"], row["development_sequence_id"], row["external_sequence_id"]))


def evaluate(
    requests_path: Path,
    contract_path: Path,
    comparisons_path: Path,
    completed_external_path: Path,
    collisions_path: Path,
    manifest_path: Path,
) -> dict[str, Any]:
    requests = _read_tsv(requests_path)
    comparisons = _read_tsv(comparisons_path)
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    expected_external = {row["sequence_id"] for row in requests if row["cohort"] == "external"}
    completed = {line.strip() for line in completed_external_path.read_text(encoding="utf-8").splitlines() if line.strip()}
    if completed != expected_external:
        missing = sorted(expected_external - completed)
        extra = sorted(completed - expected_external)
        raise DisjointnessError(f"External query completion mismatch: missing={len(missing)}, extra={len(extra)}")
    cfg = contract["comparison"]
    collisions = classify_comparisons(
        comparisons,
        requests,
        near_ani=float(cfg["near_neighbor_min_ani_percent"]),
        near_af=float(cfg["near_neighbor_min_reciprocal_aligned_fraction"]),
        duplicate_ani=float(cfg["putative_duplicate_min_ani_percent"]),
        duplicate_af=float(cfg["putative_duplicate_min_reciprocal_aligned_fraction"]),
    )
    _write_tsv(collisions_path, collisions, COLLISION_FIELDS)
    payload = {
        "schema_version": "1.0.0",
        "operation": "evaluate_genomic_disjointness",
        "generated_at_utc": utc_now(),
        "inputs": {
            "requests": {"path": str(requests_path), "sha256": sha256_file(requests_path)},
            "contract": {"path": str(contract_path), "sha256": sha256_file(contract_path)},
            "comparisons": {"path": str(comparisons_path), "sha256": sha256_file(comparisons_path)},
            "completed_external": {"path": str(completed_external_path), "sha256": sha256_file(completed_external_path)},
        },
        "counts": {
            "external_queries_completed": len(completed),
            "candidate_comparisons": len(comparisons),
            "collision_pairs": len(collisions),
            "collision_components": len({row["component_id"] for row in collisions}),
        },
        "external_lock_status": "eligible_for_lock" if not collisions else "blocked_unresolved_genomic_collisions",
        "outputs": {"collisions": {"path": str(collisions_path), "sha256": sha256_file(collisions_path)}},
    }
    _write_json(manifest_path, payload)
    return payload


CLUSTER_FIELDS = [
    "sequence_id",
    "genomic_cluster",
    "cluster_size",
]


def cluster_development(
    requests_path: Path,
    contract_path: Path,
    comparisons_path: Path,
    clusters_path: Path,
    manifest_path: Path,
    *,
    cohort: str = "development",
) -> dict[str, Any]:
    """Assign every development isolate a genomic_cluster for grouped evaluation.

    Clusters are connected components over development-versus-development
    comparisons at the frozen near-neighbour rule. No new threshold is
    introduced: the contract's near_neighbor values are reused so the grouping
    rule and the cross-cohort leakage rule are the same conservative rule.
    Singletons are emitted so every isolate carries a cluster.
    """
    if cohort not in {"development", "external"}:
        raise DisjointnessError(f"Unsupported clustering cohort: {cohort}")
    requests = _read_tsv(requests_path)
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    cfg = contract["comparison"]
    near_ani = float(cfg["near_neighbor_min_ani_percent"])
    near_af = float(cfg["near_neighbor_min_reciprocal_aligned_fraction"])

    development = [row["sequence_id"] for row in requests if row["cohort"] == cohort]
    if not development:
        raise DisjointnessError(f"The request table contains no {cohort} isolates.")
    development_set = set(development)

    parent: dict[str, str] = {sequence_id: sequence_id for sequence_id in development}

    # Iterative find: a recursive one can exceed Python's recursion limit on
    # the long chains that single-linkage produces over 10k genomes. Every
    # component's root is its smallest member either way, so clusters are
    # identical to the recursive version.
    def find(value: str) -> str:
        root = value
        while parent[root] != root:
            root = parent[root]
        while parent[value] != root:
            parent[value], value = root, parent[value]
        return root

    def union(left: str, right: str) -> None:
        a, b = find(left), find(right)
        if a != b:
            parent[max(a, b)] = min(a, b)

    linked = 0
    candidate_comparisons = 0
    comparisons_handle = comparisons_path.open("r", encoding="utf-8", newline="")
    comparisons = csv.DictReader(comparisons_handle, delimiter="\t")
    for row in comparisons:
        candidate_comparisons += 1
        left = row.get("development_sequence_id", "").strip()
        right = row.get("external_sequence_id", "").strip()
        # A development-versus-development search emits both IDs in the
        # development namespace; reject anything else rather than silently
        # clustering across cohorts.
        for sequence_id in (left, right):
            if sequence_id not in development_set:
                raise DisjointnessError(
                    f"{cohort} clustering received an out-of-cohort sequence ID: {sequence_id}"
                )
        if left == right:
            continue
        ani = _as_float(row, "ani_percent")
        af_left = _as_float(row, "aligned_fraction_development")
        af_right = _as_float(row, "aligned_fraction_external")
        if ani < 0 or ani > 100 or not (0 <= af_left <= 1 and 0 <= af_right <= 1):
            raise DisjointnessError("ANI or aligned fraction is outside its valid range")
        if ani >= near_ani and min(af_left, af_right) >= near_af:
            union(left, right)
            linked += 1
    comparisons_handle.close()

    members: dict[str, list[str]] = {}
    for sequence_id in development:
        members.setdefault(find(sequence_id), []).append(sequence_id)
    roots = sorted(members)
    cluster_ids = {root: f"{cohort}-cluster-{index:06d}" for index, root in enumerate(roots, start=1)}

    rows = []
    for root in roots:
        cluster_id = cluster_ids[root]
        size = len(members[root])
        for sequence_id in sorted(members[root]):
            rows.append(
                {
                    "sequence_id": sequence_id,
                    "genomic_cluster": cluster_id,
                    "cluster_size": str(size),
                }
            )
    rows.sort(key=lambda row: (row["genomic_cluster"], row["sequence_id"]))
    _write_tsv(clusters_path, rows, CLUSTER_FIELDS)

    sizes = sorted((len(values) for values in members.values()), reverse=True)
    payload = {
        "schema_version": "1.0.0",
        "operation": "cluster_cohort_genomes",
        "cohort": cohort,
        "generated_at_utc": utc_now(),
        "inputs": {
            "requests": {"path": str(requests_path), "sha256": sha256_file(requests_path)},
            "contract": {"path": str(contract_path), "sha256": sha256_file(contract_path)},
            "comparisons": {"path": str(comparisons_path), "sha256": sha256_file(comparisons_path)},
        },
        "thresholds": {
            "near_neighbor_min_ani_percent": near_ani,
            "near_neighbor_min_reciprocal_aligned_fraction": near_af,
            "source": "frozen disjointness contract; no new threshold introduced",
        },
        "counts": {
            "isolates": len(development),
            "candidate_comparisons": candidate_comparisons,
            "linking_edges": linked,
            "clusters": len(roots),
            "singleton_clusters": sum(1 for size in sizes if size == 1),
            "largest_cluster_size": sizes[0] if sizes else 0,
        },
        "outputs": {"clusters": {"path": str(clusters_path), "sha256": sha256_file(clusters_path)}},
        "scientific_boundary": {
            "phenotypes_read": False,
            "purpose": "grouped cross-validation and internal-holdout construction",
            "threshold_interpretation": "conservative leakage-control rule, not a biological lineage definition",
        },
    }
    _write_json(manifest_path, payload)
    return payload


def extract_duplicate_pairs(
    contract_path: Path,
    comparisons_path: Path,
    pairs_path: Path,
    manifest_path: Path,
) -> dict[str, Any]:
    """Stream a normalized comparison table and keep putative-duplicate pairs.

    Uses the frozen putative-duplicate rule from the contract (amendment 005,
    rule G). Self-comparisons are dropped; both directions of an asymmetric
    pair are kept so the consumer can take the stronger one.
    """
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    cfg = contract["comparison"]
    duplicate_ani = float(cfg["putative_duplicate_min_ani_percent"])
    duplicate_af = float(cfg["putative_duplicate_min_reciprocal_aligned_fraction"])
    scanned = kept = 0
    pairs_path.parent.mkdir(parents=True, exist_ok=True)
    with comparisons_path.open("r", encoding="utf-8", newline="") as source, pairs_path.open(
        "w", encoding="utf-8", newline=""
    ) as sink:
        reader = csv.DictReader(source, delimiter="\t")
        missing = sorted(set(NORMALIZED_COMPARISON_FIELDS).difference(reader.fieldnames or []))
        if missing:
            raise DisjointnessError(f"Comparison table is missing required columns: {missing}")
        writer = csv.writer(sink, delimiter="\t", lineterminator="\n")
        writer.writerow(NORMALIZED_COMPARISON_FIELDS)
        for row in reader:
            scanned += 1
            left = row["development_sequence_id"].strip()
            right = row["external_sequence_id"].strip()
            if left == right:
                continue
            ani = _as_float(row, "ani_percent")
            af = min(
                _as_float(row, "aligned_fraction_development"),
                _as_float(row, "aligned_fraction_external"),
            )
            if ani >= duplicate_ani and af >= duplicate_af:
                writer.writerow([row[field] for field in NORMALIZED_COMPARISON_FIELDS])
                kept += 1
    payload = {
        "schema_version": "1.0.0",
        "operation": "extract_duplicate_pairs",
        "amendment": "docs/DEDUPLICATION_POLICY.md",
        "generated_at_utc": utc_now(),
        "thresholds": {
            "putative_duplicate_min_ani_percent": duplicate_ani,
            "putative_duplicate_min_reciprocal_aligned_fraction": duplicate_af,
            "source": "frozen disjointness contract; no new threshold introduced",
        },
        "inputs": {
            "contract": {"path": str(contract_path), "sha256": sha256_file(contract_path)},
            "comparisons": {"path": str(comparisons_path), "sha256": sha256_file(comparisons_path)},
        },
        "counts": {"comparisons_scanned": scanned, "duplicate_pair_rows": kept},
        "outputs": {"pairs": {"path": str(pairs_path), "sha256": sha256_file(pairs_path)}},
        "scientific_boundary": {"phenotypes_read": False},
    }
    _write_json(manifest_path, payload)
    return payload


MEMBERSHIP_FIELDS = [
    "sequence_id",
    "isolate_id",
    "membership",
    "exclusion_rule",
    "lineage_novel",
    "n_duplicate_partners",
    "n_near_neighbor_partners",
]

EXCLUSION_RULES = {
    # Amendment 003: the duplicate rule means "same clone", which is the
    # leakage that invalidates an external evaluation. The near-neighbour rule
    # spans ordinary within-sequence-type similarity for E. coli and removes
    # lineage membership rather than leakage.
    "duplicate": ("PUTATIVE_GENOMIC_DUPLICATE",),
    "near_neighbor": ("PUTATIVE_GENOMIC_DUPLICATE", "NEAR_NEIGHBOR"),
}


def assign_external_membership(
    requests_path: Path,
    collisions_path: Path,
    membership_path: Path,
    manifest_path: Path,
    *,
    rule: str = "duplicate",
) -> dict[str, Any]:
    """Assign every external isolate to the locked set or the excluded set.

    Selection is sequence-only and phenotype-blind: it reads collision classes
    and nothing else. `lineage_novel` marks isolates with no collision of any
    class, which form the nested strictest stratum under amendment 003.
    """
    if rule not in EXCLUSION_RULES:
        raise DisjointnessError(
            f"Unknown exclusion rule {rule!r}; expected one of {sorted(EXCLUSION_RULES)}"
        )
    excluding = set(EXCLUSION_RULES[rule])
    requests = _read_tsv(requests_path)
    collisions = _read_tsv(collisions_path)

    external = [row for row in requests if row["cohort"] == "external"]
    if not external:
        raise DisjointnessError("The request table contains no external isolates.")

    duplicates: dict[str, int] = {}
    neighbors: dict[str, int] = {}
    for row in collisions:
        sequence_id = row.get("external_sequence_id", "").strip()
        collision_class = row.get("collision_class", "").strip()
        if collision_class == "PUTATIVE_GENOMIC_DUPLICATE":
            duplicates[sequence_id] = duplicates.get(sequence_id, 0) + 1
        elif collision_class == "NEAR_NEIGHBOR":
            neighbors[sequence_id] = neighbors.get(sequence_id, 0) + 1
        else:
            raise DisjointnessError(f"Unknown collision class: {collision_class!r}")

    rows = []
    for row in external:
        sequence_id = row["sequence_id"]
        n_dup = duplicates.get(sequence_id, 0)
        n_near = neighbors.get(sequence_id, 0)
        excluded = (n_dup and "PUTATIVE_GENOMIC_DUPLICATE" in excluding) or (
            n_near and "NEAR_NEIGHBOR" in excluding
        )
        rows.append(
            {
                "sequence_id": sequence_id,
                "isolate_id": row["isolate_id"],
                "membership": "external_excluded" if excluded else "external_locked",
                "exclusion_rule": rule,
                "lineage_novel": "true" if (n_dup == 0 and n_near == 0) else "false",
                "n_duplicate_partners": str(n_dup),
                "n_near_neighbor_partners": str(n_near),
            }
        )
    rows.sort(key=lambda row: row["sequence_id"])
    _write_tsv(membership_path, rows, MEMBERSHIP_FIELDS)

    locked = [row for row in rows if row["membership"] == "external_locked"]
    payload = {
        "schema_version": "1.0.0",
        "operation": "assign_external_membership",
        "generated_at_utc": utc_now(),
        "exclusion_rule": rule,
        "amendment": "docs/EXTERNAL_THRESHOLD_AMENDMENT.md",
        "inputs": {
            "requests": {"path": str(requests_path), "sha256": sha256_file(requests_path)},
            "collisions": {"path": str(collisions_path), "sha256": sha256_file(collisions_path)},
        },
        "counts": {
            "external_isolates": len(rows),
            "external_locked": len(locked),
            "external_excluded": len(rows) - len(locked),
            "lineage_novel": sum(1 for row in rows if row["lineage_novel"] == "true"),
        },
        "outputs": {
            "membership": {"path": str(membership_path), "sha256": sha256_file(membership_path)}
        },
        "scientific_boundary": {
            "phenotypes_read": False,
            "selection_basis": "sequence-only collision class",
        },
    }
    _write_json(manifest_path, payload)
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    plan = sub.add_parser("plan")
    plan.add_argument("--development-enrichment", type=Path, required=True)
    plan.add_argument("--external-sequence-manifest", type=Path, required=True)
    plan.add_argument("--requests", type=Path, required=True)
    plan.add_argument("--assembly-accessions", type=Path, required=True)
    plan.add_argument("--raw-runs", type=Path, required=True)
    plan.add_argument("--contract", type=Path, required=True)
    plan.add_argument("--expected-development", type=int, default=10334)
    plan.add_argument("--expected-external", type=int, default=3159)
    plan.add_argument("--search-min-ani", type=float, default=95.0)
    plan.add_argument("--search-min-af", type=float, default=0.50)
    plan.add_argument("--near-ani", type=float, default=99.9)
    plan.add_argument("--near-af", type=float, default=0.90)
    plan.add_argument("--putative-duplicate-ani", type=float, default=99.99)
    plan.add_argument("--putative-duplicate-af", type=float, default=0.95)
    plan.add_argument("--skani-version", default="0.3.1")
    check = sub.add_parser("evaluate")
    check.add_argument("--requests", type=Path, required=True)
    check.add_argument("--contract", type=Path, required=True)
    check.add_argument("--comparisons", type=Path, required=True)
    check.add_argument("--completed-external", type=Path, required=True)
    check.add_argument("--collisions", type=Path, required=True)
    check.add_argument("--manifest", type=Path, required=True)
    audit = sub.add_parser("audit-assemblies")
    audit.add_argument("--requests", type=Path, required=True)
    audit.add_argument("--assemblies-root", type=Path, required=True)
    audit.add_argument("--manifest", type=Path, required=True)
    audit.add_argument("--workers", type=int, default=4)
    prepare = sub.add_parser("prepare-comparison")
    prepare.add_argument("--requests", type=Path, required=True)
    prepare.add_argument("--audit-manifest", type=Path, required=True)
    prepare.add_argument("--assemblies-root", type=Path, required=True)
    prepare.add_argument("--development-list", type=Path, required=True)
    prepare.add_argument("--external-list", type=Path, required=True)
    prepare.add_argument("--path-map", type=Path, required=True)
    prepare.add_argument("--completed-external", type=Path, required=True)
    prepare.add_argument("--manifest", type=Path, required=True)
    prepare.add_argument("--workers", type=int, default=1)
    normalize = sub.add_parser("normalize-skani")
    normalize.add_argument("--raw", type=Path, required=True)
    normalize.add_argument("--path-map", type=Path, required=True)
    normalize.add_argument("--output", type=Path, required=True)
    normalize.add_argument("--manifest", type=Path, required=True)
    normalize.add_argument(
        "--query-cohort", choices=("external", "development"), default="external"
    )
    normalize.add_argument(
        "--reference-cohort", choices=("external", "development"), default="development"
    )

    duplicates = sub.add_parser("extract-duplicate-pairs")
    duplicates.add_argument("--contract", type=Path, required=True)
    duplicates.add_argument("--comparisons", type=Path, required=True)
    duplicates.add_argument("--pairs", type=Path, required=True)
    duplicates.add_argument("--manifest", type=Path, required=True)

    membership = sub.add_parser("assign-external-membership")
    membership.add_argument("--requests", type=Path, required=True)
    membership.add_argument("--collisions", type=Path, required=True)
    membership.add_argument("--membership", type=Path, required=True)
    membership.add_argument("--manifest", type=Path, required=True)
    membership.add_argument("--rule", choices=sorted(EXCLUSION_RULES), default="duplicate")

    cluster = sub.add_parser("cluster-development")
    cluster.add_argument("--requests", type=Path, required=True)
    cluster.add_argument("--contract", type=Path, required=True)
    cluster.add_argument("--comparisons", type=Path, required=True)
    cluster.add_argument("--clusters", type=Path, required=True)
    cluster.add_argument("--manifest", type=Path, required=True)
    cluster.add_argument("--cohort", choices=("development", "external"), default="development")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "plan":
        payload = write_plan(
            args.development_enrichment,
            args.external_sequence_manifest,
            args.requests,
            args.assembly_accessions,
            args.raw_runs,
            args.contract,
            expected_development=args.expected_development,
            expected_external=args.expected_external,
            search_min_ani=args.search_min_ani,
            search_min_af=args.search_min_af,
            near_ani=args.near_ani,
            near_af=args.near_af,
            duplicate_ani=args.putative_duplicate_ani,
            duplicate_af=args.putative_duplicate_af,
            skani_version=args.skani_version,
        )
    elif args.command == "evaluate":
        payload = evaluate(
            args.requests,
            args.contract,
            args.comparisons,
            args.completed_external,
            args.collisions,
            args.manifest,
        )
    elif args.command == "audit-assemblies":
        payload = audit_assemblies(
            args.requests,
            args.assemblies_root,
            args.manifest,
            workers=args.workers,
        )
    elif args.command == "prepare-comparison":
        payload = prepare_comparison_inputs(
            args.requests,
            args.audit_manifest,
            args.assemblies_root,
            args.development_list,
            args.external_list,
            args.path_map,
            args.completed_external,
            args.manifest,
            workers=args.workers,
        )
    elif args.command == "extract-duplicate-pairs":
        payload = extract_duplicate_pairs(args.contract, args.comparisons, args.pairs, args.manifest)
    elif args.command == "assign-external-membership":
        payload = assign_external_membership(
            args.requests,
            args.collisions,
            args.membership,
            args.manifest,
            rule=args.rule,
        )
    elif args.command == "cluster-development":
        payload = cluster_development(
            args.requests,
            args.contract,
            args.comparisons,
            args.clusters,
            args.manifest,
            cohort=args.cohort,
        )
    else:
        payload = normalize_skani_search(
            args.raw,
            args.path_map,
            args.output,
            args.manifest,
            query_cohort=args.query_cohort,
            reference_cohort=args.reference_cohort,
        )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 1 if args.command == "audit-assemblies" and payload["status"] == "incomplete_or_invalid" else 0


if __name__ == "__main__":
    raise SystemExit(main())
