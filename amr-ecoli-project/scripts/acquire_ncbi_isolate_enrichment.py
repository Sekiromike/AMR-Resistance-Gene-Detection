"""Acquire a versioned NCBI isolate-enrichment overlay for a frozen AST snapshot.

This command never reselects records from the mutable NCBI AST table.  It reads
the exact target accessions from a previously frozen isolate CSV, verifies that
the parent AST, isolate, and acquisition-manifest hashes still match, and then
queries only the NCBI Pathogen Detection ``isolates`` table.

The output deliberately uses ``source_*`` and ``*_raw`` names.  Collection
organization, SRA submission center, disease text, BioProject, geography, and
epidemiological type are useful provenance or audit context, but none of them
is silently promoted to AST method, laboratory site, surveillance network, or
clinical indication.

``--dry-run`` is offline: it validates the parent snapshot and writes the exact
SQL plus a planned manifest without requiring Google Cloud credentials, a
query project, or the ``bq`` CLI.  Live execution sends SQL through standard
input, is compatible with an unbilled BigQuery Sandbox project, and installs
outputs only after fail-closed identity validation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import uuid
from datetime import datetime, timezone
from pathlib import Path, PureWindowsPath
from typing import Any, Mapping, Sequence


ISOLATES_TABLE = "ncbi-pathogen-detect.pdbrowser.isolates"
SOURCE_DOCUMENTATION = (
    "https://www.ncbi.nlm.nih.gov/pathogens/docs/isolates_gcp/"
)

QUERY_NAME = "ncbi_isolate_enrichment_query.sql"
OUTPUT_NAME = "ncbi_isolate_enrichment.csv"
MANIFEST_NAME = "ncbi_isolate_enrichment_manifest.json"
SOURCE_METADATA_NAME = "ncbi_isolate_enrichment_source_metadata.json"
JOB_METADATA_NAME = "ncbi_isolate_enrichment_job_metadata.json"
PARENT_AMENDMENT_NAME = "ncbi_acquisition_manifest_amendment.json"

EXPECTED_SCIENTIFIC_NAME = "Escherichia coli"
EXPECTED_SPECIES_TAXONOMY_ID = "562"
DEFAULT_MAXIMUM_BYTES_BILLED = 10_000_000_000

TARGET_RE = re.compile(r"^PDT\d+\.\d+$")
BIOSAMPLE_RE = re.compile(r"^SAM[NED][A-Z]?\d+$")
ASSEMBLY_RE = re.compile(r"^GC[AF]_\d+\.\d+$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
JOB_ID_RE = re.compile(r"^[A-Za-z0-9_-]{1,1024}$")

OUTPUT_COLUMNS = (
    "target_acc",
    "source_target_acc",
    "source_match_count",
    "source_biosample_accession",
    "source_assembly_accession",
    "source_scientific_name",
    "source_species_taxonomy_id",
    "sra_run_accessions_raw",
    "wgs_master_accession",
    "wgs_accession_prefix",
    "sra_library_layout_raw",
    "sra_platform_raw",
    "sra_release_date_raw",
    "submitted_assembly_level_raw",
    "submitted_assembly_method_raw",
    "submitted_assembly_length_bp",
    "submitted_assembly_contig_count",
    "submitted_assembly_contig_n50",
    "alternative_isolate_identifiers_json",
    "source_bioproject_accession",
    "source_collection_date_raw",
    "source_geographic_location_raw",
    "collected_by_raw",
    "sra_submission_center_raw",
    "epidemiological_type_raw",
    "host_disease_raw",
    "outbreak_raw",
    "source_isolate_checksum",
)

SOURCE_SCHEMA = {
    "target_acc": ("STRING", "NULLABLE"),
    "biosample_acc": ("STRING", "NULLABLE"),
    "asm_acc": ("STRING", "NULLABLE"),
    "scientific_name": ("STRING", "NULLABLE"),
    "species_taxid": ("INTEGER", "NULLABLE"),
    "Run": ("STRING", "NULLABLE"),
    "wgs_master_acc": ("STRING", "NULLABLE"),
    "wgs_acc_prefix": ("STRING", "NULLABLE"),
    "LibraryLayout": ("STRING", "NULLABLE"),
    "Platform": ("STRING", "NULLABLE"),
    "sra_release_date": ("STRING", "NULLABLE"),
    "asm_level": ("STRING", "NULLABLE"),
    "assembly_method": ("STRING", "NULLABLE"),
    "asm_stats_length_bp": ("INTEGER", "NULLABLE"),
    "asm_stats_n_contig": ("INTEGER", "NULLABLE"),
    "asm_stats_contig_n50": ("INTEGER", "NULLABLE"),
    "isolate_identifiers": ("STRING", "REPEATED"),
    "bioproject_acc": ("STRING", "NULLABLE"),
    "collection_date": ("STRING", "NULLABLE"),
    "geo_loc_name": ("STRING", "NULLABLE"),
    "collected_by": ("STRING", "NULLABLE"),
    "sra_center": ("STRING", "NULLABLE"),
    "epi_type": ("STRING", "NULLABLE"),
    "host_disease": ("STRING", "NULLABLE"),
    "outbreak": ("STRING", "NULLABLE"),
    "checksum": ("STRING", "NULLABLE"),
}

FIELD_SEMANTICS: dict[str, dict[str, list[str]]] = {
    "target_acc": {
        "allowed": ["exact join to the hash-verified parent isolate target"],
        "prohibited": ["replacement by another identifier"],
    },
    "source_target_acc": {
        "allowed": ["source-side identity consistency audit"],
        "prohibited": ["identity repair when it differs from target_acc"],
    },
    "source_match_count": {
        "allowed": ["fail-closed one-to-one source cardinality validation"],
        "prohibited": ["selection of one convenient row when count is not one"],
    },
    "source_biosample_accession": {
        "allowed": ["exact BioSample identity consistency audit"],
        "prohibited": ["identity repair when it differs from the frozen BioSample"],
    },
    "source_assembly_accession": {
        "allowed": ["candidate versioned genome linkage after parent-identity review"],
        "prohibited": ["silent replacement of the frozen assembly accession"],
    },
    "sra_run_accessions_raw": {
        "allowed": ["discovery of candidate read accessions for a separately versioned read workflow"],
        "prohibited": ["use as an assembly accession", "silent splitting or normalization"],
    },
    "wgs_master_accession": {
        "allowed": ["sequence-discovery provenance"],
        "prohibited": ["use as a versioned assembly accession"],
    },
    "wgs_accession_prefix": {
        "allowed": ["sequence-discovery provenance"],
        "prohibited": ["use as a versioned assembly accession"],
    },
    "sra_library_layout_raw": {
        "allowed": ["read-acquisition planning"],
        "prohibited": ["use as an AST method"],
    },
    "sra_platform_raw": {
        "allowed": ["sequencing provenance"],
        "prohibited": ["use as an AST platform or AST method"],
    },
    "sra_release_date_raw": {
        "allowed": ["sequence-release provenance"],
        "prohibited": ["use as collection date or AST testing date"],
    },
    "sra_submission_center_raw": {
        "allowed": ["sequence-submission provenance candidate"],
        "prohibited": ["laboratory site", "AST laboratory", "surveillance network"],
    },
    "collected_by_raw": {
        "allowed": ["collection-provenance candidate requiring source verification"],
        "prohibited": ["laboratory site", "AST laboratory", "surveillance network"],
    },
    "host_disease_raw": {
        "allowed": ["source-context audit"],
        "prohibited": ["clinical indication", "breakpoint indication"],
    },
    "epidemiological_type_raw": {
        "allowed": ["source-context audit"],
        "prohibited": ["clinical indication", "laboratory site"],
    },
    "source_geographic_location_raw": {
        "allowed": ["geographic provenance and drift audit"],
        "prohibited": ["laboratory site", "clinical indication"],
    },
    "source_bioproject_accession": {
        "allowed": ["source provenance and drift audit"],
        "prohibited": ["laboratory site", "external membership without prespecification"],
    },
    "source_collection_date_raw": {
        "allowed": ["temporal provenance with submitted precision preserved"],
        "prohibited": ["laboratory site", "imputation of missing month or day"],
    },
    "alternative_isolate_identifiers_json": {
        "allowed": ["duplicate-review evidence"],
        "prohibited": ["identity replacement or deduplication without corroboration"],
    },
    "outbreak_raw": {
        "allowed": ["duplicate and relatedness review signal"],
        "prohibited": ["genomic cluster or deduplication assignment without sequence analysis"],
    },
    "submitted_assembly_level_raw": {
        "allowed": ["download triage and source audit"],
        "prohibited": ["replacement of rerun genome QC or species confirmation"],
    },
    "submitted_assembly_method_raw": {
        "allowed": ["submitted assembly provenance"],
        "prohibited": ["AST method", "replacement of rerun genome QC"],
    },
    "submitted_assembly_length_bp": {
        "allowed": ["download triage and source audit"],
        "prohibited": ["replacement of rerun genome QC"],
    },
    "submitted_assembly_contig_count": {
        "allowed": ["download triage and source audit"],
        "prohibited": ["replacement of rerun genome QC"],
    },
    "submitted_assembly_contig_n50": {
        "allowed": ["download triage and source audit"],
        "prohibited": ["replacement of rerun genome QC or species confirmation"],
    },
    "source_scientific_name": {
        "allowed": ["source identity consistency audit"],
        "prohibited": ["replacement of computational species confirmation"],
    },
    "source_species_taxonomy_id": {
        "allowed": ["source identity consistency audit"],
        "prohibited": ["replacement of computational species confirmation"],
    },
    "source_isolate_checksum": {
        "allowed": ["upstream row-version and drift audit"],
        "prohibited": ["genomic sequence identity or file-integrity claim"],
    },
}

INTENTIONALLY_UNPOPULATED_CANONICAL_FIELDS = (
    "ast_method",
    "ast_testing_date",
    "disk_content",
    "breakpoint_version",
    "site",
    "clinical_indication",
    "surveillance_network",
    "species_method",
    "genome_qc_status",
    "deduplication_group",
    "lineage_group",
    "genomic_cluster",
    "evaluation_split",
)


class EnrichmentError(RuntimeError):
    """Raised when enrichment cannot proceed without ambiguity."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def normalize_utc_timestamp(value: str) -> str:
    text = value.strip()
    if not text:
        raise EnrichmentError("The source as-of timestamp must not be blank.")
    try:
        parsed = datetime.fromisoformat(text[:-1] + "+00:00" if text.endswith("Z") else text)
    except ValueError as exc:
        raise EnrichmentError(f"Invalid ISO-8601 source as-of timestamp: {value!r}") from exc
    if parsed.tzinfo is None:
        raise EnrichmentError("The source as-of timestamp must include an explicit UTC offset.")
    return parsed.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_text(text: str) -> str:
    return sha256_bytes(text.encode("utf-8"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def target_list_sha256(targets: Sequence[str]) -> str:
    return sha256_text("".join(f"{target}\n" for target in sorted(targets)))


def executable_identity(executable: str) -> str:
    if "\\" in executable:
        return PureWindowsPath(executable).name
    return Path(executable).name


def display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(Path.cwd().resolve()).as_posix()
    except ValueError:
        return path.name


def _load_json(path: Path, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise EnrichmentError(f"{label} does not exist: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as exc:
        raise EnrichmentError(f"{label} is not valid JSON: {path}") from exc
    if not isinstance(payload, dict):
        raise EnrichmentError(f"{label} must contain a JSON object: {path}")
    return payload


def _verify_artifact(path: Path, expected_sha256: str, label: str) -> str:
    if not path.is_file():
        raise EnrichmentError(f"Parent manifest references missing {label}: {path}")
    if not SHA256_RE.fullmatch(expected_sha256):
        raise EnrichmentError(
            f"Parent manifest has an invalid {label} SHA-256: {expected_sha256!r}"
        )
    observed = sha256_file(path)
    if observed != expected_sha256:
        raise EnrichmentError(
            f"Parent {label} hash mismatch: expected {expected_sha256}, observed {observed}."
        )
    return observed


def csv_data_row_count(path: Path) -> int:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            next(reader)
        except StopIteration as exc:
            raise EnrichmentError(f"Parent CSV is empty: {path}") from exc
        return sum(1 for _ in reader)


def read_frozen_isolates(path: Path) -> tuple[list[str], dict[str, dict[str, str]], int]:
    if not path.is_file():
        raise EnrichmentError(f"Parent isolate CSV does not exist: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise EnrichmentError("Parent isolate CSV has duplicate column names.")
        required = {"target_acc", "biosample_acc", "asm_acc"}
        missing = sorted(required.difference(fields))
        if missing:
            raise EnrichmentError(f"Parent isolate CSV lacks required columns: {missing}")
        rows = list(reader)
    if not rows:
        raise EnrichmentError("Parent isolate CSV contains no data rows.")

    identities: dict[str, dict[str, str]] = {}
    duplicates: list[str] = []
    for row_number, row in enumerate(rows, start=2):
        target = (row.get("target_acc") or "").strip()
        biosample = (row.get("biosample_acc") or "").strip()
        assembly = (row.get("asm_acc") or "").strip()
        if not TARGET_RE.fullmatch(target):
            raise EnrichmentError(
                f"Invalid or missing frozen target_acc at parent row {row_number}: {target!r}"
            )
        if not BIOSAMPLE_RE.fullmatch(biosample):
            raise EnrichmentError(
                f"Invalid or missing frozen biosample_acc at parent row {row_number}: {biosample!r}"
            )
        if assembly and not ASSEMBLY_RE.fullmatch(assembly):
            raise EnrichmentError(
                f"Invalid frozen asm_acc at parent row {row_number}: {assembly!r}"
            )
        if target in identities:
            duplicates.append(target)
            continue
        identities[target] = {
            "biosample_acc": biosample,
            "asm_acc": assembly,
            "scientific_name": (row.get("scientific_name") or "").strip(),
            "species_taxid": (row.get("species_taxid") or "").strip(),
            "bioproject_acc": (row.get("bioproject_acc") or "").strip(),
            "collection_date": (row.get("collection_date") or "").strip(),
            "geo_loc_name": (row.get("geo_loc_name") or "").strip(),
            "source_row_number": str(row_number),
        }
        if identities[target]["scientific_name"] != EXPECTED_SCIENTIFIC_NAME:
            raise EnrichmentError(
                f"Frozen isolate {target} does not have exact scientific_name "
                f"{EXPECTED_SCIENTIFIC_NAME!r}."
            )
        if identities[target]["species_taxid"] != EXPECTED_SPECIES_TAXONOMY_ID:
            raise EnrichmentError(
                f"Frozen isolate {target} does not have exact species_taxid "
                f"{EXPECTED_SPECIES_TAXONOMY_ID!r}."
            )
    if duplicates:
        examples = ", ".join(sorted(set(duplicates))[:5])
        raise EnrichmentError(
            f"Frozen target_acc values are not unique; duplicate examples: {examples}"
        )
    targets = sorted(identities)
    if len(targets) != len(rows):
        raise EnrichmentError("Frozen target count does not equal parent isolate row count.")
    return targets, identities, len(rows)


def _manifest_output(manifest: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    try:
        output = manifest["retrieval"]["outputs"][name]  # type: ignore[index]
    except (KeyError, TypeError) as exc:
        raise EnrichmentError(
            f"Parent acquisition manifest lacks retrieval.outputs.{name}."
        ) from exc
    if not isinstance(output, Mapping):
        raise EnrichmentError(f"Parent retrieval.outputs.{name} must be an object.")
    return output


def _expected_sha256(output: Mapping[str, Any], label: str) -> str:
    value = str(output.get("sha256", "")).strip().lower()
    if not SHA256_RE.fullmatch(value):
        raise EnrichmentError(f"Parent manifest has an invalid {label} SHA-256: {value!r}")
    return value


def _manifest_artifact(base_dir: Path, record: Mapping[str, Any], label: str) -> Path:
    filename = str(record.get("file", "")).strip()
    if not filename or Path(filename).name != filename:
        raise EnrichmentError(
            f"Parent manifest {label} must reference one plain artifact filename."
        )
    return base_dir / filename


def _verify_parent_amendment(
    *,
    amendment_path: Path | None,
    manifest_path: Path,
    manifest: Mapping[str, Any],
    verified: Mapping[str, str],
) -> dict[str, str] | None:
    retrieval = manifest.get("retrieval") or {}
    required_corrections: dict[str, str] = {}
    selection = manifest.get("selection") or {}
    antibiotics = selection.get("antibiotics_exact_case_insensitive")
    if isinstance(antibiotics, list):
        expected_parameter = ",".join(str(value) for value in antibiotics)
        observed_parameter = str(
            ((manifest.get("queries") or {}).get("parameters") or {}).get(
                "antibiotics_csv", ""
            )
        )
        if observed_parameter != expected_parameter:
            required_corrections["/queries/parameters/antibiotics_csv"] = (
                expected_parameter
            )
    executable = str(((retrieval.get("client") or {}).get("executable")) or "")
    portable_executable = executable.replace("\\", "/").rsplit("/", 1)[-1]
    if executable and executable != portable_executable:
        required_corrections["/retrieval/client/executable"] = portable_executable

    if amendment_path is None:
        if required_corrections:
            raise EnrichmentError(
                "Parent acquisition manifest has provenance serialization errors but "
                "no verified append-only amendment."
            )
        return None

    amendment = _load_json(amendment_path, "Parent acquisition amendment")
    if amendment.get("amendment_status") != "complete":
        raise EnrichmentError("Parent acquisition amendment status is not 'complete'.")
    original = amendment.get("original_manifest") or {}
    manifest_hash = sha256_file(manifest_path)
    if original.get("sha256") != manifest_hash:
        raise EnrichmentError(
            "Parent acquisition amendment does not authenticate the supplied manifest."
        )
    if original.get("file") != manifest_path.name:
        raise EnrichmentError(
            "Parent acquisition amendment references a different manifest filename."
        )
    unchanged = amendment.get("unchanged_artifacts") or {}
    for amendment_key, verified_key in (
        ("ast_sha256", "ast"),
        ("isolates_sha256", "isolates"),
        ("source_table_metadata_sha256", "source_table_metadata"),
    ):
        if unchanged.get(amendment_key) != verified.get(verified_key):
            raise EnrichmentError(
                f"Parent acquisition amendment disagrees with verified {amendment_key}."
            )
    corrections: dict[str, str] = {}
    for item in amendment.get("corrections", []):
        if not isinstance(item, Mapping):
            raise EnrichmentError("Parent amendment corrections must be JSON objects.")
        pointer = str(item.get("json_pointer", ""))
        if not pointer or pointer in corrections:
            raise EnrichmentError(
                "Parent amendment corrections have a missing or duplicate JSON pointer."
            )
        corrections[pointer] = str(item.get("corrected_value", ""))
    for pointer, expected in required_corrections.items():
        if corrections.get(pointer) != expected:
            raise EnrichmentError(
                f"Parent acquisition amendment lacks the required correction for {pointer}."
            )
    return {
        "file": display_path(amendment_path),
        "sha256": sha256_file(amendment_path),
    }


def load_parent_context(
    *,
    isolates_path: Path,
    ast_path: Path,
    acquisition_manifest_path: Path,
    amendment_path: Path | None = None,
) -> dict[str, Any]:
    manifest = _load_json(acquisition_manifest_path, "Parent acquisition manifest")
    if manifest.get("status") != "complete":
        raise EnrichmentError("Parent acquisition manifest status is not 'complete'.")
    base_dir = acquisition_manifest_path.parent
    isolate_output = _manifest_output(manifest, "isolates")
    ast_output = _manifest_output(manifest, "ast")
    expected_isolate_hash = _expected_sha256(isolate_output, "isolate output")
    expected_ast_hash = _expected_sha256(ast_output, "AST output")
    if str(isolate_output.get("file", "")) != isolates_path.name:
        raise EnrichmentError("Parent isolate filename does not match the acquisition manifest.")
    if str(ast_output.get("file", "")) != ast_path.name:
        raise EnrichmentError("Parent AST filename does not match the acquisition manifest.")
    observed_isolate_hash = _verify_artifact(
        isolates_path, expected_isolate_hash, "isolate output"
    )
    observed_ast_hash = _verify_artifact(ast_path, expected_ast_hash, "AST output")

    targets, identities, isolate_rows = read_frozen_isolates(isolates_path)
    ast_rows = csv_data_row_count(ast_path)
    try:
        expected_isolate_rows = int(isolate_output["rows"])
        expected_ast_rows = int(ast_output["rows"])
    except (KeyError, TypeError, ValueError) as exc:
        raise EnrichmentError("Parent manifest output row counts are missing or invalid.") from exc
    if isolate_rows != expected_isolate_rows:
        raise EnrichmentError(
            "Parent isolate row count does not match the acquisition manifest."
        )
    if ast_rows != expected_ast_rows:
        raise EnrichmentError("Parent AST row count does not match the acquisition manifest.")
    for label, record, path in (
        ("isolate output", isolate_output, isolates_path),
        ("AST output", ast_output, ast_path),
    ):
        try:
            expected_bytes = int(record["bytes"])
        except (KeyError, TypeError, ValueError) as exc:
            raise EnrichmentError(
                f"Parent manifest {label} byte count is missing or invalid."
            ) from exc
        if path.stat().st_size != expected_bytes:
            raise EnrichmentError(f"Parent manifest {label} byte count does not match.")

    verified: dict[str, str] = {
        "ast": observed_ast_hash,
        "isolates": observed_isolate_hash,
    }
    verified_records: dict[str, dict[str, Any]] = {}
    for logical_name in ("ast", "isolates"):
        query_record = (manifest.get("queries") or {}).get(logical_name) or {}
        if not isinstance(query_record, Mapping):
            raise EnrichmentError(f"Parent {logical_name} query record is not an object.")
        query_path = _manifest_artifact(base_dir, query_record, f"{logical_name} query")
        query_hash = _verify_artifact(
            query_path,
            str(query_record.get("sha256", "")),
            f"{logical_name} query",
        )
        verified[f"query_{logical_name}"] = query_hash
        verified_records[f"{logical_name}_query"] = {
            "file": display_path(query_path),
            "sha256": query_hash,
        }

    metadata_record = (manifest.get("retrieval") or {}).get("source_table_metadata") or {}
    if not isinstance(metadata_record, Mapping):
        raise EnrichmentError("Parent source-table-metadata record is not an object.")
    metadata_path = _manifest_artifact(
        base_dir, metadata_record, "source table metadata"
    )
    metadata_hash = _verify_artifact(
        metadata_path,
        str(metadata_record.get("sha256", "")),
        "source table metadata",
    )
    verified["source_table_metadata"] = metadata_hash
    verified_records["source_table_metadata"] = {
        "file": display_path(metadata_path),
        "sha256": metadata_hash,
    }
    amendment_record = _verify_parent_amendment(
        amendment_path=amendment_path,
        manifest_path=acquisition_manifest_path,
        manifest=manifest,
        verified=verified,
    )

    parent_paths: dict[str, Path] = {
        "isolates": isolates_path,
        "ast": ast_path,
        "acquisition_manifest": acquisition_manifest_path,
        "query_ast": base_dir / str((manifest.get("queries") or {})["ast"]["file"]),
        "query_isolates": base_dir
        / str((manifest.get("queries") or {})["isolates"]["file"]),
        "source_table_metadata": metadata_path,
    }
    if amendment_path is not None:
        parent_paths["acquisition_amendment"] = amendment_path
    parent_hashes = {
        "isolates": observed_isolate_hash,
        "ast": observed_ast_hash,
        "acquisition_manifest": sha256_file(acquisition_manifest_path),
        "query_ast": verified["query_ast"],
        "query_isolates": verified["query_isolates"],
        "source_table_metadata": metadata_hash,
    }
    if amendment_record is not None:
        parent_hashes["acquisition_amendment"] = amendment_record["sha256"]

    return {
        "paths": parent_paths,
        "hashes": parent_hashes,
        "targets": targets,
        "identities": identities,
        "target_list_sha256": target_list_sha256(targets),
        "rows": {"isolates": isolate_rows, "ast": expected_ast_rows},
        "manifest_record": {
            "acquisition_manifest": {
                "file": display_path(acquisition_manifest_path),
                "sha256": sha256_file(acquisition_manifest_path),
            },
            "acquisition_amendment": amendment_record,
            "ast": {
                "file": display_path(ast_path),
                "rows": expected_ast_rows,
                "bytes": ast_path.stat().st_size,
                "sha256": observed_ast_hash,
            },
            "isolates": {
                "file": display_path(isolates_path),
                "rows": isolate_rows,
                "bytes": isolates_path.stat().st_size,
                "sha256": observed_isolate_hash,
            },
            "verified_supporting_artifacts": verified_records,
        },
    }


def assert_parent_unchanged(context: Mapping[str, Any]) -> None:
    paths = context["paths"]
    hashes = context["hashes"]
    for name, path in paths.items():
        if not path.is_file() or sha256_file(path) != hashes[name]:
            raise EnrichmentError(
                f"Parent {name} changed during enrichment acquisition; refusing outputs."
            )


def build_query(targets: Sequence[str], parent_isolate_sha256: str) -> str:
    if not targets:
        raise EnrichmentError("Cannot build enrichment SQL without frozen targets.")
    normalized = sorted(targets)
    if len(normalized) != len(set(normalized)):
        raise EnrichmentError("Cannot build enrichment SQL from duplicate target accessions.")
    invalid = [target for target in normalized if not TARGET_RE.fullmatch(target)]
    if invalid:
        raise EnrichmentError(f"Cannot build enrichment SQL from invalid target: {invalid[0]!r}")
    if not SHA256_RE.fullmatch(parent_isolate_sha256):
        raise EnrichmentError("Parent isolate SHA-256 is invalid.")
    literals = ",\n    ".join(f"'{target}'" for target in normalized)
    return f"""-- NCBI isolate enrichment for an immutable parent target list.
-- Parent isolate CSV SHA-256: {parent_isolate_sha256}
-- Context fields remain raw and must not fill AST method, site, or clinical indication.
WITH frozen_targets AS (
  SELECT target_acc
  FROM UNNEST([
    {literals}
  ]) AS target_acc
)
SELECT
  f.target_acc AS target_acc,
  i.target_acc AS source_target_acc,
  COUNT(i.target_acc) OVER (PARTITION BY f.target_acc) AS source_match_count,
  i.biosample_acc AS source_biosample_accession,
  i.asm_acc AS source_assembly_accession,
  i.scientific_name AS source_scientific_name,
  i.species_taxid AS source_species_taxonomy_id,
  i.Run AS sra_run_accessions_raw,
  i.wgs_master_acc AS wgs_master_accession,
  i.wgs_acc_prefix AS wgs_accession_prefix,
  i.LibraryLayout AS sra_library_layout_raw,
  i.Platform AS sra_platform_raw,
  i.sra_release_date AS sra_release_date_raw,
  i.asm_level AS submitted_assembly_level_raw,
  i.assembly_method AS submitted_assembly_method_raw,
  i.asm_stats_length_bp AS submitted_assembly_length_bp,
  i.asm_stats_n_contig AS submitted_assembly_contig_count,
  i.asm_stats_contig_n50 AS submitted_assembly_contig_n50,
  TO_JSON_STRING(i.isolate_identifiers) AS alternative_isolate_identifiers_json,
  i.bioproject_acc AS source_bioproject_accession,
  i.collection_date AS source_collection_date_raw,
  i.geo_loc_name AS source_geographic_location_raw,
  i.collected_by AS collected_by_raw,
  i.sra_center AS sra_submission_center_raw,
  i.epi_type AS epidemiological_type_raw,
  i.host_disease AS host_disease_raw,
  i.outbreak AS outbreak_raw,
  i.checksum AS source_isolate_checksum
FROM frozen_targets AS f
LEFT JOIN `{ISOLATES_TABLE}` AS i
  FOR SYSTEM_TIME AS OF @source_as_of_utc
  ON i.target_acc = f.target_acc
ORDER BY f.target_acc, i.biosample_acc, i.asm_acc
"""


def base_manifest(
    *,
    context: Mapping[str, Any],
    query: str,
    generated_at: str,
    source_as_of_utc: str,
    maximum_bytes_billed: int,
) -> dict[str, Any]:
    return {
        "schema_version": "1.0.0",
        "artifact_type": "ncbi_isolate_enrichment_overlay",
        "status": "planned",
        "generated_at_utc": generated_at,
        "source_as_of_utc": source_as_of_utc,
        "source": {
            "provider": "NCBI Pathogen Detection",
            "table": ISOLATES_TABLE,
            "documentation": SOURCE_DOCUMENTATION,
            "upstream_release_status": "alpha",
            "upstream_update_frequency": "daily",
        },
        "software": {
            "script": {
                "file": display_path(Path(__file__)),
                "sha256": sha256_file(Path(__file__)),
            },
            "python_version": platform.python_version(),
        },
        "parent_snapshot": context["manifest_record"],
        "selection": {
            "basis": "literal sorted target_acc values from the hash-verified parent isolate CSV",
            "target_count": len(context["targets"]),
            "target_list_sha256": context["target_list_sha256"],
            "target_list_serialization": "UTF-8, one sorted target_acc plus LF per record",
            "ast_table_reselected": False,
        },
        "query": {
            "file": QUERY_NAME,
            "sha256": sha256_text(query),
            "dialect": "GoogleSQL",
            "parameters": {"source_as_of_utc": source_as_of_utc},
            "maximum_bytes_billed": maximum_bytes_billed,
            "temporal_consistency": "FOR SYSTEM_TIME AS OF an explicit TIMESTAMP parameter",
        },
        "field_semantics": {
            "fields": FIELD_SEMANTICS,
            "intentionally_unpopulated_canonical_fields": list(
                INTENTIONALLY_UNPOPULATED_CANONICAL_FIELDS
            ),
        },
        "coverage": None,
        "reconciliation": None,
        "retrieval": None,
        "limitations": [
            "The enrichment is a later versioned overlay and is not the same temporal snapshot as its parent.",
            "Run and WGS identifiers are not assembly accessions.",
            "Submitted assembly statistics do not replace rerun genome QC or species confirmation.",
            "collected_by and sra_center do not establish the AST laboratory or study site.",
            "host_disease, epi_type, specimen source, BioProject, geography, and year do not establish clinical indication.",
            "No contextual field is authorized as a model feature by this acquisition.",
        ],
    }


def validate_source_metadata(metadata: Mapping[str, Any]) -> None:
    if metadata.get("type") != "TABLE":
        raise EnrichmentError(
            "The isolates source is not a native BigQuery TABLE; time-travel semantics are unavailable."
        )
    try:
        fields = metadata["schema"]["fields"]  # type: ignore[index]
    except (KeyError, TypeError) as exc:
        raise EnrichmentError("BigQuery source metadata lacks schema.fields.") from exc
    observed: dict[str, tuple[str, str]] = {}
    for field in fields:
        if not isinstance(field, Mapping) or "name" not in field:
            continue
        observed[str(field["name"])] = (
            str(field.get("type", "")),
            str(field.get("mode", "NULLABLE")),
        )
    missing = sorted(set(SOURCE_SCHEMA).difference(observed))
    if missing:
        raise EnrichmentError(f"NCBI isolates schema lacks enrichment fields: {missing}")
    changed = {
        name: {"expected": SOURCE_SCHEMA[name], "observed": observed[name]}
        for name in SOURCE_SCHEMA
        if observed[name] != SOURCE_SCHEMA[name]
    }
    if changed:
        raise EnrichmentError(f"NCBI isolates enrichment schema changed: {changed}")


def _coverage(rows: Sequence[Mapping[str, str]]) -> dict[str, dict[str, int]]:
    result: dict[str, dict[str, int]] = {}
    for field in OUTPUT_COLUMNS:
        values = [(row.get(field) or "").strip() for row in rows]
        nonmissing = [value for value in values if value]
        result[field] = {
            "present": len(nonmissing),
            "missing": len(values) - len(nonmissing),
            "unique_nonmissing": len(set(nonmissing)),
        }
    return result


def _context_comparison(source: str, parent: str) -> str:
    if not source and not parent:
        return "both_missing"
    if source == parent:
        return "same"
    if not source and parent:
        return "source_missing"
    if source and not parent:
        return "source_added"
    return "changed"


def validate_enrichment_output(
    path: Path, context: Mapping[str, Any]
) -> dict[str, Any]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise EnrichmentError("Enrichment output has duplicate column names.")
        if tuple(fields) != OUTPUT_COLUMNS:
            missing = sorted(set(OUTPUT_COLUMNS).difference(fields))
            extra = sorted(set(fields).difference(OUTPUT_COLUMNS))
            raise EnrichmentError(
                f"Enrichment output schema mismatch; missing={missing}, extra={extra}."
            )
        rows = [dict(row) for row in reader]

    expected_targets = set(context["targets"])
    if len(rows) != len(expected_targets):
        raise EnrichmentError(
            f"Enrichment row count conflict: expected {len(expected_targets)}, observed {len(rows)}."
        )
    observed_targets = [(row.get("target_acc") or "").strip() for row in rows]
    if len(set(observed_targets)) != len(observed_targets):
        raise EnrichmentError("Enrichment output contains duplicate target_acc rows.")
    if set(observed_targets) != expected_targets:
        missing = sorted(expected_targets.difference(observed_targets))[:5]
        extra = sorted(set(observed_targets).difference(expected_targets))[:5]
        raise EnrichmentError(
            f"Enrichment target-set conflict; missing={missing}, extra={extra}."
        )

    reconciliation = {
        "row_count": len(rows),
        "source_match_count_one": 0,
        "biosample_exact": 0,
        "assembly": {"same": 0, "recovered": 0, "both_missing": 0},
        "context_drift": {
            field: {
                "same": 0,
                "both_missing": 0,
                "source_missing": 0,
                "source_added": 0,
                "changed": 0,
            }
            for field in ("bioproject_acc", "collection_date", "geo_loc_name")
        },
    }
    identities = context["identities"]
    for row in rows:
        target = (row.get("target_acc") or "").strip()
        source_target = (row.get("source_target_acc") or "").strip()
        match_text = (row.get("source_match_count") or "").strip()
        try:
            match_count = int(match_text)
        except ValueError as exc:
            raise EnrichmentError(
                f"Invalid source_match_count for {target}: {match_text!r}"
            ) from exc
        if match_count != 1:
            raise EnrichmentError(
                f"Source-match conflict for {target}: expected 1, observed {match_count}."
            )
        reconciliation["source_match_count_one"] += 1
        if source_target != target:
            raise EnrichmentError(
                f"Source target identity conflict for {target}: {source_target!r}."
            )
        parent = identities[target]
        source_biosample = (row.get("source_biosample_accession") or "").strip()
        if source_biosample != parent["biosample_acc"]:
            raise EnrichmentError(
                f"BioSample identity conflict for {target}: parent={parent['biosample_acc']!r}, "
                f"source={source_biosample!r}."
            )
        reconciliation["biosample_exact"] += 1

        source_assembly = (row.get("source_assembly_accession") or "").strip()
        parent_assembly = parent["asm_acc"]
        if source_assembly and not ASSEMBLY_RE.fullmatch(source_assembly):
            raise EnrichmentError(
                f"Invalid source assembly accession for {target}: {source_assembly!r}."
            )
        if parent_assembly:
            if source_assembly != parent_assembly:
                raise EnrichmentError(
                    f"Assembly identity conflict for {target}: parent={parent_assembly!r}, "
                    f"source={source_assembly!r}."
                )
            reconciliation["assembly"]["same"] += 1
        elif source_assembly:
            reconciliation["assembly"]["recovered"] += 1
        else:
            reconciliation["assembly"]["both_missing"] += 1

        source_name = (row.get("source_scientific_name") or "").strip()
        if source_name != EXPECTED_SCIENTIFIC_NAME:
            raise EnrichmentError(
                f"Scientific-name identity conflict for {target}: expected "
                f"{EXPECTED_SCIENTIFIC_NAME!r}, observed {source_name!r}."
            )
        source_taxid = (row.get("source_species_taxonomy_id") or "").strip()
        if source_taxid != EXPECTED_SPECIES_TAXONOMY_ID:
            raise EnrichmentError(
                f"Species-taxonomy identity conflict for {target}: expected "
                f"{EXPECTED_SPECIES_TAXONOMY_ID!r}, observed {source_taxid!r}."
            )

        alternate_ids = (row.get("alternative_isolate_identifiers_json") or "").strip()
        if alternate_ids:
            try:
                parsed_ids = json.loads(alternate_ids)
            except json.JSONDecodeError as exc:
                raise EnrichmentError(
                    f"Invalid alternative-isolate-identifiers JSON for {target}."
                ) from exc
            if not isinstance(parsed_ids, list):
                raise EnrichmentError(
                    f"Alternative isolate identifiers are not a JSON list for {target}."
                )

        comparisons = {
            "bioproject_acc": (
                (row.get("source_bioproject_accession") or "").strip(),
                parent["bioproject_acc"],
            ),
            "collection_date": (
                (row.get("source_collection_date_raw") or "").strip(),
                parent["collection_date"],
            ),
            "geo_loc_name": (
                (row.get("source_geographic_location_raw") or "").strip(),
                parent["geo_loc_name"],
            ),
        }
        for field, (source_value, parent_value) in comparisons.items():
            status = _context_comparison(source_value, parent_value)
            reconciliation["context_drift"][field][status] += 1

    return {
        "rows": len(rows),
        "coverage": _coverage(rows),
        "reconciliation": reconciliation,
    }


def command_environment() -> dict[str, str]:
    environment = os.environ.copy()
    environment["CLOUDSDK_CORE_DISABLE_PROMPTS"] = "1"
    environment["CLOUDSDK_ENCODING"] = "utf-8"
    environment["PYTHONIOENCODING"] = "utf-8"
    return environment


def run_capture(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(command),
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="strict",
        env=command_environment(),
    )


def require_success(result: subprocess.CompletedProcess[str], context: str) -> None:
    if result.returncode == 0:
        return
    detail = (result.stderr or result.stdout or "no diagnostic from bq").strip()
    raise EnrichmentError(f"{context} failed. bq said: {detail}")


def bq_query_command(
    *,
    bq: str,
    project: str,
    source_as_of_utc: str,
    job_id: str,
    maximum_bytes_billed: int,
) -> list[str]:
    if not JOB_ID_RE.fullmatch(job_id):
        raise EnrichmentError(f"Invalid BigQuery job ID: {job_id!r}")
    if maximum_bytes_billed <= 0:
        raise EnrichmentError("maximum_bytes_billed must be a positive integer.")
    return [
        bq,
        f"--project_id={project}",
        f"--job_id={job_id}",
        "query",
        "--use_legacy_sql=false",
        "--format=csv",
        "--max_rows=2147483647",
        f"--maximum_bytes_billed={maximum_bytes_billed}",
        "--label=workflow:amr_ncbi_isolate_enrichment",
        f"--parameter=source_as_of_utc:TIMESTAMP:{source_as_of_utc}",
    ]


def run_query_to_file(command: Sequence[str], destination: Path, sql: str) -> None:
    with destination.open("w", encoding="utf-8", newline="") as output:
        result = subprocess.run(
            list(command),
            check=False,
            input=sql,
            stdout=output,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="strict",
            env=command_environment(),
        )
    if result.returncode != 0:
        detail = (result.stderr or "no diagnostic from bq").strip()
        raise EnrichmentError(
            "BigQuery enrichment export failed; no completed artifact was installed. "
            f"Verify authentication, project permissions, sandbox quota, and time-travel access. bq said: {detail}"
        )


def write_text_atomic(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        handle.write(content)
        temporary = Path(handle.name)
    temporary.replace(path)


def write_json_atomic(path: Path, payload: Any) -> None:
    write_text_atomic(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def check_destinations(output_dir: Path) -> None:
    if output_dir.exists():
        raise EnrichmentError(
            f"Refusing to overwrite isolate-enrichment output directory {output_dir}. "
            "Choose a new versioned output directory."
        )


def ensure_separate_output(output_dir: Path, context: Mapping[str, Any]) -> None:
    resolved = output_dir.resolve()
    parent_directories = {path.resolve().parent for path in context["paths"].values()}
    if any(parent == resolved or parent in resolved.parents for parent in parent_directories):
        raise EnrichmentError(
            "The enrichment output directory must be outside every parent snapshot directory."
        )


def write_plan(output_dir: Path, query: str, manifest: Mapping[str, Any]) -> None:
    check_destinations(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{output_dir.name}.staging-", dir=output_dir.parent
    ) as name:
        staging = Path(name)
        write_text_atomic(staging / QUERY_NAME, query)
        write_json_atomic(staging / MANIFEST_NAME, manifest)
        check_destinations(output_dir)
        staging.replace(output_dir)


def _show_table_metadata(bq: str, project: str) -> dict[str, Any]:
    result = run_capture(
        [bq, f"--project_id={project}", "show", "--format=prettyjson", ISOLATES_TABLE.replace(".", ":", 1)]
    )
    require_success(result, f"BigQuery preflight for {ISOLATES_TABLE}")
    try:
        metadata = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise EnrichmentError("bq returned invalid source-table metadata JSON.") from exc
    if not isinstance(metadata, dict):
        raise EnrichmentError("bq source-table metadata is not a JSON object.")
    validate_source_metadata(metadata)
    return metadata


def _show_job_metadata(
    bq: str, project: str, job_id: str, location: str
) -> dict[str, Any]:
    command = [bq, f"--project_id={project}"]
    if location:
        command.append(f"--location={location}")
    command.extend(["--format=prettyjson", "show", "--job=true", job_id])
    result = run_capture(command)
    require_success(result, f"BigQuery job-metadata retrieval for {job_id}")
    try:
        metadata = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise EnrichmentError("bq returned invalid job metadata JSON.") from exc
    if not isinstance(metadata, dict):
        raise EnrichmentError("bq job metadata is not a JSON object.")
    return metadata


def validate_job_metadata(
    metadata: Mapping[str, Any],
    *,
    project: str,
    job_id: str,
    query: str,
    source_as_of_utc: str,
    maximum_bytes_billed: int,
) -> dict[str, Any]:
    reference = metadata.get("jobReference") or {}
    status = metadata.get("status") or {}
    statistics = metadata.get("statistics") or {}
    query_statistics = statistics.get("query") or {}
    configuration = (metadata.get("configuration") or {}).get("query") or {}
    if reference.get("jobId") != job_id or reference.get("projectId") != project:
        raise EnrichmentError("BigQuery job metadata does not match the submitted job identity.")
    if status.get("state") != "DONE" or status.get("errorResult"):
        raise EnrichmentError("BigQuery job metadata does not describe a successful DONE job.")
    if configuration.get("query") != query:
        raise EnrichmentError("BigQuery job metadata query text differs from the frozen SQL.")
    if configuration.get("useLegacySql") is not False:
        raise EnrichmentError("BigQuery job metadata does not confirm GoogleSQL execution.")
    if str(configuration.get("maximumBytesBilled", "")) != str(maximum_bytes_billed):
        raise EnrichmentError("BigQuery job metadata has a different maximum-bytes-billed cap.")
    parameters = configuration.get("queryParameters") or []
    observed_parameter = None
    for parameter in parameters:
        if not isinstance(parameter, Mapping):
            continue
        if parameter.get("name") == "source_as_of_utc":
            observed_parameter = str(
                ((parameter.get("parameterValue") or {}).get("value")) or ""
            )
    normalized_observed = (
        normalize_utc_timestamp(observed_parameter) if observed_parameter else None
    )
    if normalized_observed != source_as_of_utc:
        raise EnrichmentError(
            "BigQuery job metadata does not confirm the source_as_of_utc parameter."
        )

    def nonnegative_integer(mapping: Mapping[str, Any], key: str) -> int:
        try:
            value = int(mapping[key])
        except (KeyError, TypeError, ValueError) as exc:
            raise EnrichmentError(f"BigQuery job metadata lacks a valid {key} value.") from exc
        if value < 0:
            raise EnrichmentError(f"BigQuery job metadata has a negative {key} value.")
        return value

    processed = nonnegative_integer(query_statistics, "totalBytesProcessed")
    billed = nonnegative_integer(query_statistics, "totalBytesBilled")
    if billed > maximum_bytes_billed:
        raise EnrichmentError("BigQuery reports billed bytes above the declared safety cap.")
    return {
        "job_id": job_id,
        "project_id": project,
        "location": str(reference.get("location", "")),
        "state": "DONE",
        "total_bytes_processed": processed,
        "total_bytes_billed": billed,
        "cache_hit": bool(query_statistics.get("cacheHit", False)),
        "created_at_ms": statistics.get("creationTime"),
        "started_at_ms": statistics.get("startTime"),
        "ended_at_ms": statistics.get("endTime"),
    }


def acquire(
    *,
    output_dir: Path,
    project: str,
    source_as_of_utc: str,
    query: str,
    manifest: dict[str, Any],
    context: Mapping[str, Any],
    maximum_bytes_billed: int,
) -> dict[str, Any]:
    check_destinations(output_dir)
    bq = shutil.which("bq")
    if bq is None:
        raise EnrichmentError(
            "The Google Cloud bq CLI is not installed or not on PATH. Use --dry-run for an offline plan."
        )

    retrieval_started = utc_now()
    job_id = f"amr_ncbi_isolate_enrichment_{uuid.uuid4().hex}"
    metadata_before = _show_table_metadata(bq, project)
    source_location = str(metadata_before.get("location", "")).strip()
    version_result = run_capture([bq, "version"])
    bq_version = (
        (version_result.stdout or version_result.stderr).strip()
        if version_result.returncode == 0
        else None
    )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{output_dir.name}.staging-", dir=output_dir.parent
    ) as name:
        staging = Path(name)
        output_path = staging / OUTPUT_NAME
        run_query_to_file(
            bq_query_command(
                bq=bq,
                project=project,
                source_as_of_utc=source_as_of_utc,
                job_id=job_id,
                maximum_bytes_billed=maximum_bytes_billed,
            ),
            output_path,
            query,
        )
        validation = validate_enrichment_output(output_path, context)
        job_metadata = _show_job_metadata(
            bq, project, job_id, source_location
        )
        job_summary = validate_job_metadata(
            job_metadata,
            project=project,
            job_id=job_id,
            query=query,
            source_as_of_utc=source_as_of_utc,
            maximum_bytes_billed=maximum_bytes_billed,
        )
        metadata_after = _show_table_metadata(bq, project)
        assert_parent_unchanged(context)

        metadata_payload = {
            "retrieved_at_utc": retrieval_started,
            "source_as_of_utc": source_as_of_utc,
            "table": ISOLATES_TABLE,
            "preflight": metadata_before,
            "postflight": metadata_after,
        }
        metadata_path = staging / SOURCE_METADATA_NAME
        write_json_atomic(metadata_path, metadata_payload)
        job_metadata_path = staging / JOB_METADATA_NAME
        job_metadata_payload = {
            "schema_version": "1.0.0",
            "artifact_type": "sanitized_bigquery_job_provenance",
            "job": job_summary,
            "query": {
                "sha256": sha256_text(query),
                "dialect": "GoogleSQL",
                "source_as_of_utc": source_as_of_utc,
                "maximum_bytes_billed": maximum_bytes_billed,
            },
            "privacy_note": (
                "Allowlisted from bq job metadata; account email, self-links, and "
                "other client-private fields are intentionally omitted."
            ),
        }
        write_json_atomic(job_metadata_path, job_metadata_payload)

        completed = utc_now()
        manifest["status"] = "complete"
        manifest["coverage"] = validation["coverage"]
        manifest["reconciliation"] = validation["reconciliation"]
        manifest["retrieval"] = {
            "started_at_utc": retrieval_started,
            "completed_at_utc": completed,
            "query_project": project,
            "client": {
                "basename": executable_identity(bq),
                "version": bq_version,
            },
            "job": {
                **job_summary,
                "metadata_file": JOB_METADATA_NAME,
                "metadata_sha256": sha256_file(job_metadata_path),
                "maximum_bytes_billed": maximum_bytes_billed,
            },
            "source_table_metadata": {
                "file": SOURCE_METADATA_NAME,
                "sha256": sha256_file(metadata_path),
            },
            "output": {
                "file": OUTPUT_NAME,
                "rows": validation["rows"],
                "bytes": output_path.stat().st_size,
                "sha256": sha256_file(output_path),
            },
        }
        write_text_atomic(staging / QUERY_NAME, query)
        write_json_atomic(staging / MANIFEST_NAME, manifest)
        check_destinations(output_dir)
        staging.replace(output_dir)
    return manifest


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--parent-isolates",
        type=Path,
        default=Path("data/source/ncbi_ast_snapshot/ncbi_isolates.csv"),
    )
    parser.add_argument(
        "--parent-ast",
        type=Path,
        default=Path("data/source/ncbi_ast_snapshot/ncbi_ast.csv"),
    )
    parser.add_argument(
        "--parent-manifest",
        type=Path,
        default=Path("data/source/ncbi_ast_snapshot/ncbi_acquisition_manifest.json"),
    )
    parser.add_argument(
        "--parent-amendment",
        type=Path,
        help=(
            "Append-only parent acquisition correction; auto-detected beside "
            "--parent-manifest when present."
        ),
    )
    parser.add_argument(
        "--source-as-of-utc",
        required=True,
        help="Explicit ISO-8601 timestamp for BigQuery FOR SYSTEM_TIME AS OF.",
    )
    parser.add_argument(
        "--maximum-bytes-billed",
        type=int,
        default=DEFAULT_MAXIMUM_BYTES_BILLED,
        help="Fail the query above this byte cap (default: 10,000,000,000).",
    )
    parser.add_argument(
        "--project",
        help="Google Cloud query project (or GOOGLE_CLOUD_PROJECT); billing is optional in BigQuery Sandbox.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Offline plan only: validate parents and emit exact SQL plus planned manifest.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        source_as_of_utc = normalize_utc_timestamp(args.source_as_of_utc)
        if args.maximum_bytes_billed <= 0:
            raise EnrichmentError("--maximum-bytes-billed must be positive.")
        amendment_path = args.parent_amendment
        inferred_amendment = args.parent_manifest.parent / PARENT_AMENDMENT_NAME
        if amendment_path is None and inferred_amendment.is_file():
            amendment_path = inferred_amendment
        context = load_parent_context(
            isolates_path=args.parent_isolates,
            ast_path=args.parent_ast,
            acquisition_manifest_path=args.parent_manifest,
            amendment_path=amendment_path,
        )
        ensure_separate_output(args.output_dir, context)
        check_destinations(args.output_dir)
        query = build_query(context["targets"], context["hashes"]["isolates"])
        manifest = base_manifest(
            context=context,
            query=query,
            generated_at=utc_now(),
            source_as_of_utc=source_as_of_utc,
            maximum_bytes_billed=args.maximum_bytes_billed,
        )
        if args.dry_run:
            write_plan(args.output_dir, query, manifest)
            print(f"Offline isolate-enrichment plan written to {args.output_dir}")
            return 0

        project = (args.project or os.environ.get("GOOGLE_CLOUD_PROJECT") or "").strip()
        if not project:
            raise EnrichmentError(
                "A Google Cloud query project is required for live acquisition. "
                "Pass --project PROJECT_ID or set GOOGLE_CLOUD_PROJECT."
            )
        acquire(
            output_dir=args.output_dir,
            project=project,
            source_as_of_utc=source_as_of_utc,
            query=query,
            manifest=manifest,
            context=context,
            maximum_bytes_billed=args.maximum_bytes_billed,
        )
        print(f"NCBI isolate enrichment complete: {args.output_dir / MANIFEST_NAME}")
        return 0
    except (EnrichmentError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
