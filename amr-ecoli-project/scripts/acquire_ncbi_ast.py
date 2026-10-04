"""Acquire raw E. coli AST and isolate metadata from NCBI BigQuery.

The command is deliberately an acquisition layer, not a harmonization layer.
It preserves every AST column returned by NCBI and exports a documented set of
scalar Isolates Browser columns. It never fills missing assay, breakpoint, or
sample metadata.

``--dry-run`` is completely offline: it writes the exact SQL and a planned
provenance manifest without requiring Google Cloud credentials or the ``bq``
CLI. A real acquisition requires an explicit query project, an authenticated
``bq`` installation, and writes all outputs atomically. The query project may
use the no-billing BigQuery Sandbox when the workload stays within its quotas.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


AST_TABLE = "ncbi-pathogen-detect.pdbrowser.ast"
ISOLATES_TABLE = "ncbi-pathogen-detect.pdbrowser.isolates"
SOURCE_DOCUMENTATION = {
    "ast": "https://www.ncbi.nlm.nih.gov/pathogens/docs/ast_gcp/",
    "ast_context": "https://www.ncbi.nlm.nih.gov/pathogens/docs/ast/",
    "bigquery_access": "https://www.ncbi.nlm.nih.gov/pathogens/docs/getting_started_bigquery/",
    "isolates": "https://www.ncbi.nlm.nih.gov/pathogens/docs/isolates_gcp/",
}

AST_QUERY_NAME = "ncbi_ast_query.sql"
ISOLATES_QUERY_NAME = "ncbi_isolates_query.sql"
MANIFEST_NAME = "ncbi_acquisition_manifest.json"
SOURCE_METADATA_NAME = "ncbi_source_tables.json"
AST_OUTPUT_NAME = "ncbi_ast.csv"
ISOLATES_OUTPUT_NAME = "ncbi_isolates.csv"


class AcquisitionError(RuntimeError):
    """Raised when acquisition cannot complete without ambiguous partial data."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


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


def normalize_antibiotics(values: Sequence[str]) -> list[str]:
    normalized: list[str] = []
    for value in values:
        drug = value.strip().lower()
        if not drug:
            raise ValueError("Antibiotic names must not be blank.")
        if any(delimiter in drug for delimiter in (",", "|")) or any(
            character in drug for character in "\r\n\x00"
        ):
            raise ValueError(
                "Antibiotic names must not contain ',', '|', newlines, or NUL bytes."
            )
        if drug not in normalized:
            normalized.append(drug)
    return sorted(normalized)


def serialize_antibiotics_parameter(antibiotics: Sequence[str]) -> str:
    """Serialize normalized antibiotic names for the shared BigQuery parameter."""
    return ",".join(antibiotics)


def executable_identity(executable: str) -> str:
    """Return a portable executable basename without retaining a local path."""
    return executable.replace("\\", "/").rsplit("/", 1)[-1]


def build_queries() -> dict[str, str]:
    predicate = """LOWER(TRIM(a.scientific_name)) = LOWER(@species)
  AND (
    @antibiotics_csv = ''
    OR LOWER(TRIM(a.antibiotic)) IN UNNEST(SPLIT(@antibiotics_csv, ','))
  )"""
    ast = f"""-- Raw, submitter-provided AST rows. No breakpoint or method is inferred.
SELECT a.*
FROM `{AST_TABLE}` AS a
WHERE {predicate}
ORDER BY a.target_acc, a.antibiotic
"""
    isolates = f"""-- Scalar isolate metadata for exactly the AST targets selected above.
WITH selected_targets AS (
  SELECT DISTINCT a.target_acc
  FROM `{AST_TABLE}` AS a
  WHERE {predicate}
)
SELECT
  i.target_acc,
  i.biosample_acc,
  i.asm_acc,
  i.scientific_name,
  i.species_taxid,
  i.taxid,
  i.strain,
  i.bioproject_acc,
  i.collection_date,
  i.geo_loc_name,
  i.host,
  i.isolation_source,
  i.source_type,
  i.creation_date AS target_creation_date,
  i.erd_group,
  i.minsame,
  i.mindiff
FROM `{ISOLATES_TABLE}` AS i
INNER JOIN selected_targets AS selected USING (target_acc)
ORDER BY i.target_acc
"""
    return {"ast": ast, "isolates": isolates}


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


def count_csv_rows(path: Path, required_columns: set[str]) -> int:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise AcquisitionError(f"BigQuery produced an empty file: {path}") from exc
        missing = required_columns - set(header)
        if missing:
            raise AcquisitionError(
                f"BigQuery output {path.name} lacks required source columns: {sorted(missing)}. "
                "The upstream NCBI schema may have changed."
            )
        return sum(1 for _ in reader)


def base_manifest(
    *, species: str, antibiotics: Sequence[str], queries: dict[str, str], generated_at: str
) -> dict[str, Any]:
    return {
        "schema_version": "1.0.0",
        "status": "planned",
        "generated_at_utc": generated_at,
        "source": {
            "provider": "NCBI Pathogen Detection",
            "tables": [AST_TABLE, ISOLATES_TABLE],
            "documentation": SOURCE_DOCUMENTATION,
            "upstream_release_status": "alpha",
            "upstream_update_frequency": "daily",
            "scientific_caveat": (
                "AST values and field relationships are submitter-provided and are not "
                "verified by NCBI beyond basic quality control."
            ),
        },
        "selection": {
            "species_exact_case_insensitive": species,
            "antibiotics_exact_case_insensitive": list(antibiotics),
        },
        "queries": {
            "ast": {"file": AST_QUERY_NAME, "sha256": sha256_text(queries["ast"])},
            "isolates": {
                "file": ISOLATES_QUERY_NAME,
                "sha256": sha256_text(queries["isolates"]),
            },
            "dialect": "Google Standard SQL",
            "parameters": {
                "species": species,
                "antibiotics_csv": serialize_antibiotics_parameter(antibiotics),
            },
        },
        "retrieval": None,
        "limitations": [
            "No missing assay method, breakpoint authority/version, or measurement is inferred.",
            "A testing-standard name does not establish a breakpoint edition.",
            "The source tables are mutable daily snapshots; retained exports and hashes are the audit record.",
        ],
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
        errors="replace",
        env=command_environment(),
    )


def require_success(result: subprocess.CompletedProcess[str], context: str) -> None:
    if result.returncode == 0:
        return
    detail = (result.stderr or result.stdout or "no diagnostic from bq").strip()
    raise AcquisitionError(
        f"{context} failed. Authenticate non-interactively with `gcloud auth login "
        "--no-launch-browser` (or configure an approved service account), verify the query "
        f"project, and retry. bq said: {detail}"
    )


def bq_query_command(
    *, bq: str, project: str, species: str, antibiotics: Sequence[str]
) -> list[str]:
    return [
        bq,
        f"--project_id={project}",
        "query",
        "--use_legacy_sql=false",
        "--format=csv",
        "--max_rows=2147483647",
        f"--parameter=species:STRING:{species}",
        "--parameter=antibiotics_csv:STRING:"
        f"{serialize_antibiotics_parameter(antibiotics)}",
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
            errors="replace",
            env=command_environment(),
        )
    if result.returncode != 0:
        detail = (result.stderr or "no diagnostic from bq").strip()
        raise AcquisitionError(
            "BigQuery export failed; no completed raw export was installed. Verify Google Cloud "
            f"authentication, project permissions, quota, and the current NCBI schema. bq said: {detail}"
        )


def check_destinations(output_dir: Path, force: bool) -> None:
    names = {
        AST_QUERY_NAME,
        ISOLATES_QUERY_NAME,
        MANIFEST_NAME,
        SOURCE_METADATA_NAME,
        AST_OUTPUT_NAME,
        ISOLATES_OUTPUT_NAME,
    }
    existing = sorted(name for name in names if (output_dir / name).exists())
    if existing and not force:
        raise AcquisitionError(
            f"Refusing to overwrite existing acquisition artifacts in {output_dir}: {existing}. "
            "Choose a new output directory or pass --force explicitly."
        )


def write_plan(output_dir: Path, queries: dict[str, str], manifest: dict[str, Any]) -> None:
    write_text_atomic(output_dir / AST_QUERY_NAME, queries["ast"])
    write_text_atomic(output_dir / ISOLATES_QUERY_NAME, queries["isolates"])
    write_json_atomic(output_dir / MANIFEST_NAME, manifest)


def acquire(
    *,
    output_dir: Path,
    project: str,
    species: str,
    antibiotics: Sequence[str],
    queries: dict[str, str],
    manifest: dict[str, Any],
) -> dict[str, Any]:
    bq = shutil.which("bq")
    if bq is None:
        raise AcquisitionError(
            "The Google Cloud `bq` CLI is not installed or not on PATH. Install the Google Cloud "
            "CLI and authenticate it, or use --dry-run to emit an offline query plan."
        )

    retrieval_started = utc_now()
    metadata: dict[str, Any] = {"retrieved_at_utc": retrieval_started, "tables": {}}
    for logical_name, table in (("ast", AST_TABLE), ("isolates", ISOLATES_TABLE)):
        result = run_capture(
            [bq, f"--project_id={project}", "show", "--format=prettyjson", table.replace(".", ":", 1)]
        )
        require_success(result, f"BigQuery preflight for {table}")
        try:
            metadata["tables"][logical_name] = json.loads(result.stdout)
        except json.JSONDecodeError as exc:
            raise AcquisitionError(f"bq returned invalid table metadata JSON for {table}.") from exc

    version_result = run_capture([bq, "version"])
    bq_version = (version_result.stdout or version_result.stderr).strip() if version_result.returncode == 0 else None

    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ncbi-acquisition-", dir=output_dir) as temporary_name:
        temporary = Path(temporary_name)
        metadata_path = temporary / SOURCE_METADATA_NAME
        write_json_atomic(metadata_path, metadata)
        ast_path = temporary / AST_OUTPUT_NAME
        isolates_path = temporary / ISOLATES_OUTPUT_NAME
        run_query_to_file(
            bq_query_command(
                bq=bq,
                project=project,
                species=species,
                antibiotics=antibiotics,
            ),
            ast_path,
            queries["ast"],
        )
        run_query_to_file(
            bq_query_command(
                bq=bq,
                project=project,
                species=species,
                antibiotics=antibiotics,
            ),
            isolates_path,
            queries["isolates"],
        )

        ast_rows = count_csv_rows(ast_path, {"target_acc", "scientific_name", "antibiotic"})
        isolate_rows = count_csv_rows(isolates_path, {"target_acc", "biosample_acc", "asm_acc"})
        completed = utc_now()
        manifest["status"] = "complete"
        manifest["retrieval"] = {
            "started_at_utc": retrieval_started,
            "completed_at_utc": completed,
            "query_project": project,
            "client": {
                "executable": executable_identity(bq),
                "version": bq_version,
            },
            "source_table_metadata": {
                "file": SOURCE_METADATA_NAME,
                "sha256": sha256_file(metadata_path),
            },
            "outputs": {
                "ast": {
                    "file": AST_OUTPUT_NAME,
                    "rows": ast_rows,
                    "bytes": ast_path.stat().st_size,
                    "sha256": sha256_file(ast_path),
                },
                "isolates": {
                    "file": ISOLATES_OUTPUT_NAME,
                    "rows": isolate_rows,
                    "bytes": isolates_path.stat().st_size,
                    "sha256": sha256_file(isolates_path),
                },
            },
        }

        # Install data only after both queries and structural checks succeed.
        metadata_path.replace(output_dir / SOURCE_METADATA_NAME)
        ast_path.replace(output_dir / AST_OUTPUT_NAME)
        isolates_path.replace(output_dir / ISOLATES_OUTPUT_NAME)

    write_plan(output_dir, queries, manifest)
    return manifest


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        # Frozen source snapshots live under data/source/ and are never
        # overwritten; a second retrieval requires a new dated directory.
        default=Path("data/source/ncbi_ast"),
        help="Directory for raw exports, exact SQL, source metadata, and manifest.",
    )
    parser.add_argument(
        "--project",
        help="Google Cloud query project (or set GOOGLE_CLOUD_PROJECT). Required for execution; billing is optional within BigQuery Sandbox quotas.",
    )
    parser.add_argument("--species", default="Escherichia coli", help="Exact case-insensitive species name.")
    parser.add_argument(
        "--antibiotic",
        action="append",
        default=[],
        help="Exact case-insensitive antibiotic name; repeat as needed. Empty means all antibiotics.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Offline plan only: emit SQL and manifest without checking bq, auth, or project.",
    )
    parser.add_argument("--force", action="store_true", help="Explicitly replace artifacts in output-dir.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        species = args.species.strip()
        if not species:
            raise AcquisitionError("--species must not be blank.")
        antibiotics = normalize_antibiotics(args.antibiotic)
        queries = build_queries()
        generated_at = utc_now()
        manifest = base_manifest(
            species=species,
            antibiotics=antibiotics,
            queries=queries,
            generated_at=generated_at,
        )
        check_destinations(args.output_dir, args.force)
        if args.dry_run:
            args.output_dir.mkdir(parents=True, exist_ok=True)
            write_plan(args.output_dir, queries, manifest)
            print(f"Offline query plan written to {args.output_dir}")
            return 0

        project = (args.project or os.environ.get("GOOGLE_CLOUD_PROJECT") or "").strip()
        if not project:
            raise AcquisitionError(
                "A Google Cloud query project is required for acquisition. Pass --project PROJECT_ID "
                "or set GOOGLE_CLOUD_PROJECT. Use --dry-run to emit SQL without credentials."
            )
        acquire(
            output_dir=args.output_dir,
            project=project,
            species=species,
            antibiotics=antibiotics,
            queries=queries,
            manifest=manifest,
        )
        print(f"NCBI acquisition complete: {args.output_dir / MANIFEST_NAME}")
        return 0
    except (AcquisitionError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
