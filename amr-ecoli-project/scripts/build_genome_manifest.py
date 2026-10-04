"""Build reproducible genome request, staging, and analysis manifests.

This module deliberately separates three operations:

``request``
    Deduplicate the curated isolate-antibiotic cohort into one record per
    assembly while rejecting ambiguous isolate/assembly mappings.
``stage``
    Validate an extracted and rehydrated NCBI Datasets genome package and
    hard-link (or copy) one genomic FASTA per requested assembly.
``finalize``
    Consolidate QUAST, seven-locus MLST, and AMRFinderPlus results with file
    checksums, explicit exclusions, and non-excluding typing warnings.

MLST is reported only as a typing covariate. It is not used here as a
whole-genome relatedness cluster or as a phylogeny.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ASSEMBLY_RE = re.compile(r"^GC[AF]_\d+\.\d+$")
FASTA_SUFFIXES = ("_genomic.fna", "_genomic.fna.gz")

REQUEST_FIELDS = [
    "assembly_accession",
    "isolate_id",
    "biosample_accession",
    "n_ast_records",
]
STAGE_FIELDS = [
    "assembly_accession",
    "isolate_id",
    "biosample_accession",
    "n_ast_records",
    "source_fasta",
    "staged_fasta",
    "fasta_sha256",
    "fasta_size_bytes",
    "fasta_contigs",
    "fasta_total_bases",
    "fasta_ambiguous_bases",
    "stage_status",
    "exclusion_reason",
]
FINAL_FIELDS = [
    "assembly_accession",
    "isolate_id",
    "biosample_accession",
    "fasta_path",
    "fasta_sha256",
    "quast_exit_code",
    "quast_total_length",
    "quast_contigs",
    "quast_n50",
    "quast_gc_percent",
    "quast_ns_per_100kb",
    "mlst_exit_code",
    "mlst_scheme",
    "mlst_sequence_type",
    "mlst_status",
    "mlst_score",
    "mlst_alleles",
    "amrfinder_exit_code",
    "amrfinder_hit_count",
    "quast_report_sha256",
    "mlst_report_sha256",
    "amrfinder_report_sha256",
    "eligible_for_modeling",
    "overall_status",
    "exclusion_reasons",
    "warning_reasons",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_tsv(path: Path, rows: Iterable[Mapping[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _join_values(values: Iterable[str]) -> str:
    return ";".join(sorted({value.strip() for value in values if value and value.strip()}))


def build_requests(cohort: Path) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """Return unambiguous one-assembly records and row-level exclusions."""
    with cohort.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"isolate_id", "assembly_accession"}
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Cohort is missing required columns: {missing}")
        rows = list(reader)

    exclusions: list[dict[str, str]] = []
    valid_rows: list[dict[str, str]] = []
    for index, row in enumerate(rows, start=2):
        isolate = (row.get("isolate_id") or "").strip()
        accession = (row.get("assembly_accession") or "").strip()
        if not isolate:
            exclusions.append(
                {
                    "assembly_accession": accession,
                    "isolate_id": isolate,
                    "reason": "MISSING_ISOLATE_ID",
                    "source_row": str(index),
                }
            )
        elif not ASSEMBLY_RE.fullmatch(accession):
            exclusions.append(
                {
                    "assembly_accession": accession,
                    "isolate_id": isolate,
                    "reason": "INVALID_ASSEMBLY_ACCESSION",
                    "source_row": str(index),
                }
            )
        else:
            row["_source_row"] = str(index)
            valid_rows.append(row)

    isolate_accessions: dict[str, set[str]] = {}
    accession_isolates: dict[str, set[str]] = {}
    for row in valid_rows:
        isolate = row["isolate_id"].strip()
        accession = row["assembly_accession"].strip()
        isolate_accessions.setdefault(isolate, set()).add(accession)
        accession_isolates.setdefault(accession, set()).add(isolate)

    ambiguous_isolates = {key for key, values in isolate_accessions.items() if len(values) != 1}
    ambiguous_accessions = {key for key, values in accession_isolates.items() if len(values) != 1}
    accepted_rows: list[dict[str, str]] = []
    for row in valid_rows:
        isolate = row["isolate_id"].strip()
        accession = row["assembly_accession"].strip()
        reasons = []
        if isolate in ambiguous_isolates:
            reasons.append("ISOLATE_MULTIPLE_ASSEMBLIES")
        if accession in ambiguous_accessions:
            reasons.append("ASSEMBLY_MULTIPLE_ISOLATES")
        if reasons:
            exclusions.append(
                {
                    "assembly_accession": accession,
                    "isolate_id": isolate,
                    "reason": ";".join(reasons),
                    "source_row": row["_source_row"],
                }
            )
        else:
            accepted_rows.append(row)

    grouped: dict[str, list[dict[str, str]]] = {}
    for row in accepted_rows:
        grouped.setdefault(row["assembly_accession"].strip(), []).append(row)

    requests: list[dict[str, Any]] = []
    for accession in sorted(grouped):
        group = grouped[accession]
        requests.append(
            {
                "assembly_accession": accession,
                "isolate_id": _join_values(row["isolate_id"] for row in group),
                "biosample_accession": _join_values(row.get("biosample_accession", "") for row in group),
                "n_ast_records": len(group),
            }
        )
    exclusions.sort(key=lambda row: (row["assembly_accession"], row["isolate_id"], row["reason"]))
    return requests, exclusions


def write_request_outputs(
    cohort: Path,
    accessions: Path,
    mapping: Path,
    exclusions_path: Path,
    manifest_path: Path,
) -> dict[str, Any]:
    requests, exclusions = build_requests(cohort)
    accessions.parent.mkdir(parents=True, exist_ok=True)
    accessions.write_text(
        "".join(f"{row['assembly_accession']}\n" for row in requests), encoding="ascii"
    )
    _write_tsv(mapping, requests, REQUEST_FIELDS)
    _write_tsv(
        exclusions_path,
        exclusions,
        ["assembly_accession", "isolate_id", "reason", "source_row"],
    )
    payload = {
        "schema_version": "1.0.0",
        "operation": "genome_request",
        "generated_at_utc": utc_now(),
        "cohort": str(cohort),
        "cohort_sha256": sha256_file(cohort),
        "requested_assemblies": len(requests),
        "excluded_rows": len(exclusions),
        "outputs": {
            "accessions": {"path": str(accessions), "sha256": sha256_file(accessions)},
            "mapping": {"path": str(mapping), "sha256": sha256_file(mapping)},
            "exclusions": {"path": str(exclusions_path), "sha256": sha256_file(exclusions_path)},
        },
    }
    _write_json(manifest_path, payload)
    return payload


def _walk_scalars(value: Any) -> Iterable[str]:
    if isinstance(value, str):
        yield value
    elif isinstance(value, Mapping):
        for child in value.values():
            yield from _walk_scalars(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk_scalars(child)


def report_accessions(report_path: Path) -> set[str]:
    accessions: set[str] = set()
    with report_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(f"Invalid JSON in {report_path} line {line_number}: {error}") from error
            accessions.update(value for value in _walk_scalars(record) if ASSEMBLY_RE.fullmatch(value))
    return accessions


def inspect_fasta(path: Path) -> dict[str, int]:
    opener = gzip.open if path.suffix == ".gz" else open
    contigs = 0
    total_bases = 0
    ambiguous_bases = 0
    saw_sequence = False
    with opener(path, "rt", encoding="ascii", errors="strict") as handle:
        for line_number, line in enumerate(handle, start=1):
            text = line.strip()
            if not text:
                continue
            if text.startswith(">"):
                if len(text) == 1:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                contigs += 1
                continue
            if contigs == 0:
                raise ValueError(f"Sequence precedes the first FASTA header at line {line_number}")
            sequence = text.upper()
            invalid = set(sequence).difference("ACGTNRYKMSWBDHV.-")
            if invalid:
                raise ValueError(f"Invalid FASTA characters at line {line_number}: {sorted(invalid)}")
            sequence = sequence.replace("-", "").replace(".", "")
            total_bases += len(sequence)
            ambiguous_bases += sum(base not in "ACGT" for base in sequence)
            saw_sequence = saw_sequence or bool(sequence)
    if contigs == 0 or not saw_sequence or total_bases == 0:
        raise ValueError("FASTA contains no sequence")
    return {
        "contigs": contigs,
        "total_bases": total_bases,
        "ambiguous_bases": ambiguous_bases,
    }


def _safe_package_file(package_root: Path, candidate: Path) -> bool:
    try:
        candidate.resolve(strict=True).relative_to(package_root.resolve(strict=True))
    except (FileNotFoundError, ValueError):
        return False
    return candidate.is_file() and not candidate.is_symlink()


def _stage_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if sha256_file(source) == sha256_file(destination):
            return
        destination.unlink()
    try:
        os.link(source, destination)
    except OSError:
        shutil.copy2(source, destination)


def stage_package(
    requests_path: Path,
    package_dir: Path,
    stage_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, str]], dict[str, Any]]:
    requests = _read_tsv(requests_path)
    data_dir = package_dir / "ncbi_dataset" / "data"
    report_path = data_dir / "assembly_data_report.jsonl"
    catalog_path = data_dir / "dataset_catalog.json"
    if not data_dir.is_dir():
        raise ValueError(f"Not an extracted NCBI Datasets genome package: missing {data_dir}")
    if not report_path.is_file():
        raise ValueError(f"NCBI Datasets assembly report is missing: {report_path}")
    reported = report_accessions(report_path)

    staged: list[dict[str, Any]] = []
    exclusions: list[dict[str, str]] = []
    for request in requests:
        accession = request["assembly_accession"]
        row: dict[str, Any] = {
            **request,
            "source_fasta": "",
            "staged_fasta": "",
            "fasta_sha256": "",
            "fasta_size_bytes": "",
            "fasta_contigs": "",
            "fasta_total_bases": "",
            "fasta_ambiguous_bases": "",
            "stage_status": "EXCLUDE",
            "exclusion_reason": "",
        }
        reason = ""
        assembly_dir = data_dir / accession
        candidates = []
        if assembly_dir.is_dir():
            candidates = sorted(
                path
                for path in assembly_dir.iterdir()
                if path.name.endswith(FASTA_SUFFIXES) and _safe_package_file(package_dir, path)
            )
        if accession not in reported:
            reason = "ACCESSION_ABSENT_FROM_ASSEMBLY_REPORT"
        elif not assembly_dir.is_dir():
            reason = "ACCESSION_DIRECTORY_MISSING"
        elif not candidates:
            reason = "GENOMIC_FASTA_MISSING"
        elif len(candidates) > 1:
            reason = "MULTIPLE_GENOMIC_FASTA_FILES"
        else:
            source = candidates[0]
            try:
                fasta_metrics = inspect_fasta(source)
            except (OSError, UnicodeError, ValueError) as error:
                reason = f"INVALID_GENOMIC_FASTA:{type(error).__name__}"
            else:
                suffix = ".fna.gz" if source.name.endswith(".fna.gz") else ".fna"
                destination = stage_dir / accession / f"{accession}{suffix}"
                _stage_file(source, destination)
                row.update(
                    {
                        "source_fasta": str(source),
                        "staged_fasta": str(destination),
                        "fasta_sha256": sha256_file(destination),
                        "fasta_size_bytes": destination.stat().st_size,
                        "fasta_contigs": fasta_metrics["contigs"],
                        "fasta_total_bases": fasta_metrics["total_bases"],
                        "fasta_ambiguous_bases": fasta_metrics["ambiguous_bases"],
                        "stage_status": "STAGED",
                    }
                )
        if reason:
            row["exclusion_reason"] = reason
            exclusions.append(
                {
                    "assembly_accession": accession,
                    "isolate_id": request.get("isolate_id", ""),
                    "reason": reason,
                }
            )
        staged.append(row)

    package_metadata = {
        "package_directory": str(package_dir),
        "assembly_report": str(report_path),
        "assembly_report_sha256": sha256_file(report_path),
        "dataset_catalog": str(catalog_path) if catalog_path.exists() else None,
        "dataset_catalog_sha256": sha256_file(catalog_path) if catalog_path.exists() else None,
        "reported_assembly_accessions": len(reported),
    }
    return staged, exclusions, package_metadata


def write_stage_outputs(
    requests: Path,
    package_dir: Path,
    stage_dir: Path,
    table: Path,
    exclusions_path: Path,
    manifest_path: Path,
) -> dict[str, Any]:
    staged, exclusions, package_metadata = stage_package(requests, package_dir, stage_dir)
    _write_tsv(table, staged, STAGE_FIELDS)
    _write_tsv(exclusions_path, exclusions, ["assembly_accession", "isolate_id", "reason"])
    payload = {
        "schema_version": "1.0.0",
        "operation": "ncbi_datasets_genome_staging",
        "generated_at_utc": utc_now(),
        "requests": {"path": str(requests), "sha256": sha256_file(requests)},
        "package": package_metadata,
        "requested_assemblies": len(staged),
        "staged_assemblies": sum(row["stage_status"] == "STAGED" for row in staged),
        "excluded_assemblies": len(exclusions),
        "outputs": {
            "table": {"path": str(table), "sha256": sha256_file(table)},
            "exclusions": {"path": str(exclusions_path), "sha256": sha256_file(exclusions_path)},
        },
    }
    _write_json(manifest_path, payload)
    return payload


def _read_exit_code(path: Path) -> int | None:
    try:
        return int(path.read_text(encoding="ascii").strip())
    except (FileNotFoundError, ValueError):
        return None


def parse_quast_report(path: Path) -> dict[str, float]:
    metrics: dict[str, float] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 2:
                continue
            key = row[0].strip()
            value = row[1].strip().replace(",", "")
            try:
                metrics[key] = float(value)
            except ValueError:
                continue
    return metrics


def parse_mlst_report(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one MLST row, found {len(rows)}")
    row = {str(key).strip().upper(): str(value).strip() for key, value in rows[0].items()}
    return {
        "scheme": row.get("SCHEME", ""),
        "sequence_type": row.get("ST", ""),
        "status": row.get("STATUS", ""),
        "score": row.get("SCORE", ""),
        "alleles": row.get("ALLELES", ""),
    }


def count_amrfinder_hits(path: Path) -> int:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise ValueError("AMRFinderPlus output is empty and has no header")
    return max(0, len(rows) - 1)


def _quast_value(metrics: Mapping[str, float], *keys: str) -> float | None:
    for key in keys:
        if key in metrics:
            return metrics[key]
    return None


def finalize_analyses(
    stage_table: Path,
    quast_root: Path,
    mlst_root: Path,
    amrfinder_root: Path,
    *,
    min_total_length: int,
    max_total_length: int,
    max_contigs: int,
    min_n50: int,
    max_ns_per_100kb: float,
) -> list[dict[str, Any]]:
    final_rows: list[dict[str, Any]] = []
    for staged in _read_tsv(stage_table):
        accession = staged["assembly_accession"]
        exclusions = [value for value in staged.get("exclusion_reason", "").split(";") if value]
        warnings: list[str] = []
        row: dict[str, Any] = {
            "assembly_accession": accession,
            "isolate_id": staged.get("isolate_id", ""),
            "biosample_accession": staged.get("biosample_accession", ""),
            "fasta_path": staged.get("staged_fasta", ""),
            "fasta_sha256": staged.get("fasta_sha256", ""),
        }
        if staged.get("stage_status") == "STAGED":
            quast_dir = quast_root / accession
            quast_report = quast_dir / "report.tsv"
            quast_exit = _read_exit_code(quast_dir / "exit_code.txt")
            row["quast_exit_code"] = "" if quast_exit is None else quast_exit
            metrics: dict[str, float] = {}
            if quast_exit != 0:
                exclusions.append("QUAST_EXECUTION_FAILED")
            elif not quast_report.is_file():
                exclusions.append("QUAST_REPORT_MISSING")
            else:
                metrics = parse_quast_report(quast_report)
                row["quast_report_sha256"] = sha256_file(quast_report)
            total_length = _quast_value(metrics, "Total length", "Total length (>= 0 bp)")
            contigs = _quast_value(metrics, "# contigs", "# contigs (>= 0 bp)")
            n50 = _quast_value(metrics, "N50")
            gc_percent = _quast_value(metrics, "GC (%)")
            ns_per_100kb = _quast_value(metrics, "# N's per 100 kbp")
            row.update(
                {
                    "quast_total_length": "" if total_length is None else int(total_length),
                    "quast_contigs": "" if contigs is None else int(contigs),
                    "quast_n50": "" if n50 is None else int(n50),
                    "quast_gc_percent": "" if gc_percent is None else gc_percent,
                    "quast_ns_per_100kb": "" if ns_per_100kb is None else ns_per_100kb,
                }
            )
            if quast_exit == 0:
                missing_metrics = [
                    name
                    for name, value in {
                        "total_length": total_length,
                        "contigs": contigs,
                        "n50": n50,
                        "ns_per_100kb": ns_per_100kb,
                    }.items()
                    if value is None
                ]
                if missing_metrics:
                    exclusions.append("QUAST_REQUIRED_METRICS_MISSING:" + ",".join(missing_metrics))
                else:
                    if total_length < min_total_length:
                        exclusions.append("ASSEMBLY_LENGTH_BELOW_MIN")
                    if total_length > max_total_length:
                        exclusions.append("ASSEMBLY_LENGTH_ABOVE_MAX")
                    if contigs > max_contigs:
                        exclusions.append("ASSEMBLY_CONTIG_COUNT_ABOVE_MAX")
                    if n50 < min_n50:
                        exclusions.append("ASSEMBLY_N50_BELOW_MIN")
                    if ns_per_100kb > max_ns_per_100kb:
                        exclusions.append("ASSEMBLY_NS_ABOVE_MAX")

            mlst_dir = mlst_root / accession
            mlst_report = mlst_dir / "mlst.tsv"
            mlst_exit = _read_exit_code(mlst_dir / "exit_code.txt")
            row["mlst_exit_code"] = "" if mlst_exit is None else mlst_exit
            if mlst_exit != 0:
                warnings.append("MLST_EXECUTION_FAILED")
            elif not mlst_report.is_file():
                warnings.append("MLST_REPORT_MISSING")
            else:
                try:
                    mlst = parse_mlst_report(mlst_report)
                except ValueError:
                    warnings.append("MLST_REPORT_INVALID")
                else:
                    row.update(
                        {
                            "mlst_scheme": mlst["scheme"],
                            "mlst_sequence_type": mlst["sequence_type"],
                            "mlst_status": mlst["status"],
                            "mlst_score": mlst["score"],
                            "mlst_alleles": mlst["alleles"],
                            "mlst_report_sha256": sha256_file(mlst_report),
                        }
                    )
                    if mlst["scheme"] != "ecoli_achtman_4":  # Achtman; "ecoli" is Pasteur in mlst >= 2.23
                        warnings.append("UNEXPECTED_MLST_SCHEME")
                    if mlst["status"] != "PERFECT":
                        warnings.append("MLST_NOT_PERFECT")

            amr_dir = amrfinder_root / accession
            amr_report = amr_dir / "amrfinder.tsv"
            amr_exit = _read_exit_code(amr_dir / "exit_code.txt")
            row["amrfinder_exit_code"] = "" if amr_exit is None else amr_exit
            if amr_exit != 0:
                exclusions.append("AMRFINDER_EXECUTION_FAILED")
            elif not amr_report.is_file():
                exclusions.append("AMRFINDER_REPORT_MISSING")
            else:
                try:
                    row["amrfinder_hit_count"] = count_amrfinder_hits(amr_report)
                except ValueError:
                    exclusions.append("AMRFINDER_REPORT_INVALID")
                else:
                    row["amrfinder_report_sha256"] = sha256_file(amr_report)

        exclusions = sorted(set(exclusions))
        warnings = sorted(set(warnings))
        row["eligible_for_modeling"] = "false" if exclusions else "true"
        row["overall_status"] = "EXCLUDE" if exclusions else ("WARN" if warnings else "PASS")
        row["exclusion_reasons"] = ";".join(exclusions)
        row["warning_reasons"] = ";".join(warnings)
        final_rows.append(row)
    return final_rows


def write_final_outputs(args: argparse.Namespace) -> dict[str, Any]:
    rows = finalize_analyses(
        args.stage_table,
        args.quast_root,
        args.mlst_root,
        args.amrfinder_root,
        min_total_length=args.min_total_length,
        max_total_length=args.max_total_length,
        max_contigs=args.max_contigs,
        min_n50=args.min_n50,
        max_ns_per_100kb=args.max_ns_per_100kb,
    )
    _write_tsv(args.output, rows, FINAL_FIELDS)
    exclusions = [
        {
            "assembly_accession": row["assembly_accession"],
            "isolate_id": row["isolate_id"],
            "reason": row["exclusion_reasons"],
        }
        for row in rows
        if row["exclusion_reasons"]
    ]
    _write_tsv(args.exclusions, exclusions, ["assembly_accession", "isolate_id", "reason"])
    tool_metadata: dict[str, Any] = {}
    for name, path in {
        "tool_versions": args.tool_versions,
        "mlst_scheme_info": args.mlst_scheme_info,
        "amrfinder_database_metadata": args.amrfinder_database_metadata,
        "amrfinder_database_checksums": args.amrfinder_database_checksums,
    }.items():
        tool_metadata[name] = {
            "path": str(path),
            "sha256": sha256_file(path) if path.is_file() else None,
        }
    payload = {
        "schema_version": "1.0.0",
        "operation": "genome_qc_typing_amr_annotation",
        "generated_at_utc": utc_now(),
        "status": "PASS" if rows and all(not row["exclusion_reasons"] for row in rows) else "REVIEW",
        "assembly_count": len(rows),
        "eligible_assembly_count": sum(row["eligible_for_modeling"] == "true" for row in rows),
        "excluded_assembly_count": len(exclusions),
        "quality_thresholds": {
            "min_total_length": args.min_total_length,
            "max_total_length": args.max_total_length,
            "max_contigs": args.max_contigs,
            "min_n50": args.min_n50,
            "max_ns_per_100kb": args.max_ns_per_100kb,
        },
        "typing_scope": (
            "Seven-locus MLST is a typing covariate only; it is not a whole-genome "
            "relatedness cluster and is not used as a phylogeny."
        ),
        "stage_table": {"path": str(args.stage_table), "sha256": sha256_file(args.stage_table)},
        "tool_and_database_provenance": tool_metadata,
        "outputs": {
            "table": {"path": str(args.output), "sha256": sha256_file(args.output)},
            "exclusions": {"path": str(args.exclusions), "sha256": sha256_file(args.exclusions)},
        },
    }
    _write_json(args.manifest, payload)
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    request = subparsers.add_parser("request", help="Build unique assembly requests from a cohort")
    request.add_argument("--cohort", type=Path, required=True)
    request.add_argument("--accessions", type=Path, required=True)
    request.add_argument("--mapping", type=Path, required=True)
    request.add_argument("--exclusions", type=Path, required=True)
    request.add_argument("--manifest", type=Path, required=True)

    stage = subparsers.add_parser("stage", help="Validate and stage a rehydrated NCBI package")
    stage.add_argument("--requests", type=Path, required=True)
    stage.add_argument("--package-dir", type=Path, required=True)
    stage.add_argument("--stage-dir", type=Path, required=True)
    stage.add_argument("--table", type=Path, required=True)
    stage.add_argument("--exclusions", type=Path, required=True)
    stage.add_argument("--manifest", type=Path, required=True)

    final = subparsers.add_parser("finalize", help="Consolidate assembly analysis outputs")
    final.add_argument("--stage-table", type=Path, required=True)
    final.add_argument("--quast-root", type=Path, required=True)
    final.add_argument("--mlst-root", type=Path, required=True)
    final.add_argument("--amrfinder-root", type=Path, required=True)
    final.add_argument("--tool-versions", type=Path, required=True)
    final.add_argument("--mlst-scheme-info", type=Path, required=True)
    final.add_argument("--amrfinder-database-metadata", type=Path, required=True)
    final.add_argument("--amrfinder-database-checksums", type=Path, required=True)
    final.add_argument("--output", type=Path, required=True)
    final.add_argument("--exclusions", type=Path, required=True)
    final.add_argument("--manifest", type=Path, required=True)
    final.add_argument("--min-total-length", type=int, required=True)
    final.add_argument("--max-total-length", type=int, required=True)
    final.add_argument("--max-contigs", type=int, required=True)
    final.add_argument("--min-n50", type=int, required=True)
    final.add_argument("--max-ns-per-100kb", type=float, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "request":
        payload = write_request_outputs(
            args.cohort, args.accessions, args.mapping, args.exclusions, args.manifest
        )
        print(
            f"Genome requests: {payload['requested_assemblies']} assemblies; "
            f"{payload['excluded_rows']} excluded cohort rows"
        )
    elif args.command == "stage":
        payload = write_stage_outputs(
            args.requests,
            args.package_dir,
            args.stage_dir,
            args.table,
            args.exclusions,
            args.manifest,
        )
        print(
            f"Genome staging: {payload['staged_assemblies']}/{payload['requested_assemblies']} staged"
        )
    else:
        payload = write_final_outputs(args)
        print(
            f"Genome analyses: {payload['eligible_assembly_count']}/{payload['assembly_count']} eligible; "
            f"status={payload['status']}"
        )


if __name__ == "__main__":
    main()
