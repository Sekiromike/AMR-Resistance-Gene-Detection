"""Acquire and audit BioSample antibiograms for a frozen NCBI AST snapshot.

The command reads the exact BioSample accessions from the hash-verified parent
snapshot.  Live mode archives raw EFetch XML in bounded batches, emits every
BioSample antibiogram row without filling missing metadata, and cross-checks
the rows against the frozen native AST records.  A composite match is audit
evidence only: BioSample antibiograms do not carry the Pathogen Detection AST
record identifier, so ambiguous matches are retained rather than selected.

``--dry-run`` is completely offline.  It verifies the parent snapshot and
writes a request plan and planned manifest without contacting NCBI.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import http.client
import json
import os
import platform
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import acquire_ncbi_isolate_enrichment as parent_snapshot


EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BIOSAMPLE_DOCUMENTATION = "https://www.ncbi.nlm.nih.gov/biosample/docs/"
ANTIBIOGRAM_DOCUMENTATION = (
    "https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/"
)
EUTILITIES_DOCUMENTATION = "https://www.ncbi.nlm.nih.gov/books/NBK25497/"
TOOL_NAME = "amr_ecoli_biosample_overlay"

PLAN_NAME = "ncbi_biosample_antibiogram_requests.json"
OVERLAY_NAME = "ncbi_biosample_antibiograms.csv"
RECONCILIATION_NAME = "ncbi_biosample_ast_reconciliation.csv"
MANIFEST_NAME = "ncbi_biosample_antibiogram_manifest.json"
RAW_DIRECTORY_NAME = "raw_xml"

DEFAULT_BATCH_SIZE = 200
MAX_BATCH_SIZE = 200
DEFAULT_REQUESTS_PER_SECOND = 2.5
DEFAULT_RETRIES = 4

HEADER_ALIASES = {
    "antibiotic": "antibiotic_raw",
    "resistance phenotype": "resistance_phenotype_raw",
    "measurement sign": "measurement_sign_raw",
    "measurement": "measurement_raw",
    "measurement unit": "measurement_unit_raw",
    "measurement units": "measurement_unit_raw",
    "laboratory typing method": "laboratory_typing_method_raw",
    "laboratory typing platform": "laboratory_typing_platform_raw",
    "vendor": "vendor_raw",
    "laboratory typing method version or reagent": (
        "laboratory_typing_method_version_or_reagent_raw"
    ),
    "testing standard": "testing_standard_raw",
    "testing standard version": "testing_standard_version_raw",
    "ast testing date": "ast_testing_date_raw",
    "testing date": "ast_testing_date_raw",
}

OVERLAY_COLUMNS = (
    "biosample_accession",
    "biosample_record_submission_date_raw",
    "biosample_record_last_update_raw",
    "antibiogram_table_index",
    "antibiogram_row_index",
    "antibiogram_row_id",
    "antibiotic_raw",
    "resistance_phenotype_raw",
    "measurement_sign_raw",
    "measurement_raw",
    "measurement_unit_raw",
    "laboratory_typing_method_raw",
    "laboratory_typing_platform_raw",
    "vendor_raw",
    "laboratory_typing_method_version_or_reagent_raw",
    "testing_standard_raw",
    "testing_standard_version_raw",
    "ast_testing_date_raw",
    "header_cells_json",
    "row_cells_json",
    "unmapped_cells_json",
    "biosample_record_xml_sha256",
    "raw_batch_file",
)

RECONCILIATION_COLUMNS = (
    "source_ast_record_id",
    "source_record_checksum",
    "target_acc",
    "biosample_accession",
    "antibiotic_raw",
    "phenotype_raw",
    "measurement_type",
    "measurement_sign_raw",
    "measurement_raw",
    "measurement_unit",
    "submitted_standard_raw",
    "ast_platform_raw",
    "ast_vendor_raw",
    "ast_reagent_raw",
    "core_match_count",
    "core_match_row_ids_json",
    "reconciliation_status",
    "unique_core_match_row_id",
    "standard_exact_when_present",
    "platform_exact_when_present",
    "vendor_exact_when_present",
    "reagent_exact_when_present",
)

MISSING_TOKENS = {"", "missing", "not collected", "not applicable", "na", "n/a"}


class BioSampleOverlayError(RuntimeError):
    """Raised when acquisition or validation cannot proceed unambiguously."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(Path.cwd().resolve()).as_posix()
    except ValueError:
        return path.name


def normalized_text(value: str) -> str:
    return " ".join(value.strip().casefold().split())


def is_present(value: str) -> bool:
    return normalized_text(value) not in MISSING_TOKENS


def batches(values: Sequence[str], batch_size: int) -> list[list[str]]:
    if batch_size < 1 or batch_size > MAX_BATCH_SIZE:
        raise BioSampleOverlayError(
            f"Batch size must be between 1 and {MAX_BATCH_SIZE}."
        )
    return [list(values[index : index + batch_size]) for index in range(0, len(values), batch_size)]


def accession_list_sha256(accessions: Iterable[str]) -> str:
    payload = "".join(f"{value}\n" for value in accessions).encode("utf-8")
    return sha256_bytes(payload)


def build_request_plan(accessions: Sequence[str], batch_size: int) -> dict[str, Any]:
    ordered = sorted(accessions)
    if not ordered or len(ordered) != len(set(ordered)):
        raise BioSampleOverlayError("BioSample accessions must be nonempty and unique.")
    invalid = [value for value in ordered if not parent_snapshot.BIOSAMPLE_RE.fullmatch(value)]
    if invalid:
        raise BioSampleOverlayError(f"Invalid BioSample accession: {invalid[0]!r}")
    records = []
    for index, values in enumerate(batches(ordered, batch_size), start=1):
        records.append(
            {
                "batch_index": index,
                "raw_file": f"{RAW_DIRECTORY_NAME}/batch_{index:05d}.xml",
                "accession_count": len(values),
                "accessions_sha256": accession_list_sha256(values),
                "accessions": values,
            }
        )
    return {
        "schema_version": "1.0.0",
        "endpoint": EFETCH_URL,
        "http_method": "POST",
        "parameters_without_private_values": {
            "db": "biosample",
            "retmode": "xml",
            "tool": TOOL_NAME,
            "email_supplied_at_runtime": True,
            "api_key_supplied_from_environment_when_available": True,
        },
        "batch_size": batch_size,
        "batch_count": len(records),
        "accession_count": len(ordered),
        "accession_list_serialization": "UTF-8, sorted accession plus LF per record",
        "accession_list_sha256": accession_list_sha256(ordered),
        "batches": records,
    }


def _element_text(element: ET.Element | None) -> str:
    if element is None:
        return ""
    return "".join(element.itertext()).strip()


def _record_accession(record: ET.Element) -> str:
    candidates = []
    for identifier in record.findall("./Ids/Id"):
        value = _element_text(identifier)
        database = normalized_text(
            identifier.attrib.get("db", "") or identifier.attrib.get("db_label", "")
        )
        if database == "biosample" or parent_snapshot.BIOSAMPLE_RE.fullmatch(value):
            candidates.append(value)
    unique = sorted(set(candidates))
    if len(unique) != 1:
        raise BioSampleOverlayError(
            f"BioSample XML record has {len(unique)} primary accession candidates."
        )
    return unique[0]


def _row_id(accession: str, table_index: int, row_index: int, cells: Sequence[str]) -> str:
    payload = json.dumps(
        [accession, table_index, row_index, list(cells)],
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return f"{accession}:antibiogram:{sha256_bytes(payload)[:20]}"


def parse_biosample_xml(
    payload: bytes,
    *,
    expected_accessions: Sequence[str],
    raw_batch_file: str,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    try:
        root = ET.fromstring(payload)
    except ET.ParseError as exc:
        raise BioSampleOverlayError("NCBI EFetch response is not valid XML.") from exc
    if root.tag.rsplit("}", 1)[-1] != "BioSampleSet":
        raise BioSampleOverlayError(f"Unexpected EFetch root element: {root.tag!r}")
    records = list(root.findall("./BioSample"))
    by_accession: dict[str, ET.Element] = {}
    duplicates: list[str] = []
    for record in records:
        accession = _record_accession(record)
        if accession in by_accession:
            duplicates.append(accession)
        by_accession[accession] = record
    expected = set(expected_accessions)
    observed = set(by_accession)
    if duplicates or expected != observed:
        raise BioSampleOverlayError(
            "BioSample batch identity mismatch: "
            f"missing={sorted(expected - observed)[:5]}, "
            f"extra={sorted(observed - expected)[:5]}, "
            f"duplicates={sorted(set(duplicates))[:5]}."
        )

    overlay_rows: list[dict[str, str]] = []
    accessions_with_antibiograms = 0
    unknown_headers: Counter[str] = Counter()
    for accession in sorted(by_accession):
        record = by_accession[accession]
        record_sha = sha256_bytes(ET.tostring(record, encoding="utf-8"))
        record_rows = 0
        tables = [
            table
            for table in record.findall("./Description/Comment/Table")
            if normalized_text(_element_text(table.find("./Caption"))) == "antibiogram"
        ]
        for table_index, table in enumerate(tables, start=1):
            headers = [_element_text(cell) for cell in table.findall("./Header/Cell")]
            if not headers:
                raise BioSampleOverlayError(
                    f"BioSample {accession} has an antibiogram table without headers."
                )
            mapped_headers: list[str | None] = []
            for header in headers:
                mapped = HEADER_ALIASES.get(normalized_text(header))
                mapped_headers.append(mapped)
                if mapped is None:
                    unknown_headers[header] += 1
            for row_index, row in enumerate(table.findall("./Body/Row"), start=1):
                cells = [_element_text(cell) for cell in row.findall("./Cell")]
                if len(cells) != len(headers):
                    raise BioSampleOverlayError(
                        f"BioSample {accession} antibiogram row has {len(cells)} cells "
                        f"for {len(headers)} headers."
                    )
                output = {column: "" for column in OVERLAY_COLUMNS}
                output.update(
                    {
                        "biosample_accession": accession,
                        "biosample_record_submission_date_raw": record.attrib.get(
                            "submission_date", ""
                        ),
                        "biosample_record_last_update_raw": record.attrib.get(
                            "last_update", ""
                        ),
                        "antibiogram_table_index": str(table_index),
                        "antibiogram_row_index": str(row_index),
                        "antibiogram_row_id": _row_id(
                            accession, table_index, row_index, cells
                        ),
                        "header_cells_json": json.dumps(headers, ensure_ascii=False),
                        "row_cells_json": json.dumps(cells, ensure_ascii=False),
                        "biosample_record_xml_sha256": record_sha,
                        "raw_batch_file": raw_batch_file,
                    }
                )
                unmapped: dict[str, str] = {}
                for header, mapped, value in zip(headers, mapped_headers, cells):
                    if mapped is None:
                        unmapped[header] = value
                    elif output[mapped] and output[mapped] != value:
                        raise BioSampleOverlayError(
                            f"BioSample {accession} has duplicate conflicting header {header!r}."
                        )
                    else:
                        output[mapped] = value
                output["unmapped_cells_json"] = json.dumps(
                    unmapped, ensure_ascii=False, sort_keys=True
                )
                overlay_rows.append(output)
                record_rows += 1
        if record_rows:
            accessions_with_antibiograms += 1
    return overlay_rows, {
        "biosample_records": len(records),
        "biosamples_with_antibiograms": accessions_with_antibiograms,
        "antibiogram_rows": len(overlay_rows),
        "unknown_headers": dict(sorted(unknown_headers.items())),
    }


def read_parent_ast(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        required = {
            "id",
            "checksum",
            "target_acc",
            "biosample_acc",
            "antibiotic",
            "phenotype",
            "measurement_sign",
            "mic",
            "disk_diffusion",
            "standard",
            "platform",
            "vendor",
            "reagent",
        }
        missing = sorted(required.difference(fields))
        if missing:
            raise BioSampleOverlayError(f"Parent AST CSV lacks columns: {missing}")
        rows = list(reader)
    if not rows:
        raise BioSampleOverlayError("Parent AST CSV has no records.")
    return rows


def _decimal_equal(left: str, right: str) -> bool:
    try:
        return Decimal(left.strip()) == Decimal(right.strip())
    except (InvalidOperation, ValueError):
        return False


def _parent_measurement(row: Mapping[str, str]) -> tuple[str, str, str]:
    mic = row.get("mic", "").strip()
    zone = row.get("disk_diffusion", "").strip()
    if bool(mic) == bool(zone):
        return "", "", ""
    if mic:
        return "MIC", mic, "mg/L"
    return "disk diffusion", zone, "mm"


def _same_when_parent_present(parent: str, source: str) -> str:
    if not is_present(parent):
        return "not_asserted_parent_missing"
    return str(normalized_text(parent) == normalized_text(source))


def reconcile_ast_rows(
    parent_rows: Sequence[Mapping[str, str]],
    overlay_rows: Sequence[Mapping[str, str]],
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    by_biosample: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in overlay_rows:
        by_biosample[row["biosample_accession"]].append(row)

    output_rows: list[dict[str, str]] = []
    reverse_matches: Counter[str] = Counter()
    statuses: Counter[str] = Counter()
    for parent in parent_rows:
        method, measurement, unit = _parent_measurement(parent)
        matches = []
        for candidate in by_biosample.get(parent["biosample_acc"], []):
            if normalized_text(parent["antibiotic"]) != normalized_text(
                candidate["antibiotic_raw"]
            ):
                continue
            if normalized_text(parent["phenotype"]) != normalized_text(
                candidate["resistance_phenotype_raw"]
            ):
                continue
            if parent["measurement_sign"].strip() != candidate[
                "measurement_sign_raw"
            ].strip():
                continue
            if not _decimal_equal(measurement, candidate["measurement_raw"]):
                continue
            if normalized_text(method) != normalized_text(
                candidate["laboratory_typing_method_raw"]
            ):
                continue
            if normalized_text(unit) != normalized_text(candidate["measurement_unit_raw"]):
                continue
            matches.append(candidate)

        match_ids = sorted(candidate["antibiogram_row_id"] for candidate in matches)
        for match_id in match_ids:
            reverse_matches[match_id] += 1
        if len(matches) == 1:
            status = "unique_core_match"
            unique = matches[0]
        elif matches:
            status = "ambiguous_core_match"
            unique = None
        else:
            status = "no_core_match"
            unique = None
        statuses[status] += 1
        output = {
            "source_ast_record_id": parent["id"],
            "source_record_checksum": parent["checksum"],
            "target_acc": parent["target_acc"],
            "biosample_accession": parent["biosample_acc"],
            "antibiotic_raw": parent["antibiotic"],
            "phenotype_raw": parent["phenotype"],
            "measurement_type": method,
            "measurement_sign_raw": parent["measurement_sign"],
            "measurement_raw": measurement,
            "measurement_unit": unit,
            "submitted_standard_raw": parent["standard"],
            "ast_platform_raw": parent["platform"],
            "ast_vendor_raw": parent["vendor"],
            "ast_reagent_raw": parent["reagent"],
            "core_match_count": str(len(matches)),
            "core_match_row_ids_json": json.dumps(match_ids),
            "reconciliation_status": status,
            "unique_core_match_row_id": (
                unique["antibiogram_row_id"] if unique is not None else ""
            ),
            "standard_exact_when_present": (
                _same_when_parent_present(parent["standard"], unique["testing_standard_raw"])
                if unique is not None
                else "not_assessed_nonunique"
            ),
            "platform_exact_when_present": (
                _same_when_parent_present(
                    parent["platform"], unique["laboratory_typing_platform_raw"]
                )
                if unique is not None
                else "not_assessed_nonunique"
            ),
            "vendor_exact_when_present": (
                _same_when_parent_present(parent["vendor"], unique["vendor_raw"])
                if unique is not None
                else "not_assessed_nonunique"
            ),
            "reagent_exact_when_present": (
                _same_when_parent_present(
                    parent["reagent"],
                    unique["laboratory_typing_method_version_or_reagent_raw"],
                )
                if unique is not None
                else "not_assessed_nonunique"
            ),
        }
        output_rows.append(output)

    overlay_ids = {row["antibiogram_row_id"] for row in overlay_rows}
    return output_rows, {
        "parent_ast_rows": len(parent_rows),
        "status_counts": dict(sorted(statuses.items())),
        "overlay_rows_with_no_parent_core_match": sum(
            1 for row_id in overlay_ids if reverse_matches[row_id] == 0
        ),
        "overlay_rows_with_one_parent_core_match": sum(
            1 for row_id in overlay_ids if reverse_matches[row_id] == 1
        ),
        "overlay_rows_with_multiple_parent_core_matches": sum(
            1 for row_id in overlay_ids if reverse_matches[row_id] > 1
        ),
    }


def field_coverage(rows: Sequence[Mapping[str, str]]) -> dict[str, dict[str, int]]:
    fields = (
        "antibiotic_raw",
        "resistance_phenotype_raw",
        "measurement_sign_raw",
        "measurement_raw",
        "measurement_unit_raw",
        "laboratory_typing_method_raw",
        "laboratory_typing_platform_raw",
        "vendor_raw",
        "laboratory_typing_method_version_or_reagent_raw",
        "testing_standard_raw",
        "testing_standard_version_raw",
        "ast_testing_date_raw",
    )
    coverage = {}
    for field in fields:
        present = [row[field] for row in rows if is_present(row[field])]
        coverage[field] = {
            "present": len(present),
            "missing": len(rows) - len(present),
            "unique_nonmissing": len(set(present)),
        }
    return coverage


def _write_csv(path: Path, columns: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _check_output_dir(path: Path) -> None:
    if path.exists():
        raise BioSampleOverlayError(f"Output path already exists: {path}")
    if path.resolve() in {Path.cwd().resolve(), Path.cwd().resolve().parent}:
        raise BioSampleOverlayError("Output directory is too broad.")


def _install_directory(staging: Path, output_dir: Path) -> None:
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging.replace(output_dir)


def base_manifest(
    *, context: Mapping[str, Any], plan: Mapping[str, Any], generated_at: str
) -> dict[str, Any]:
    return {
        "schema_version": "1.0.0",
        "artifact_type": "ncbi_biosample_antibiogram_overlay",
        "status": "planned",
        "generated_at_utc": generated_at,
        "source": {
            "provider": "NCBI BioSample",
            "endpoint": EFETCH_URL,
            "documentation": BIOSAMPLE_DOCUMENTATION,
            "antibiogram_documentation": ANTIBIOGRAM_DOCUMENTATION,
            "eutilities_policy": EUTILITIES_DOCUMENTATION,
            "temporal_semantics": (
                "live mutable records retrieved at recorded UTC times; BioSample EFetch "
                "does not provide historical time-travel"
            ),
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
            "basis": "sorted unique BioSample accessions from the hash-verified parent isolate CSV",
            "biosample_count": len(context["identities"]),
            "biosample_list_sha256": accession_list_sha256(
                sorted(identity["biosample_acc"] for identity in context["identities"].values())
            ),
            "parent_ast_rows": context["rows"]["ast"],
            "mutable_ast_table_reselected": False,
        },
        "request_plan": {
            "file": PLAN_NAME,
            "sha256": "pending",
            "batch_size": plan["batch_size"],
            "batch_count": plan["batch_count"],
            "http_method": "POST",
            "private_values_persisted": False,
        },
        "retrieval": None,
        "coverage": None,
        "reconciliation": None,
        "limitations": [
            "BioSample records are mutable and EFetch has no historical time-travel; raw XML hashes freeze only this retrieval.",
            "BioSample antibiograms do not expose the NCBI Pathogen Detection native AST record identifier.",
            "Composite matches are retained as ambiguous when more than one source row qualifies.",
            "A BioSample record submission or update timestamp is not an AST testing date.",
            "A collection date is not an AST testing date.",
            "Testing-standard name without an explicit version does not establish breakpoint version.",
            "Platform, vendor, and reagent do not establish laboratory site or clinical indication.",
            "No overlay field is authorized as a model feature by this acquisition.",
        ],
        "evidence_tier": "SOURCE_COHORT_ONLY",
    }


def write_plan(
    output_dir: Path, plan: dict[str, Any], manifest: dict[str, Any]
) -> dict[str, Any]:
    _check_output_dir(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="amr-biosample-plan-", dir=output_dir.parent))
    try:
        _write_json(staging / PLAN_NAME, plan)
        manifest["request_plan"]["sha256"] = sha256_file(staging / PLAN_NAME)
        _write_json(staging / MANIFEST_NAME, manifest)
        _install_directory(staging, output_dir)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return manifest


def fetch_batch(
    accessions: Sequence[str], *, email: str, api_key: str, timeout: float
) -> bytes:
    parameters = {
        "db": "biosample",
        "id": ",".join(accessions),
        "retmode": "xml",
        "tool": TOOL_NAME,
        "email": email,
    }
    if api_key:
        parameters["api_key"] = api_key
    request = urllib.request.Request(
        EFETCH_URL,
        data=urllib.parse.urlencode(parameters).encode("ascii"),
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "User-Agent": f"{TOOL_NAME}/1.0",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = response.read()
    except (
        urllib.error.URLError,
        http.client.HTTPException,
        ConnectionError,
        TimeoutError,
    ) as exc:
        raise BioSampleOverlayError(f"NCBI EFetch request failed: {exc}") from exc
    if not payload:
        raise BioSampleOverlayError("NCBI EFetch returned an empty response.")
    return payload


def acquire(
    *,
    output_dir: Path,
    context: Mapping[str, Any],
    plan: dict[str, Any],
    manifest: dict[str, Any],
    email: str,
    api_key: str,
    requests_per_second: float,
    retries: int,
    timeout: float,
) -> dict[str, Any]:
    _check_output_dir(output_dir)
    if not email.strip() or "@" not in email:
        raise BioSampleOverlayError(
            "Live retrieval requires a contact email via --email or NCBI_EMAIL."
        )
    maximum_rate = 10.0 if api_key else 3.0
    if not (0 < requests_per_second <= maximum_rate):
        raise BioSampleOverlayError(
            f"Request rate must be >0 and <= {maximum_rate:g}/s for this credential mode."
        )
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="amr-biosample-live-", dir=output_dir.parent))
    raw_dir = staging / RAW_DIRECTORY_NAME
    raw_dir.mkdir()
    started = utc_now()
    all_rows: list[dict[str, str]] = []
    batch_records = []
    try:
        _write_json(staging / PLAN_NAME, plan)
        manifest["request_plan"]["sha256"] = sha256_file(staging / PLAN_NAME)
        for record in plan["batches"]:
            payload = b""
            last_error: Exception | None = None
            for attempt in range(1, retries + 1):
                try:
                    payload = fetch_batch(
                        record["accessions"], email=email, api_key=api_key, timeout=timeout
                    )
                    last_error = None
                    break
                except BioSampleOverlayError as exc:
                    last_error = exc
                    if attempt < retries:
                        time.sleep(min(30.0, 2.0 ** (attempt - 1)))
            if last_error is not None:
                raise last_error
            raw_path = staging / record["raw_file"]
            raw_path.write_bytes(payload)
            rows, validation = parse_biosample_xml(
                payload,
                expected_accessions=record["accessions"],
                raw_batch_file=record["raw_file"],
            )
            all_rows.extend(rows)
            batch_records.append(
                {
                    "batch_index": record["batch_index"],
                    "file": record["raw_file"],
                    "bytes": len(payload),
                    "sha256": sha256_bytes(payload),
                    **validation,
                }
            )
            time.sleep(1.0 / requests_per_second)

        parent_snapshot.assert_parent_unchanged(context)
        parent_rows = read_parent_ast(context["paths"]["ast"])
        reconciliation_rows, reconciliation = reconcile_ast_rows(parent_rows, all_rows)
        _write_csv(staging / OVERLAY_NAME, OVERLAY_COLUMNS, all_rows)
        _write_csv(
            staging / RECONCILIATION_NAME,
            RECONCILIATION_COLUMNS,
            reconciliation_rows,
        )
        manifest["status"] = "complete"
        manifest["coverage"] = {
            "biosample_records": len(context["identities"]),
            "biosamples_with_antibiograms": len(
                {row["biosample_accession"] for row in all_rows}
            ),
            "biosamples_without_antibiograms": len(context["identities"])
            - len({row["biosample_accession"] for row in all_rows}),
            "antibiogram_rows": len(all_rows),
            "fields": field_coverage(all_rows),
        }
        manifest["reconciliation"] = reconciliation
        manifest["retrieval"] = {
            "started_at_utc": started,
            "completed_at_utc": utc_now(),
            "client": "Python urllib.request",
            "contact_email_supplied": True,
            "contact_email_persisted": False,
            "api_key_supplied": bool(api_key),
            "api_key_persisted": False,
            "requests_per_second": requests_per_second,
            "raw_batches": batch_records,
            "outputs": {
                "overlay": {
                    "file": OVERLAY_NAME,
                    "rows": len(all_rows),
                    "bytes": (staging / OVERLAY_NAME).stat().st_size,
                    "sha256": sha256_file(staging / OVERLAY_NAME),
                },
                "reconciliation": {
                    "file": RECONCILIATION_NAME,
                    "rows": len(reconciliation_rows),
                    "bytes": (staging / RECONCILIATION_NAME).stat().st_size,
                    "sha256": sha256_file(staging / RECONCILIATION_NAME),
                },
            },
        }
        _write_json(staging / MANIFEST_NAME, manifest)
        _install_directory(staging, output_dir)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
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
    parser.add_argument("--parent-amendment", type=Path)
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""))
    parser.add_argument("--api-key-env", default="NCBI_API_KEY")
    parser.add_argument(
        "--requests-per-second", type=float, default=DEFAULT_REQUESTS_PER_SECOND
    )
    parser.add_argument("--retries", type=int, default=DEFAULT_RETRIES)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        amendment = args.parent_amendment
        inferred = args.parent_manifest.parent / parent_snapshot.PARENT_AMENDMENT_NAME
        if amendment is None and inferred.is_file():
            amendment = inferred
        context = parent_snapshot.load_parent_context(
            isolates_path=args.parent_isolates,
            ast_path=args.parent_ast,
            acquisition_manifest_path=args.parent_manifest,
            amendment_path=amendment,
        )
        accessions = sorted(
            identity["biosample_acc"] for identity in context["identities"].values()
        )
        plan = build_request_plan(accessions, args.batch_size)
        manifest = base_manifest(context=context, plan=plan, generated_at=utc_now())
        if args.dry_run:
            write_plan(args.output_dir, plan, manifest)
            print(
                f"Planned {plan['accession_count']} BioSamples in "
                f"{plan['batch_count']} immutable POST batches at {args.output_dir}."
            )
            return 0
        api_key = os.environ.get(args.api_key_env, "") if args.api_key_env else ""
        completed = acquire(
            output_dir=args.output_dir,
            context=context,
            plan=plan,
            manifest=manifest,
            email=args.email,
            api_key=api_key,
            requests_per_second=args.requests_per_second,
            retries=args.retries,
            timeout=args.timeout,
        )
        print(
            f"Acquired {completed['coverage']['biosample_records']} BioSamples and "
            f"{completed['coverage']['antibiogram_rows']} antibiogram rows at "
            f"{args.output_dir}."
        )
        return 0
    except (
        BioSampleOverlayError,
        parent_snapshot.EnrichmentError,
        OSError,
        ValueError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
