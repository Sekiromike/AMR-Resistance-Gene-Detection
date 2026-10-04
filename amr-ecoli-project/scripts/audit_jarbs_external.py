"""Audit the frozen JARBS supplement without interpreting missing fields.

The workbook reader intentionally uses only the Python standard library.  It
is small enough to keep the source-to-cell mapping reviewable and does not
modify the publisher workbook.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Iterable


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PACKAGE_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
CELL_REFERENCE_RE = re.compile(r"([A-Z]+)([0-9]+)")
RUN_ACCESSION_RE = re.compile(r"(?:SRR|ERR|DRR)[0-9]+")
MIC_RE = re.compile(r"^(<=|>=|<|>)?([0-9]+(?:\.[0-9]+)?)$")
TARGET_DRUG_COLUMNS = {
    "ceftriaxone": "CTRX",
    "ciprofloxacin": "CPFX",
    "gentamicin": "GM",
}


class JarbsAuditError(RuntimeError):
    """Raised when a source artifact cannot be interpreted unambiguously."""


def _column_index(cell_reference: str) -> int:
    match = CELL_REFERENCE_RE.fullmatch(cell_reference)
    if match is None:
        raise JarbsAuditError(f"Invalid XLSX cell reference: {cell_reference!r}")
    result = 0
    for character in match.group(1):
        result = result * 26 + ord(character) - ord("A") + 1
    return result - 1


def _shared_strings(archive: zipfile.ZipFile) -> list[str]:
    try:
        payload = archive.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    root = ET.fromstring(payload)
    return [
        "".join(node.text or "" for node in item.iter(f"{{{MAIN_NS}}}t"))
        for item in root.findall(f"{{{MAIN_NS}}}si")
    ]


def _sheet_paths(archive: zipfile.ZipFile) -> list[tuple[str, str]]:
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relations = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {
        relation.attrib["Id"]: relation.attrib["Target"]
        for relation in relations.findall(f"{{{PACKAGE_REL_NS}}}Relationship")
    }
    result = []
    for sheet in workbook.findall(f".//{{{MAIN_NS}}}sheet"):
        relation_id = sheet.attrib[f"{{{REL_NS}}}id"]
        target = targets[relation_id].replace("\\", "/").lstrip("/")
        if not target.startswith("xl/"):
            target = f"xl/{target}"
        result.append((sheet.attrib["name"], target))
    return result


def _cell_value(cell: ET.Element, shared_strings: list[str]) -> str:
    cell_type = cell.attrib.get("t", "")
    if cell_type == "inlineStr":
        return "".join(
            node.text or "" for node in cell.iter(f"{{{MAIN_NS}}}t")
        )
    value = cell.find(f"{{{MAIN_NS}}}v")
    raw = "" if value is None or value.text is None else value.text
    if cell_type == "s" and raw:
        try:
            return shared_strings[int(raw)]
        except (IndexError, ValueError) as exc:
            raise JarbsAuditError(f"Invalid shared-string index: {raw!r}") from exc
    if cell_type == "b":
        return "TRUE" if raw == "1" else "FALSE"
    return raw


def read_xlsx(path: Path) -> dict[str, list[list[str]]]:
    """Return visible cell values for every worksheet in workbook order."""
    if not path.is_file():
        raise JarbsAuditError(f"Workbook not found: {path}")
    if not zipfile.is_zipfile(path):
        raise JarbsAuditError(f"Workbook is not an OOXML ZIP container: {path}")
    with zipfile.ZipFile(path) as archive:
        shared_strings = _shared_strings(archive)
        sheets: dict[str, list[list[str]]] = {}
        for sheet_name, sheet_path in _sheet_paths(archive):
            root = ET.fromstring(archive.read(sheet_path))
            rows: list[list[str]] = []
            for row in root.findall(f".//{{{MAIN_NS}}}sheetData/{{{MAIN_NS}}}row"):
                indexed_values = {
                    _column_index(cell.attrib["r"]): _cell_value(
                        cell, shared_strings
                    )
                    for cell in row.findall(f"{{{MAIN_NS}}}c")
                }
                width = max(indexed_values, default=-1) + 1
                rows.append([indexed_values.get(index, "") for index in range(width)])
            sheets[sheet_name] = rows
    return sheets


def workbook_inventory(sheets: dict[str, list[list[str]]]) -> dict[str, Any]:
    result: dict[str, Any] = {"sheet_count": len(sheets), "sheets": []}
    for name, rows in sheets.items():
        width = max((len(row) for row in rows), default=0)
        result["sheets"].append(
            {
                "name": name,
                "row_count_including_header": len(rows),
                "maximum_column_count": width,
                "header": rows[0] if rows else [],
                "first_data_rows": rows[1:4],
            }
        )
    return result


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def _as_records(rows: list[list[str]]) -> list[dict[str, str]]:
    if not rows:
        raise JarbsAuditError("The supplement worksheet is empty.")
    header = rows[0]
    if not header or len(header) != len(set(header)):
        raise JarbsAuditError("The supplement header is empty or contains duplicates.")
    records = []
    for row in rows[1:]:
        padded = row + [""] * (len(header) - len(row))
        records.append(dict(zip(header, padded[: len(header)])))
    return records


def _mic_parts(raw_value: str) -> tuple[str, str]:
    match = MIC_RE.fullmatch(raw_value.strip())
    if match is None:
        raise JarbsAuditError(f"Unexpected target-drug MIC value: {raw_value!r}")
    comparator, number = match.groups()
    try:
        normalized = format(Decimal(number), "f")
    except InvalidOperation as exc:  # pragma: no cover - guarded by the regex
        raise JarbsAuditError(f"Invalid MIC numeric value: {raw_value!r}") from exc
    return comparator or "=", normalized


def audit_supplement(
    sheets: dict[str, list[list[str]]],
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    if list(sheets) != ["Dataset 6"]:
        raise JarbsAuditError(f"Unexpected workbook sheets: {list(sheets)!r}")
    records = _as_records(sheets["Dataset 6"])
    footnotes = [row for row in records if not row.get("isolate ID", "").strip()]
    isolate_rows = [row for row in records if row.get("isolate ID", "").strip()]
    required = {"isolate ID", "genome_species", *TARGET_DRUG_COLUMNS.values()}
    missing_columns = sorted(required.difference(records[0] if records else {}))
    if missing_columns:
        raise JarbsAuditError(f"Missing supplement columns: {missing_columns}")
    isolate_ids = [row["isolate ID"].strip() for row in isolate_rows]
    duplicates = sorted(
        value for value, count in Counter(isolate_ids).items() if count > 1
    )
    if duplicates:
        raise JarbsAuditError(f"Duplicate supplement isolate ID: {duplicates[0]}")

    species_counts = Counter(row["genome_species"].strip() for row in isolate_rows)
    ecoli_rows = [
        row for row in isolate_rows if row["genome_species"].strip() == "Escherichia coli"
    ]
    target_summary: dict[str, Any] = {}
    for drug, column in TARGET_DRUG_COLUMNS.items():
        present_values = [row[column].strip() for row in ecoli_rows if row[column].strip()]
        parsed = [_mic_parts(value) for value in present_values]
        target_summary[drug] = {
            "column": column,
            "present": len(present_values),
            "missing": len(ecoli_rows) - len(present_values),
            "censored": sum(comparator != "=" for comparator, _ in parsed),
            "comparator_counts": dict(sorted(Counter(x[0] for x in parsed).items())),
            "distinct_numeric_values": len({x[1] for x in parsed}),
        }
    all_three = sum(
        all(row[column].strip() for column in TARGET_DRUG_COLUMNS.values())
        for row in ecoli_rows
    )
    summary = {
        "worksheet": "Dataset 6",
        "worksheet_rows_including_header": len(sheets["Dataset 6"]),
        "isolate_rows": len(isolate_rows),
        "footnote_rows": len(footnotes),
        "footnote_text": [
            value
            for row in footnotes
            for value in row.values()
            if value.strip()
        ],
        "unique_isolate_ids": len(set(isolate_ids)),
        "species_counts": dict(sorted(species_counts.items())),
        "ecoli_isolates": len(ecoli_rows),
        "ecoli_all_three_target_mics_present": all_three,
        "ecoli_target_drugs": target_summary,
        "explicit_field_presence": {
            "isolate_id": True,
            "genome_qc": True,
            "species": True,
            "quantitative_mic": True,
            "mlst": True,
            "resistance_genes": True,
            "biosample_accession": False,
            "sra_run_accession": False,
            "collection_date": False,
            "ast_testing_date": False,
            "hospital_or_collection_site": False,
            "specimen_source": False,
            "clinical_indication": False,
            "ast_method": False,
            "breakpoint_standard_version": False,
        },
    }
    return isolate_rows, summary


def _external_identifier(element: ET.Element, namespace: str) -> str:
    values = [
        (node.text or "").strip()
        for node in element.findall("./IDENTIFIERS/EXTERNAL_ID")
        if node.attrib.get("namespace") == namespace and (node.text or "").strip()
    ]
    unique = sorted(set(values))
    if len(unique) != 1:
        raise JarbsAuditError(
            f"Expected one {namespace} external identifier, found {unique!r}."
        )
    return unique[0]


def parse_sra_xml(path: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    root = ET.parse(path).getroot()
    run_rows: list[dict[str, str]] = []
    for package in root.findall("./EXPERIMENT_PACKAGE"):
        experiment = package.find("./EXPERIMENT")
        sample = package.find("./SAMPLE")
        submission = package.find("./SUBMISSION")
        if experiment is None or sample is None or submission is None:
            raise JarbsAuditError("Incomplete SRA experiment package.")
        attributes = {
            (item.findtext("TAG") or "").strip(): (item.findtext("VALUE") or "").strip()
            for item in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")
        }
        isolate_candidates = {
            attributes.get("strain", "").strip(),
            attributes.get("sample_name", "").strip(),
        }
        isolate_candidates.discard("")
        if len(isolate_candidates) != 1:
            raise JarbsAuditError(
                f"Conflicting SRA isolate identifiers: {sorted(isolate_candidates)!r}"
            )
        isolate_id = next(iter(isolate_candidates))
        platform_element = experiment.find("./PLATFORM")
        platform_children = [] if platform_element is None else list(platform_element)
        if len(platform_children) != 1:
            raise JarbsAuditError(f"Unexpected platform structure for {isolate_id}.")
        platform = platform_children[0].tag
        instrument = (platform_children[0].findtext("INSTRUMENT_MODEL") or "").strip()
        layout_element = experiment.find("./DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT")
        layout_children = [] if layout_element is None else list(layout_element)
        layout = layout_children[0].tag if len(layout_children) == 1 else ""
        biosample = _external_identifier(sample, "BioSample")
        bioproject = _external_identifier(package.find("./STUDY") or ET.Element("x"), "BioProject")
        for run in package.findall("./RUN_SET/RUN"):
            sra_file = run.find("./SRAFiles/SRAFile")
            run_rows.append(
                {
                    "isolate_id": isolate_id,
                    "scientific_name": (sample.findtext("./SAMPLE_NAME/SCIENTIFIC_NAME") or "").strip(),
                    "taxid": (sample.findtext("./SAMPLE_NAME/TAXON_ID") or "").strip(),
                    "biosample_accession": biosample,
                    "sra_sample_accession": sample.attrib.get("accession", ""),
                    "experiment_accession": experiment.attrib.get("accession", ""),
                    "run_accession": run.attrib.get("accession", ""),
                    "bioproject_accession": bioproject,
                    "submission_accession": submission.attrib.get("accession", ""),
                    "platform": platform,
                    "instrument_model": instrument,
                    "library_layout": layout,
                    "collection_date_raw": attributes.get("collection_date", ""),
                    "geo_loc_name_raw": attributes.get("geo_loc_name", ""),
                    "isolation_source_raw": attributes.get("isolation_source", ""),
                    "run_published_raw": run.attrib.get("published", ""),
                    "run_total_bases": run.attrib.get("total_bases", ""),
                    "run_file_md5": "" if sra_file is None else sra_file.attrib.get("md5", ""),
                }
            )
    run_accessions = [row["run_accession"] for row in run_rows]
    if not run_rows or len(run_accessions) != len(set(run_accessions)):
        raise JarbsAuditError("SRA run accessions are empty or duplicated.")
    isolate_counts = Counter(row["isolate_id"] for row in run_rows)
    summary = {
        "run_rows": len(run_rows),
        "unique_run_accessions": len(set(run_accessions)),
        "unique_isolate_ids": len(isolate_counts),
        "runs_per_isolate_counts": dict(sorted(Counter(isolate_counts.values()).items())),
        "platform_counts": dict(sorted(Counter(row["platform"] for row in run_rows).items())),
        "biosample_count": len({row["biosample_accession"] for row in run_rows}),
        "bioproject_counts": dict(sorted(Counter(row["bioproject_accession"] for row in run_rows).items())),
        "collection_date_present_runs": sum(bool(row["collection_date_raw"]) for row in run_rows),
        "geo_location_present_runs": sum(bool(row["geo_loc_name_raw"]) for row in run_rows),
        "isolation_source_present_runs": sum(bool(row["isolation_source_raw"]) for row in run_rows),
    }
    return run_rows, summary


def build_isolate_sequence_manifest(
    supplement_rows: list[dict[str, str]], run_rows: list[dict[str, str]]
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    by_isolate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in run_rows:
        by_isolate[row["isolate_id"]].append(row)
    supplement_ids = {row["isolate ID"] for row in supplement_rows}
    sra_ids = set(by_isolate)
    manifest = []
    for row in supplement_rows:
        isolate_id = row["isolate ID"]
        runs = sorted(by_isolate.get(isolate_id, []), key=lambda item: item["run_accession"])
        short_reads = [item for item in runs if item["platform"] == "ILLUMINA"]
        long_reads = [item for item in runs if item["platform"] != "ILLUMINA"]
        identity_fields = [
            "biosample_accession",
            "sra_sample_accession",
            "collection_date_raw",
            "geo_loc_name_raw",
            "isolation_source_raw",
        ]
        identities = {
            field: sorted({item[field] for item in runs if item[field]})
            for field in identity_fields
        }
        conflicts = {key: value for key, value in identities.items() if len(value) > 1}
        status = "linked_one_short_read_run"
        if len(short_reads) != 1:
            status = f"short_read_run_count_{len(short_reads)}"
        if conflicts:
            status = "conflicting_sra_sample_metadata"
        manifest.append(
            {
                "isolate_id": isolate_id,
                "supplement_species": row["genome_species"],
                "biosample_accession": identities["biosample_accession"][0] if len(identities["biosample_accession"]) == 1 else "",
                "sra_sample_accession": identities["sra_sample_accession"][0] if len(identities["sra_sample_accession"]) == 1 else "",
                "collection_date_raw": identities["collection_date_raw"][0] if len(identities["collection_date_raw"]) == 1 else "",
                "geo_loc_name_raw": identities["geo_loc_name_raw"][0] if len(identities["geo_loc_name_raw"]) == 1 else "",
                "isolation_source_raw": identities["isolation_source_raw"][0] if len(identities["isolation_source_raw"]) == 1 else "",
                "primary_short_read_run": short_reads[0]["run_accession"] if len(short_reads) == 1 else "",
                "primary_short_read_experiment": short_reads[0]["experiment_accession"] if len(short_reads) == 1 else "",
                "primary_short_read_md5": short_reads[0]["run_file_md5"] if len(short_reads) == 1 else "",
                "long_read_runs_json": json.dumps([item["run_accession"] for item in long_reads], separators=(",", ":")),
                "all_runs_json": json.dumps([item["run_accession"] for item in runs], separators=(",", ":")),
                "sequence_link_status": status,
            }
        )
    status_counts = Counter(row["sequence_link_status"] for row in manifest)
    summary = {
        "supplement_isolates_missing_from_sra": sorted(supplement_ids - sra_ids),
        "sra_isolates_missing_from_supplement": sorted(sra_ids - supplement_ids),
        "link_status_counts": dict(sorted(status_counts.items())),
        "ecoli_link_status_counts": dict(
            sorted(
                Counter(
                    row["sequence_link_status"]
                    for row in manifest
                    if row["supplement_species"] == "Escherichia coli"
                ).items()
            )
        ),
        "linked_isolates": len(manifest),
    }
    return manifest, summary


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_assembly_report(path: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    reports = payload.get("reports")
    if not isinstance(reports, list):
        raise JarbsAuditError("NCBI Datasets assembly report has no reports list.")
    rows = []
    for report in reports:
        organism = report.get("organism") or {}
        assembly_info = report.get("assembly_info") or {}
        infraspecific = organism.get("infraspecific_names") or {}
        rows.append(
            {
                "assembly_accession": str(report.get("accession") or ""),
                "paired_assembly_accession": str(report.get("paired_accession") or ""),
                "source_database": str(report.get("source_database") or ""),
                "organism_name": str(organism.get("organism_name") or ""),
                "taxid": str(organism.get("tax_id") or ""),
                "strain_raw": str(infraspecific.get("strain") or ""),
                "assembly_level": str(assembly_info.get("assembly_level") or ""),
                "assembly_status": str(assembly_info.get("assembly_status") or ""),
                "release_date_raw": str(assembly_info.get("release_date") or ""),
                "bioproject_accession": str(assembly_info.get("bioproject_accession") or ""),
                "sequencing_tech_raw": str(assembly_info.get("sequencing_tech") or ""),
            }
        )
    accessions = [row["assembly_accession"] for row in rows]
    if not rows or any(not value for value in accessions) or len(accessions) != len(set(accessions)):
        raise JarbsAuditError("Registered assembly accessions are empty or duplicated.")
    summary = {
        "registered_assembly_accessions": len(rows),
        "registered_assembly_strains": len({row["strain_raw"] for row in rows if row["strain_raw"]}),
        "source_database_counts": dict(sorted(Counter(row["source_database"] for row in rows).items())),
        "organism_counts": dict(sorted(Counter(row["organism_name"] for row in rows).items())),
        "ecoli_assembly_accessions": sum(row["organism_name"] == "Escherichia coli" for row in rows),
        "coverage_boundary": "registered_complete_assemblies_only_not_all_supplement_draft_genomes",
    }
    return rows, summary


def audit_development_overlap(
    sequence_manifest: list[dict[str, str]],
    assembly_rows: list[dict[str, str]],
    development_isolates_path: Path,
    development_enrichment_path: Path,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    isolates = _read_csv(development_isolates_path)
    enrichment = _read_csv(development_enrichment_path)
    targets_by_biosample: dict[str, set[str]] = defaultdict(set)
    targets_by_run: dict[str, set[str]] = defaultdict(set)
    targets_by_assembly: dict[str, set[str]] = defaultdict(set)
    targets_by_isolate_id: dict[str, set[str]] = defaultdict(set)
    for row in isolates:
        target = row.get("target_acc", "")
        if row.get("biosample_acc", ""):
            targets_by_biosample[row["biosample_acc"]].add(target)
        if row.get("strain", ""):
            targets_by_isolate_id[row["strain"]].add(target)
        if row.get("asm_acc", ""):
            targets_by_assembly[row["asm_acc"]].add(target)
    for row in enrichment:
        target = row.get("target_acc", "")
        for run in RUN_ACCESSION_RE.findall(row.get("sra_run_accessions_raw", "")):
            targets_by_run[run].add(target)
        raw_alternatives = row.get("alternative_isolate_identifiers_json", "")
        if raw_alternatives:
            try:
                alternatives = json.loads(raw_alternatives)
            except json.JSONDecodeError as exc:
                raise JarbsAuditError("Invalid alternative-isolate JSON.") from exc
            for value in alternatives:
                if value:
                    targets_by_isolate_id[str(value)].add(target)

    external_by_biosample: dict[str, set[str]] = defaultdict(set)
    external_by_run: dict[str, set[str]] = defaultdict(set)
    external_by_isolate_id: dict[str, set[str]] = defaultdict(set)
    external_by_assembly: dict[str, set[str]] = defaultdict(set)
    for row in sequence_manifest:
        isolate_id = row["isolate_id"]
        external_by_isolate_id[isolate_id].add(isolate_id)
        if row["biosample_accession"]:
            external_by_biosample[row["biosample_accession"]].add(isolate_id)
        for run in json.loads(row["all_runs_json"]):
            external_by_run[run].add(isolate_id)
    for row in assembly_rows:
        external_by_assembly[row["assembly_accession"]].add(row["strain_raw"])

    overlaps = []
    for identifier_type, development, external in (
        ("biosample_accession", targets_by_biosample, external_by_biosample),
        ("sra_run_accession", targets_by_run, external_by_run),
        ("assembly_accession", targets_by_assembly, external_by_assembly),
        ("isolate_identifier", targets_by_isolate_id, external_by_isolate_id),
    ):
        for identifier in sorted(set(development).intersection(external)):
            overlaps.append(
                {
                    "identifier_type": identifier_type,
                    "identifier": identifier,
                    "development_target_accessions_json": json.dumps(sorted(development[identifier]), separators=(",", ":")),
                    "external_isolate_ids_json": json.dumps(sorted(external[identifier]), separators=(",", ":")),
                }
            )
    counts = Counter(row["identifier_type"] for row in overlaps)
    summary = {
        "development_isolate_rows": len(isolates),
        "development_enrichment_rows": len(enrichment),
        "overlap_identifier_counts": dict(sorted(counts.items())),
        "overlap_rows": len(overlaps),
        "assembly_accession_test": "complete_for_registered_assemblies_not_supplement_wide",
        "registered_external_assembly_accessions_tested": len(assembly_rows),
        "read_hash_test": "not_testable_development_snapshot_has_no_read_hashes",
        "genomic_near_neighbor_test": "pending_sequence_based_comparison",
        "external_lock_status": (
            "blocked_exact_identifier_overlap"
            if overlaps
            else "blocked_genomic_disjointness_pending"
        ),
    }
    return overlaps, summary


def build_ecoli_ast_rows(
    supplement_rows: list[dict[str, str]],
    sequence_manifest: list[dict[str, str]],
) -> list[dict[str, str]]:
    sequence_by_id = {row["isolate_id"]: row for row in sequence_manifest}
    result = []
    for row in supplement_rows:
        if row["genome_species"] != "Escherichia coli":
            continue
        sequence = sequence_by_id[row["isolate ID"]]
        result.append(
            {
                "isolate_id": row["isolate ID"],
                "genome_species": row["genome_species"],
                "biosample_accession": sequence["biosample_accession"],
                "primary_short_read_run": sequence["primary_short_read_run"],
                "collection_date_raw": sequence["collection_date_raw"],
                "geo_loc_name_raw": sequence["geo_loc_name_raw"],
                "isolation_source_raw": sequence["isolation_source_raw"],
                "mean_total_coverage_raw": row["Mean_total_coverage"],
                "num_contigs_raw": row["Num_contigs"],
                "total_bases_raw": row["Total_bases"],
                "ceftriaxone_mic_raw": row["CTRX"],
                "ciprofloxacin_mic_raw": row["CPFX"],
                "gentamicin_mic_raw": row["GM"],
                "mlst_raw": row["MLST"],
                "sequence_link_status": sequence["sequence_link_status"],
            }
        )
    return result


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--supplement-xlsx", type=Path, required=True)
    parser.add_argument(
        "--inventory-only",
        action="store_true",
        help="Print workbook schema and the first three data rows as JSON.",
    )
    parser.add_argument("--sra-xml", type=Path)
    parser.add_argument("--sra-runinfo", type=Path)
    parser.add_argument("--sra-esearch", type=Path)
    parser.add_argument("--ncbi-datasets-assembly-report", type=Path)
    parser.add_argument("--development-isolates", type=Path)
    parser.add_argument("--development-enrichment", type=Path)
    parser.add_argument("--source-archive", type=Path)
    parser.add_argument("--pmc-oa-record", type=Path)
    parser.add_argument("--pmc-fulltext", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--snapshot-id")
    parser.add_argument("--retrieved-at")
    parser.add_argument(
        "--contact-email-supplied-at-runtime",
        action="store_true",
        help="Record that NCBI requests used a runtime contact without storing it.",
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    sheets = read_xlsx(args.supplement_xlsx)
    if args.inventory_only:
        print(json.dumps(workbook_inventory(sheets), indent=2, ensure_ascii=False))
        return 0
    required = {
        "sra_xml": args.sra_xml,
        "sra_runinfo": args.sra_runinfo,
        "sra_esearch": args.sra_esearch,
        "ncbi_datasets_assembly_report": args.ncbi_datasets_assembly_report,
        "development_isolates": args.development_isolates,
        "development_enrichment": args.development_enrichment,
        "source_archive": args.source_archive,
        "pmc_oa_record": args.pmc_oa_record,
        "pmc_fulltext": args.pmc_fulltext,
        "output_dir": args.output_dir,
        "snapshot_id": args.snapshot_id,
        "retrieved_at": args.retrieved_at,
    }
    missing = sorted(name for name, value in required.items() if not value)
    if missing:
        raise JarbsAuditError(f"Missing required audit arguments: {missing}")

    expected_supplement_bytes = 1493480
    expected_supplement_md5 = "00659b377a4830745048ba1643656276"
    if args.supplement_xlsx.stat().st_size != expected_supplement_bytes:
        raise JarbsAuditError("Supplement workbook size differs from the PMC record.")
    if md5_file(args.supplement_xlsx) != expected_supplement_md5:
        raise JarbsAuditError("Supplement workbook MD5 differs from the PMC record.")
    supplement_rows, supplement_summary = audit_supplement(sheets)
    run_rows, sra_summary = parse_sra_xml(args.sra_xml)
    runinfo_rows = _read_csv(args.sra_runinfo)
    runinfo_accessions = {row.get("Run", "") for row in runinfo_rows}
    xml_accessions = {row["run_accession"] for row in run_rows}
    if runinfo_accessions != xml_accessions:
        raise JarbsAuditError("SRA RunInfo and SRA XML accession sets differ.")
    search_payload = json.loads(args.sra_esearch.read_text(encoding="utf-8"))
    search_count = int(search_payload["esearchresult"]["count"])
    if search_count != len(run_rows):
        raise JarbsAuditError("SRA ESearch count does not match fetched run records.")

    sequence_manifest, linkage_summary = build_isolate_sequence_manifest(
        supplement_rows, run_rows
    )
    assembly_rows, assembly_summary = parse_assembly_report(
        args.ncbi_datasets_assembly_report
    )
    overlaps, overlap_summary = audit_development_overlap(
        sequence_manifest,
        assembly_rows,
        args.development_isolates,
        args.development_enrichment,
    )
    ecoli_rows = build_ecoli_ast_rows(supplement_rows, sequence_manifest)
    output_dir: Path = args.output_dir
    run_path = output_dir / "jarbs_sra_run_manifest.csv"
    isolate_path = output_dir / "jarbs_isolate_sequence_manifest.csv"
    ecoli_path = output_dir / "jarbs_ecoli_target_ast.csv"
    overlap_path = output_dir / "jarbs_development_exact_overlaps.csv"
    assembly_path = output_dir / "jarbs_registered_assembly_manifest.csv"
    audit_path = output_dir / "jarbs_external_audit.json"
    _write_csv(run_path, run_rows, list(run_rows[0]))
    _write_csv(isolate_path, sequence_manifest, list(sequence_manifest[0]))
    _write_csv(ecoli_path, ecoli_rows, list(ecoli_rows[0]))
    _write_csv(assembly_path, assembly_rows, list(assembly_rows[0]))
    _write_csv(
        overlap_path,
        overlaps,
        [
            "identifier_type",
            "identifier",
            "development_target_accessions_json",
            "external_isolate_ids_json",
        ],
    )
    audit = {
        "schema_version": "1.0.0",
        "snapshot_id": args.snapshot_id,
        "retrieved_at": args.retrieved_at,
        "source_identity": {
            "article_doi": "10.1038/s41467-023-43516-4",
            "pmcid": "PMC10698200",
            "bioproject": "PRJDB10842",
            "article_license": "CC BY 4.0",
            "supplement_publisher_md5": expected_supplement_md5,
            "supplement_publisher_size_bytes": expected_supplement_bytes,
        },
        "supplement": supplement_summary,
        "sra": sra_summary,
        "registered_assemblies": assembly_summary,
        "supplement_sra_linkage": linkage_summary,
        "development_overlap": overlap_summary,
        "scientific_boundary": {
            "role": "conditional_external_challenge_candidate",
            "phenotype_population": "resistance_enriched_not_prevalence_representative",
            "ast_testing_date": "absent",
            "hospital_site": "absent_from_public_isolate_level_files",
            "clinical_indication": "absent",
            "method_provenance": "paper_level_central_broth_microdilution_not_row_level",
            "breakpoint_version": "paper_reports_CLSI_2021_with_PIP_TAZ_2022_exception",
            "external_lock_status": overlap_summary["external_lock_status"],
        },
    }
    _write_json(audit_path, audit)

    artifact_paths = [
        args.supplement_xlsx,
        args.source_archive,
        args.pmc_oa_record,
        args.pmc_fulltext,
        args.sra_esearch,
        args.sra_runinfo,
        args.sra_xml,
        args.ncbi_datasets_assembly_report,
        args.development_isolates,
        args.development_enrichment,
        run_path,
        isolate_path,
        ecoli_path,
        assembly_path,
        overlap_path,
        audit_path,
    ]
    manifest = {
        "schema_version": "1.0.0",
        "snapshot_id": args.snapshot_id,
        "retrieved_at": args.retrieved_at,
        "contact_email_supplied_at_runtime": args.contact_email_supplied_at_runtime,
        "contact_email_persisted": False,
        "artifacts": [
            {
                "path": path.resolve().relative_to(Path.cwd().resolve()).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in artifact_paths
        ],
    }
    manifest_path = output_dir / "jarbs_source_manifest.json"
    _write_json(manifest_path, manifest)
    print(json.dumps(audit, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
