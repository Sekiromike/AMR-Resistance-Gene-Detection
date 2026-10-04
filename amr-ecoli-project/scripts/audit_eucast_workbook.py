"""Audit selected EUCAST v16.1 OOXML cells without modifying the workbook.

The output preserves visible values, rich-text run properties, cell styles,
merged ranges, comments, ATU columns, and section-note anchors.  It is review
evidence only and never emits an approved breakpoint-rule table.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from typing import Any, Iterable


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PACKAGE_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
CELL_RE = re.compile(r"^([A-Z]+)([1-9][0-9]*)$")
TARGET_SHEET = "Enterobacterales"
TARGET_ROWS = (35, 36, 37, 55, 56, 79, 80, 81, 82, 83, 93, 94, 95, 97, 98)


class EucastAuditError(RuntimeError):
    """Raised when the frozen workbook cannot be authenticated or read exactly."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _column_number(name: str) -> int:
    value = 0
    for character in name:
        value = value * 26 + ord(character) - ord("A") + 1
    return value


def _cell_coordinates(reference: str) -> tuple[int, int]:
    match = CELL_RE.fullmatch(reference)
    if match is None:
        raise EucastAuditError(f"Invalid OOXML cell reference: {reference!r}")
    return _column_number(match.group(1)), int(match.group(2))


def _target_cells() -> set[str]:
    return {f"{column}{row}" for row in TARGET_ROWS for column in "ABCDEFGHI"}


def _rich_text(element: ET.Element) -> dict[str, Any]:
    runs = []
    run_elements = element.findall(f"{{{MAIN_NS}}}r")
    for run in run_elements:
        text = "".join(node.text or "" for node in run.iter(f"{{{MAIN_NS}}}t"))
        properties = run.find(f"{{{MAIN_NS}}}rPr")
        details: dict[str, Any] = {"text": text}
        if properties is not None:
            for tag, key in (("b", "bold"), ("i", "italic"), ("strike", "strike")):
                if properties.find(f"{{{MAIN_NS}}}{tag}") is not None:
                    details[key] = True
            for tag, key in (
                ("vertAlign", "vertical_alignment"),
                ("sz", "size"),
                ("rFont", "font"),
                ("family", "family"),
                ("charset", "charset"),
            ):
                node = properties.find(f"{{{MAIN_NS}}}{tag}")
                if node is not None and "val" in node.attrib:
                    details[key] = node.attrib["val"]
            color = properties.find(f"{{{MAIN_NS}}}color")
            if color is not None:
                details["color"] = dict(sorted(color.attrib.items()))
        runs.append(details)
    text = "".join(node.text or "" for node in element.iter(f"{{{MAIN_NS}}}t"))
    return {"text": text, "runs": runs}


def _shared_strings(archive: zipfile.ZipFile) -> list[dict[str, Any]]:
    try:
        root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    except KeyError:
        return []
    return [_rich_text(item) for item in root.findall(f"{{{MAIN_NS}}}si")]


def _sheet_path(archive: zipfile.ZipFile, sheet_name: str) -> str:
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {
        relation.attrib["Id"]: relation.attrib["Target"]
        for relation in relationships.findall(f"{{{PACKAGE_REL_NS}}}Relationship")
    }
    for sheet in workbook.findall(f".//{{{MAIN_NS}}}sheet"):
        if sheet.attrib.get("name") != sheet_name:
            continue
        relationship_id = sheet.attrib[f"{{{REL_NS}}}id"]
        target = targets[relationship_id].replace("\\", "/").lstrip("/")
        return target if target.startswith("xl/") else f"xl/{target}"
    raise EucastAuditError(f"Worksheet not found: {sheet_name}")


def _styles(archive: zipfile.ZipFile) -> list[dict[str, Any]]:
    root = ET.fromstring(archive.read("xl/styles.xml"))
    custom_formats = {
        node.attrib["numFmtId"]: node.attrib["formatCode"]
        for node in root.findall(f".//{{{MAIN_NS}}}numFmts/{{{MAIN_NS}}}numFmt")
    }
    cell_xfs = root.find(f"{{{MAIN_NS}}}cellXfs")
    if cell_xfs is None:
        raise EucastAuditError("Workbook has no cellXfs style table")
    result = []
    for index, xf in enumerate(cell_xfs.findall(f"{{{MAIN_NS}}}xf")):
        style: dict[str, Any] = {"style_id": index, **dict(sorted(xf.attrib.items()))}
        number_format_id = xf.attrib.get("numFmtId")
        if number_format_id in custom_formats:
            style["number_format_code"] = custom_formats[number_format_id]
        alignment = xf.find(f"{{{MAIN_NS}}}alignment")
        if alignment is not None:
            style["alignment"] = dict(sorted(alignment.attrib.items()))
        protection = xf.find(f"{{{MAIN_NS}}}protection")
        if protection is not None:
            style["protection"] = dict(sorted(protection.attrib.items()))
        result.append(style)
    return result


def _cell_record(
    cell: ET.Element,
    shared_strings: list[dict[str, Any]],
    styles: list[dict[str, Any]],
) -> dict[str, Any]:
    reference = cell.attrib["r"]
    cell_type = cell.attrib.get("t", "")
    value_node = cell.find(f"{{{MAIN_NS}}}v")
    raw = "" if value_node is None or value_node.text is None else value_node.text
    rich_text: dict[str, Any] | None = None
    if cell_type == "s" and raw:
        try:
            rich_text = shared_strings[int(raw)]
        except (IndexError, ValueError) as exc:
            raise EucastAuditError(f"Invalid shared-string index in {reference}: {raw!r}") from exc
        value = rich_text["text"]
    elif cell_type == "inlineStr":
        inline = cell.find(f"{{{MAIN_NS}}}is")
        rich_text = _rich_text(inline) if inline is not None else {"text": "", "runs": []}
        value = rich_text["text"]
    elif cell_type == "b":
        value = "TRUE" if raw == "1" else "FALSE"
    else:
        value = raw
    style_id = int(cell.attrib.get("s", "0"))
    if style_id >= len(styles):
        raise EucastAuditError(f"Invalid style ID in {reference}: {style_id}")
    formula = cell.find(f"{{{MAIN_NS}}}f")
    record: dict[str, Any] = {
        "reference": reference,
        "value": value,
        "cell_type": cell_type or "numeric_or_formula",
        "style": styles[style_id],
    }
    if formula is not None:
        record["formula"] = formula.text or ""
    if rich_text is not None and rich_text["runs"]:
        record["rich_text_runs"] = rich_text["runs"]
    return record


def _range_contains(cell_reference: str, range_reference: str) -> bool:
    start, _, end = range_reference.partition(":")
    end = end or start
    cell_column, cell_row = _cell_coordinates(cell_reference)
    start_column, start_row = _cell_coordinates(start)
    end_column, end_row = _cell_coordinates(end)
    return start_column <= cell_column <= end_column and start_row <= cell_row <= end_row


def _comments(archive: zipfile.ZipFile, sheet_path: str) -> list[dict[str, str]]:
    sheet_file = Path(sheet_path)
    relationships_path = (
        sheet_file.parent / "_rels" / f"{sheet_file.name}.rels"
    ).as_posix()
    try:
        relationships = ET.fromstring(archive.read(relationships_path))
    except KeyError:
        return []
    comments_target = None
    for relation in relationships.findall(f"{{{PACKAGE_REL_NS}}}Relationship"):
        if relation.attrib.get("Type", "").endswith("/comments"):
            comments_target = relation.attrib["Target"].replace("\\", "/")
            break
    if comments_target is None:
        return []
    comments_path = (sheet_file.parent / comments_target).as_posix()
    comments_root = ET.fromstring(archive.read(comments_path))
    authors = [
        node.text or "" for node in comments_root.findall(f".//{{{MAIN_NS}}}authors/{{{MAIN_NS}}}author")
    ]
    result = []
    for comment in comments_root.findall(f".//{{{MAIN_NS}}}commentList/{{{MAIN_NS}}}comment"):
        author_id = int(comment.attrib.get("authorId", "0"))
        result.append(
            {
                "reference": comment.attrib["ref"],
                "author": authors[author_id] if author_id < len(authors) else "",
                "text": "".join(
                    node.text or "" for node in comment.iter(f"{{{MAIN_NS}}}t")
                ),
            }
        )
    return result


def audit_workbook(workbook_path: Path, manifest_path: Path) -> dict[str, Any]:
    if not workbook_path.is_file() or not zipfile.is_zipfile(workbook_path):
        raise EucastAuditError(f"Not a readable OOXML workbook: {workbook_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    observed_hash = sha256_file(workbook_path)
    observed_bytes = workbook_path.stat().st_size
    if manifest.get("artifact") != workbook_path.name:
        raise EucastAuditError("Acquisition manifest names a different workbook")
    if manifest.get("sha256") != observed_hash or manifest.get("bytes") != observed_bytes:
        raise EucastAuditError("Workbook hash or byte count does not match acquisition manifest")

    targets = _target_cells()
    with zipfile.ZipFile(workbook_path) as archive:
        bad_member = archive.testzip()
        if bad_member is not None:
            raise EucastAuditError(f"Workbook ZIP member failed integrity testing: {bad_member}")
        expected_members = manifest.get("integrity", {}).get("xlsx_member_count")
        if expected_members != len(archive.namelist()):
            raise EucastAuditError("Workbook member count does not match acquisition manifest")
        sheet_path = _sheet_path(archive, TARGET_SHEET)
        sheet = ET.fromstring(archive.read(sheet_path))
        shared_strings = _shared_strings(archive)
        styles = _styles(archive)
        cell_records = {
            cell.attrib["r"]: _cell_record(cell, shared_strings, styles)
            for cell in sheet.findall(f".//{{{MAIN_NS}}}sheetData/{{{MAIN_NS}}}row/{{{MAIN_NS}}}c")
            if cell.attrib.get("r") in targets
        }
        for reference in sorted(targets, key=lambda value: _cell_coordinates(value)[::-1]):
            cell_records.setdefault(
                reference,
                {"reference": reference, "value": "", "cell_type": "absent", "style": None},
            )
        merged_ranges = [
            node.attrib["ref"]
            for node in sheet.findall(f".//{{{MAIN_NS}}}mergeCells/{{{MAIN_NS}}}mergeCell")
            if any(_range_contains(reference, node.attrib["ref"]) for reference in targets)
        ]
        comments = [
            comment for comment in _comments(archive, sheet_path) if comment["reference"] in targets
        ]

    values = {reference: record["value"] for reference, record in cell_records.items()}
    required_values = {
        "B55": "1",
        "C55": "2",
        "F55": "27",
        "G55": "24",
        "B56": "1",
        "C56": "1",
        "D82": "0.5",
        "H82": "22-24",
        "F83": "NoteB",
        "G83": "NoteB",
        "B97": "(2)1",
        "F97": "(17)A",
        "B98": "2",
        "F98": "17",
    }
    mismatches = {
        reference: {"expected": expected, "observed": values.get(reference, "")}
        for reference, expected in required_values.items()
        if values.get(reference, "") != expected
    }
    if mismatches:
        raise EucastAuditError(f"Targeted EUCAST cell contract changed: {mismatches}")

    ordered_cells = [
        cell_records[reference]
        for reference in sorted(cell_records, key=lambda value: _cell_coordinates(value)[::-1])
    ]
    return {
        "schema_version": "1.0.0",
        "scientific_status": "CANDIDATE_TRANSCRIPTION_AUDIT_ONLY",
        "approval_status": "INDEPENDENT_REVIEW_REQUIRED",
        "source": {
            "workbook": str(workbook_path),
            "sha256": observed_hash,
            "bytes": observed_bytes,
            "acquisition_manifest": str(manifest_path),
            "acquisition_manifest_sha256": sha256_file(manifest_path),
            "xlsx_member_count": expected_members,
            "zip_integrity": "PASS",
            "worksheet": TARGET_SHEET,
            "worksheet_path": sheet_path,
        },
        "target_rows": list(TARGET_ROWS),
        "cells": ordered_cells,
        "merged_ranges_intersecting_targets": sorted(merged_ranges),
        "comments_on_target_cells": comments,
        "findings": {
            "ceftriaxone_indication_rows_are_distinct": values["A55"] != values["A56"],
            "ciprofloxacin_non_meningitis_has_mic_atu": values["D82"],
            "ciprofloxacin_non_meningitis_has_zone_atu": values["H82"],
            "ciprofloxacin_meningitis_requires_note_b_disk_method": {
                "disk_content": values["E83"],
                "zone_s": values["F83"],
                "zone_r": values["G83"],
                "section_notes": values["I81"],
            },
            "gentamicin_systemic_values_are_bracketed": {
                "mic_s": values["B97"],
                "mic_r": values["C97"],
                "zone_s": values["F97"],
                "zone_r": values["G97"],
                "section_notes": values["I95"],
            },
            "gentamicin_urinary_values_are_unbracketed": {
                "mic_s": values["B98"],
                "mic_r": values["C98"],
                "zone_s": values["F98"],
                "zone_r": values["G98"],
            },
        },
        "limitations": [
            "This is a deterministic OOXML audit, not an independent human transcription approval.",
            "No breakpoint rule CSV is created or approved.",
            "Clinical indication, AST method, test date, and disk potency remain source-data requirements and are not inferred.",
            "ATU and bracketed values are retained as review constraints and are not converted to ordinary S/I/R thresholds.",
        ],
    }


def write_audit(payload: dict[str, Any], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(output_path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--acquisition-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        payload = audit_workbook(args.workbook, args.acquisition_manifest)
        write_audit(payload, args.output)
    except (EucastAuditError, OSError, json.JSONDecodeError, zipfile.BadZipFile) as exc:
        raise SystemExit(f"EUCAST workbook audit failed: {exc}") from exc
    print(
        f"EUCAST workbook audit complete: {len(payload['cells'])} target cells; "
        f"status={payload['approval_status']}; output={args.output}"
    )


if __name__ == "__main__":
    main()
