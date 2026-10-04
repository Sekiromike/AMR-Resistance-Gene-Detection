"""Independently interpret quantitative AST using frozen breakpoint rules.

The normalized AST input is never modified in place. Each emitted category is
linked to a unique machine-readable rule and to the SHA-256 of the original
breakpoint artifact from which that rule was transcribed. Censored measurements
are represented as intervals; a row is interpreted only when the entire interval
belongs to one S, I, or R region.

The accompanying normalization manifest is mandatory and authenticates the AST
file, acquisition provenance, mapping version, transform script, and declared
source-stage eligibility count. Rows marked ineligible cannot reach rule matching.

Required AST columns:
    source_ast_record_id, antibiotic, submitted_ast_category,
    ast_measurement_type, measurement_sign, ast_value, ast_unit, ast_method,
    disk_content, ast_testing_date, breakpoint_eligible,
    breakpoint_ineligibility_reasons

Required rule columns:
    breakpoint_rule_id, breakpoint_artifact_sha256, breakpoint_standard,
    breakpoint_version, antibiotic, ast_measurement_type, ast_unit,
    susceptible_breakpoint, resistant_breakpoint

Optional rule constraints are ``ast_method``, ``disk_content``,
``clinical_indication``, and ``scientific_name``.
When a constraint is populated in a rule, the AST row must contain an exact
case-insensitive, whitespace-normalized match. No units or metadata are inferred.

Breakpoint semantics:
    MIC:  S <= susceptible_breakpoint; R > resistant_breakpoint
    ZONE: S >= susceptible_breakpoint; R < resistant_breakpoint

Values strictly between those regions are I. The two breakpoint thresholds may
be equal, in which case the intermediate region is empty.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import NamedTuple


AST_REQUIRED_COLUMNS = (
    "source_ast_record_id",
    "antibiotic",
    "submitted_ast_category",
    "ast_measurement_type",
    "measurement_sign",
    "ast_value",
    "ast_unit",
    "ast_method",
    "disk_content",
    "ast_testing_date",
    "breakpoint_eligible",
    "breakpoint_ineligibility_reasons",
)

RULE_REQUIRED_COLUMNS = (
    "breakpoint_rule_id",
    "breakpoint_artifact_sha256",
    "breakpoint_standard",
    "breakpoint_version",
    "antibiotic",
    "ast_measurement_type",
    "ast_unit",
    "susceptible_breakpoint",
    "resistant_breakpoint",
)

OPTIONAL_CONSTRAINTS = (
    "ast_method",
    "disk_content",
    "clinical_indication",
    "scientific_name",
)
CONTEXT_REQUIRED_COLUMNS = (
    "isolate_id",
    "biosample_accession",
    "assembly_accession",
    "clinical_indication",
)
CONTEXT_ADDED_COLUMNS = ("clinical_indication",)
VALID_COMPARATORS = {"=", "<", "<=", ">", ">="}
VALID_MEASUREMENT_TYPES = {"MIC", "ZONE"}
MIC_AST_METHODS = {
    "agar dilution",
    "broth macrodilution",
    "broth microdilution",
    "gradient diffusion",
}
ZONE_AST_METHODS = {"disk diffusion"}
DISK_CONTENT_PATTERN = re.compile(
    r"^\d+(?:\.\d+)?(?:\s*/\s*\d+(?:\.\d+)?)?\s*(?:ug|mcg|µg|μg)$",
    flags=re.IGNORECASE,
)
ADDED_OUTPUT_COLUMNS = (
    "ast_category",
    "submitted_category_concordance",
    "breakpoint_rule_id",
    "breakpoint_rule_row_sha256",
    "breakpoint_artifact_sha256",
    "breakpoint_standard",
    "breakpoint_version",
)
EXCLUSION_COLUMNS = ("exclusion_reason", "exclusion_detail")


class BreakpointError(ValueError):
    """Raised when an input contract or frozen-artifact invariant is violated."""


class Interval(NamedTuple):
    lower: Decimal | None
    lower_inclusive: bool
    upper: Decimal | None
    upper_inclusive: bool


class BreakpointRule(NamedTuple):
    rule_id: str
    artifact_sha256: str
    standard: str
    version: str
    antibiotic: str
    measurement_type: str
    unit: str
    susceptible_breakpoint: Decimal
    resistant_breakpoint: Decimal
    constraints: dict[str, str]
    row_sha256: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _blank(value: object) -> bool:
    return value is None or str(value).strip().lower() in {"", "nan", "none", "null"}


def _match_text(value: object) -> str:
    """Normalize only case and whitespace; do not map biological concepts."""
    return " ".join(str(value).strip().casefold().split())


def _valid_ast_method(value: object, measurement_type: object) -> bool:
    method = _match_text(value)
    kind = str(measurement_type or "").strip().upper()
    if kind == "MIC":
        return method in MIC_AST_METHODS
    if kind == "ZONE":
        return method in ZONE_AST_METHODS
    return False


def _valid_disk_content(value: object) -> bool:
    text = str(value or "").strip()
    if not DISK_CONTENT_PATTERN.fullmatch(text):
        return False
    return all(float(part) > 0 for part in re.findall(r"\d+(?:\.\d+)?", text))


def _valid_iso_date_or_datetime(value: object) -> bool:
    text = str(value or "").strip()
    try:
        datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return False
    return bool(re.match(r"^\d{4}-\d{2}-\d{2}(?:T|$)", text))


def _decimal(value: object, field: str) -> Decimal:
    if _blank(value):
        raise BreakpointError(f"Missing {field}.")
    try:
        parsed = Decimal(str(value).strip())
    except InvalidOperation as exc:
        raise BreakpointError(f"Invalid numeric {field}: {value!r}.") from exc
    if not parsed.is_finite():
        raise BreakpointError(f"Non-finite {field}: {value!r}.")
    return parsed


def _canonical_rule_sha256(row: dict[str, str]) -> str:
    canonical = {
        key: "" if value is None else str(value).strip()
        for key, value in sorted(row.items())
        if key != "breakpoint_rule_row_sha256"
    }
    payload = json.dumps(canonical, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise BreakpointError(f"CSV input does not exist: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise BreakpointError(f"CSV has no header: {path}")
        fields = list(reader.fieldnames)
        if len(fields) != len(set(fields)):
            raise BreakpointError(f"CSV has duplicate column names: {path}")
        rows = [dict(row) for row in reader]
    return fields, rows


def verify_normalization_manifest(
    manifest_path: Path,
    ast_path: Path,
    ast_rows: list[dict[str, str]],
    *,
    allow_unverified: bool = False,
) -> dict[str, object]:
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BreakpointError(f"Cannot read normalization manifest {manifest_path}: {exc}") from exc
    scientific_status = str(manifest.get("scientific_status") or "")
    allowed_statuses = {"SOURCE_NORMALIZATION_ONLY"}
    if allow_unverified:
        allowed_statuses.add("SMOKE_TEST_ONLY_UNVERIFIED_INPUT")
    if scientific_status not in allowed_statuses:
        raise BreakpointError(
            "Normalization manifest is not authenticated publication-source normalization."
        )
    mapping_version = str(manifest.get("mapping_version") or "")
    if not re.fullmatch(r"\d+\.\d+\.\d+", mapping_version):
        raise BreakpointError("Normalization manifest has no semantic mapping version.")
    output = ((manifest.get("outputs") or {}).get("normalized") or {})
    expected_hash = str(output.get("sha256") or "")
    expected_rows = output.get("rows")
    expected_bytes = output.get("bytes")
    if not expected_hash or expected_rows is None or expected_bytes is None:
        raise BreakpointError("Normalization manifest lacks normalized output hash/rows/bytes.")
    observed_hash = sha256_file(ast_path)
    if expected_hash != observed_hash:
        raise BreakpointError(
            f"Normalized AST hash mismatch: manifest={expected_hash}, observed={observed_hash}"
        )
    if int(expected_rows) != len(ast_rows) or int(expected_bytes) != ast_path.stat().st_size:
        raise BreakpointError("Normalization manifest row or byte count does not match AST input.")
    observed_eligible = sum(
        str(row.get("breakpoint_eligible", "") or "").strip().casefold() == "true"
        for row in ast_rows
    )
    if int(manifest.get("n_breakpoint_eligible_rows", -1)) != observed_eligible:
        raise BreakpointError("Normalization manifest breakpoint-eligibility count mismatch.")
    normalizer = manifest.get("normalizer") or {}
    script_hash = str(normalizer.get("script_sha256") or "")
    if not re.fullmatch(r"[0-9a-f]{64}", script_hash):
        raise BreakpointError("Normalization manifest lacks a valid normalizer script hash.")
    acquisition = manifest.get("acquisition_provenance")
    acquisition_hash: str | None = None
    if scientific_status == "SOURCE_NORMALIZATION_ONLY":
        if not isinstance(acquisition, dict) or not re.fullmatch(
            r"[0-9a-f]{64}", str(acquisition.get("sha256") or "")
        ):
            raise BreakpointError("Normalization manifest lacks authenticated acquisition provenance.")
        acquisition_hash = str(acquisition["sha256"])
    return {
        "path": str(manifest_path),
        "sha256": sha256_file(manifest_path),
        "mapping_version": mapping_version,
        "normalizer_script_sha256": script_hash,
        "normalized_ast_sha256": observed_hash,
        "n_breakpoint_eligible_rows": observed_eligible,
        "scientific_status": scientific_status,
        "acquisition_manifest_sha256": acquisition_hash,
    }


def enrich_with_context(
    ast_rows: list[dict[str, str]], context_rows: list[dict[str, str]]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Join indication context without repairing or inferring isolate identity."""
    by_isolate: dict[str, list[dict[str, str]]] = {}
    for context in context_rows:
        isolate_id = str(context.get("isolate_id", "")).strip()
        if not isolate_id:
            raise BreakpointError("Context metadata contains a blank isolate_id.")
        by_isolate.setdefault(isolate_id, []).append(context)
    duplicate_ids = sorted(key for key, rows in by_isolate.items() if len(rows) != 1)
    if duplicate_ids:
        raise BreakpointError(f"Context metadata isolate_id is not unique: {duplicate_ids[:20]}")

    enriched: list[dict[str, str]] = []
    excluded: list[dict[str, str]] = []
    for source in ast_rows:
        row = dict(source)
        isolate_id = str(row.get("isolate_id", "")).strip()
        matches = by_isolate.get(isolate_id, [])
        if not matches:
            row["exclusion_reason"] = "breakpoint_context_not_found"
            row["exclusion_detail"] = "No exact isolate_id exists in context metadata."
            excluded.append(row)
            continue
        context = matches[0]
        mismatches = [
            field
            for field in ("biosample_accession", "assembly_accession")
            if str(row.get(field, "")).strip() != str(context.get(field, "")).strip()
        ]
        if mismatches:
            row["exclusion_reason"] = "breakpoint_context_identity_mismatch"
            row["exclusion_detail"] = "AST and context disagree for: " + ",".join(mismatches)
            excluded.append(row)
            continue
        for field in CONTEXT_ADDED_COLUMNS:
            observed = row.get(field, "")
            supplied = context.get(field, "")
            if not _blank(observed) and _match_text(observed) != _match_text(supplied):
                row["exclusion_reason"] = "breakpoint_context_value_mismatch"
                row["exclusion_detail"] = f"AST and context disagree for: {field}"
                excluded.append(row)
                break
            row[field] = supplied
        else:
            enriched.append(row)
    return enriched, excluded


def _require_columns(fields: list[str], required: tuple[str, ...], label: str) -> None:
    missing = sorted(set(required).difference(fields))
    if missing:
        raise BreakpointError(f"{label} is missing required columns: {missing}")


def load_breakpoint_rules(rules_path: Path, artifact_path: Path) -> list[BreakpointRule]:
    """Load rules only after verifying every row against the frozen artifact."""
    if not artifact_path.is_file():
        raise BreakpointError(f"Frozen breakpoint artifact does not exist: {artifact_path}")
    observed_artifact_sha = sha256_file(artifact_path)
    fields, rows = _read_csv(rules_path)
    _require_columns(fields, RULE_REQUIRED_COLUMNS, "Breakpoint rules CSV")
    if not rows:
        raise BreakpointError("Breakpoint rules CSV has no rows.")

    loaded: list[BreakpointRule] = []
    seen_rule_ids: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        rule_id = str(row["breakpoint_rule_id"] or "").strip()
        if not rule_id:
            raise BreakpointError(f"Rule line {line_number} has no breakpoint_rule_id.")
        if rule_id in seen_rule_ids:
            raise BreakpointError(f"Duplicate breakpoint_rule_id: {rule_id}")
        seen_rule_ids.add(rule_id)

        declared_sha = str(row["breakpoint_artifact_sha256"] or "").strip().lower()
        if declared_sha != observed_artifact_sha:
            raise BreakpointError(
                f"Rule {rule_id} artifact SHA-256 mismatch: "
                f"declared={declared_sha or '<missing>'}, observed={observed_artifact_sha}"
            )

        standard = str(row["breakpoint_standard"] or "").strip()
        version = str(row["breakpoint_version"] or "").strip()
        antibiotic = str(row["antibiotic"] or "").strip()
        measurement_type = str(row["ast_measurement_type"] or "").strip().upper()
        unit = str(row["ast_unit"] or "").strip()
        for field, value in (
            ("breakpoint_standard", standard),
            ("breakpoint_version", version),
            ("antibiotic", antibiotic),
            ("ast_unit", unit),
        ):
            if not value:
                raise BreakpointError(f"Rule {rule_id} has no {field}.")
        if measurement_type not in VALID_MEASUREMENT_TYPES:
            raise BreakpointError(
                f"Rule {rule_id} has unsupported ast_measurement_type {measurement_type!r}."
            )

        susceptible = _decimal(row["susceptible_breakpoint"], "susceptible_breakpoint")
        resistant = _decimal(row["resistant_breakpoint"], "resistant_breakpoint")
        if measurement_type == "MIC" and resistant < susceptible:
            raise BreakpointError(
                f"Rule {rule_id} has MIC resistant breakpoint below susceptible breakpoint."
            )
        if measurement_type == "ZONE" and susceptible < resistant:
            raise BreakpointError(
                f"Rule {rule_id} has ZONE susceptible breakpoint below resistant breakpoint."
            )

        constraints = {
            field: str(row.get(field, "") or "").strip()
            for field in OPTIONAL_CONSTRAINTS
            if not _blank(row.get(field, ""))
        }
        loaded.append(
            BreakpointRule(
                rule_id=rule_id,
                artifact_sha256=declared_sha,
                standard=standard,
                version=version,
                antibiotic=antibiotic,
                measurement_type=measurement_type,
                unit=unit,
                susceptible_breakpoint=susceptible,
                resistant_breakpoint=resistant,
                constraints=constraints,
                row_sha256=_canonical_rule_sha256(row),
            )
        )
    return loaded


def measurement_interval(comparator: str, value: Decimal) -> Interval:
    if comparator == "=":
        return Interval(value, True, value, True)
    if comparator == "<":
        return Interval(None, False, value, False)
    if comparator == "<=":
        return Interval(None, False, value, True)
    if comparator == ">":
        return Interval(value, False, None, False)
    if comparator == ">=":
        return Interval(value, True, None, False)
    raise BreakpointError(f"Unsupported comparator: {comparator!r}")


def _is_subset(inner: Interval, outer: Interval) -> bool:
    if outer.lower is not None:
        if inner.lower is None or inner.lower < outer.lower:
            return False
        if (
            inner.lower == outer.lower
            and inner.lower_inclusive
            and not outer.lower_inclusive
        ):
            return False
    if outer.upper is not None:
        if inner.upper is None or inner.upper > outer.upper:
            return False
        if (
            inner.upper == outer.upper
            and inner.upper_inclusive
            and not outer.upper_inclusive
        ):
            return False
    return True


def interpret_interval(interval: Interval, rule: BreakpointRule) -> str | None:
    """Return S/I/R only if the complete measurement interval implies it."""
    susceptible = rule.susceptible_breakpoint
    resistant = rule.resistant_breakpoint
    if rule.measurement_type == "MIC":
        regions = {
            "S": Interval(None, False, susceptible, True),
            "I": Interval(susceptible, False, resistant, True),
            "R": Interval(resistant, False, None, False),
        }
    else:
        regions = {
            "R": Interval(None, False, resistant, False),
            "I": Interval(resistant, True, susceptible, False),
            "S": Interval(susceptible, True, None, False),
        }
    matches = [category for category, region in regions.items() if _is_subset(interval, region)]
    return matches[0] if len(matches) == 1 else None


def _base_rule_matches(row: dict[str, str], rule: BreakpointRule) -> bool:
    return (
        _match_text(row.get("antibiotic", "")) == _match_text(rule.antibiotic)
        and str(row.get("ast_measurement_type", "")).strip().upper() == rule.measurement_type
        and _match_text(row.get("ast_unit", "")) == _match_text(rule.unit)
    )


def _select_rule(
    row: dict[str, str], rules: list[BreakpointRule]
) -> tuple[BreakpointRule | None, str | None, str]:
    base = [rule for rule in rules if _base_rule_matches(row, rule)]
    if not base:
        return None, "no_breakpoint_rule", "No rule matches antibiotic, measurement type, and unit."

    matches: list[BreakpointRule] = []
    missing_constraints: set[str] = set()
    for rule in base:
        rule_matches = True
        for field, expected in rule.constraints.items():
            observed = row.get(field, "")
            if _blank(observed):
                missing_constraints.add(field)
                rule_matches = False
            elif _match_text(observed) != _match_text(expected):
                rule_matches = False
        if rule_matches:
            matches.append(rule)

    if len(matches) == 1:
        return matches[0], None, ""
    if len(matches) > 1:
        ids = ",".join(sorted(rule.rule_id for rule in matches))
        return None, "ambiguous_breakpoint_rules", f"Multiple rules match: {ids}"
    if missing_constraints:
        fields = ",".join(sorted(missing_constraints))
        return (
            None,
            "missing_breakpoint_constraint_metadata",
            f"Rule-constrained metadata are absent: {fields}",
        )
    return None, "no_matching_breakpoint_rule", "Assay constraints do not match any candidate rule."


def _concordance(submitted: object, interpreted: str) -> str:
    if _blank(submitted):
        return "NOT_REPORTED"
    return "CONCORDANT" if str(submitted).strip().upper() == interpreted else "DISCORDANT"


def interpret_rows(
    rows: list[dict[str, str]], rules: list[BreakpointRule]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    interpreted: list[dict[str, str]] = []
    excluded: list[dict[str, str]] = []

    for source_row in rows:
        row = dict(source_row)
        error: tuple[str, str] | None = None
        eligibility = str(row.get("breakpoint_eligible", "") or "").strip().casefold()
        eligibility_reasons = str(
            row.get("breakpoint_ineligibility_reasons", "") or ""
        ).strip()
        if eligibility not in {"true", "false"}:
            error = (
                "invalid_breakpoint_eligibility_flag",
                "breakpoint_eligible must be exactly true or false.",
            )
        elif eligibility == "false":
            error = (
                "source_stage_breakpoint_ineligible",
                eligibility_reasons or "Source-stage eligibility failed without a reason code.",
            )
        elif eligibility_reasons:
            error = (
                "inconsistent_breakpoint_eligibility_state",
                "An eligible row must not contain breakpoint ineligibility reasons.",
            )
        elif not _valid_ast_method(
            row.get("ast_method", ""), row.get("ast_measurement_type", "")
        ):
            error = (
                "invalid_breakpoint_eligibility_state",
                "AST method is absent, unsupported, or inconsistent with measurement type.",
            )
        elif not _valid_iso_date_or_datetime(row.get("ast_testing_date", "")):
            error = (
                "invalid_breakpoint_eligibility_state",
                "AST testing date is absent or is not an ISO date/datetime.",
            )
        elif (
            str(row.get("ast_measurement_type", "") or "").strip().upper() == "ZONE"
            and not _valid_disk_content(row.get("disk_content", ""))
        ):
            error = (
                "invalid_breakpoint_eligibility_state",
                "A zone measurement requires a positive disk potency with an explicit mass unit.",
            )
        for field in ("source_ast_record_id", "antibiotic", "ast_measurement_type", "ast_unit"):
            if error is None and _blank(row.get(field, "")):
                error = (f"missing_{field}", f"Required AST field {field} is absent.")
                break

        comparator = str(row.get("measurement_sign", "") or "").strip()
        if error is None and not comparator:
            error = ("missing_measurement_sign", "A comparator was not reported.")
        elif error is None and comparator not in VALID_COMPARATORS:
            error = ("invalid_measurement_sign", f"Unsupported comparator: {comparator!r}")

        value: Decimal | None = None
        if error is None:
            try:
                value = _decimal(row.get("ast_value", ""), "ast_value")
            except BreakpointError as exc:
                error = ("invalid_ast_value", str(exc))

        rule: BreakpointRule | None = None
        if error is None:
            rule, reason, detail = _select_rule(row, rules)
            if reason is not None:
                error = (reason, detail)

        category: str | None = None
        if error is None:
            assert value is not None and rule is not None
            category = interpret_interval(measurement_interval(comparator, value), rule)
            if category is None:
                error = (
                    "censored_measurement_crosses_breakpoints",
                    "The complete censored interval does not imply one S/I/R category.",
                )

        if error is not None:
            row["exclusion_reason"], row["exclusion_detail"] = error
            excluded.append(row)
            continue

        assert rule is not None and category is not None
        row.update(
            {
                "ast_category": category,
                "submitted_category_concordance": _concordance(
                    row.get("submitted_ast_category", ""), category
                ),
                "breakpoint_rule_id": rule.rule_id,
                "breakpoint_rule_row_sha256": rule.row_sha256,
                "breakpoint_artifact_sha256": rule.artifact_sha256,
                "breakpoint_standard": rule.standard,
                "breakpoint_version": rule.version,
            }
        )
        interpreted.append(row)
    return interpreted, excluded


def _write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def run_interpretation(
    ast_path: Path,
    normalization_manifest_path: Path,
    rules_path: Path,
    artifact_path: Path,
    output_path: Path,
    exclusions_path: Path,
    manifest_path: Path,
    context_path: Path | None = None,
    allow_unverified_normalization: bool = False,
) -> dict[str, object]:
    """Validate inputs, interpret rows, and atomically write audit artifacts."""
    resolved_paths = [
        path.resolve()
        for path in (ast_path, normalization_manifest_path, rules_path, artifact_path)
    ]
    if context_path is not None:
        resolved_paths.append(context_path.resolve())
    output_paths = [path.resolve() for path in (output_path, exclusions_path, manifest_path)]
    if len(set(output_paths)) != len(output_paths):
        raise BreakpointError("Output, exclusions, and manifest paths must be distinct.")
    if set(resolved_paths).intersection(output_paths):
        raise BreakpointError("An output path must not overwrite an input or frozen artifact.")

    ast_fields, ast_rows = _read_csv(ast_path)
    source_ast_row_count = len(ast_rows)
    _require_columns(ast_fields, AST_REQUIRED_COLUMNS, "Normalized AST CSV")
    normalization_manifest = verify_normalization_manifest(
        normalization_manifest_path,
        ast_path,
        ast_rows,
        allow_unverified=allow_unverified_normalization,
    )
    collisions = sorted(set(ast_fields).intersection(ADDED_OUTPUT_COLUMNS + EXCLUSION_COLUMNS))
    if collisions:
        raise BreakpointError(f"Normalized AST already contains interpreter-owned columns: {collisions}")

    context_manifest = None
    context_excluded: list[dict[str, str]] = []
    if context_path is not None:
        _require_columns(
            ast_fields,
            ("isolate_id", "biosample_accession", "assembly_accession"),
            "Normalized AST CSV used with context metadata",
        )
        context_fields, context_rows = _read_csv(context_path)
        _require_columns(context_fields, CONTEXT_REQUIRED_COLUMNS, "Breakpoint context metadata")
        collisions = sorted(set(CONTEXT_ADDED_COLUMNS).intersection(ast_fields))
        if collisions:
            raise BreakpointError(
                f"Normalized AST already contains context-owned columns: {collisions}"
            )
        ast_rows, context_excluded = enrich_with_context(ast_rows, context_rows)
        ast_fields = ast_fields + list(CONTEXT_ADDED_COLUMNS)
        context_manifest = {
            "path": str(context_path),
            "sha256": sha256_file(context_path),
            "rows": len(context_rows),
        }

    rules = load_breakpoint_rules(rules_path, artifact_path)
    interpreted, interpretation_excluded = interpret_rows(ast_rows, rules)
    excluded = context_excluded + interpretation_excluded
    output_fields = ast_fields + list(ADDED_OUTPUT_COLUMNS)
    exclusion_fields = ast_fields + list(EXCLUSION_COLUMNS)
    _write_csv(output_path, output_fields, interpreted)
    _write_csv(exclusions_path, exclusion_fields, excluded)

    reason_counts: dict[str, int] = {}
    for row in excluded:
        reason = row["exclusion_reason"]
        reason_counts[reason] = reason_counts.get(reason, 0) + 1
    category_counts = {
        category: sum(row["ast_category"] == category for row in interpreted)
        for category in ("S", "I", "R")
    }
    concordance_counts: dict[str, int] = {}
    for row in interpreted:
        key = row["submitted_category_concordance"]
        concordance_counts[key] = concordance_counts.get(key, 0) + 1

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_status": (
            "INDEPENDENT_BREAKPOINT_INTERPRETATION"
            if normalization_manifest["scientific_status"] == "SOURCE_NORMALIZATION_ONLY"
            else "SMOKE_TEST_ONLY"
        ),
        "algorithm": {
            "name": "censoring-aware-interval-interpretation",
            "version": "1.0.0",
            "mic": "S <= susceptible_breakpoint; I between thresholds; R > resistant_breakpoint",
            "zone": "R < resistant_breakpoint; I between thresholds; S >= susceptible_breakpoint",
            "censoring_policy": "interpret only when the full measurement interval implies one category",
            "constraint_policy": "exact case-insensitive whitespace-normalized matching; no inference or unit conversion",
        },
        "sources": {
            "normalized_ast": {
                "path": str(ast_path),
                "sha256": sha256_file(ast_path),
                "rows": source_ast_row_count,
            },
            "normalization_manifest": normalization_manifest,
            "breakpoint_context_metadata": context_manifest,
            "breakpoint_rules": {
                "path": str(rules_path),
                "sha256": sha256_file(rules_path),
                "rows": len(rules),
                "rule_row_sha256": {rule.rule_id: rule.row_sha256 for rule in rules},
            },
            "frozen_breakpoint_artifact": {
                "path": str(artifact_path),
                "sha256": sha256_file(artifact_path),
                "bytes": artifact_path.stat().st_size,
            },
        },
        "outputs": {
            "interpreted": {
                "path": str(output_path),
                "sha256": sha256_file(output_path),
                "rows": len(interpreted),
            },
            "exclusions": {
                "path": str(exclusions_path),
                "sha256": sha256_file(exclusions_path),
                "rows": len(excluded),
            },
        },
        "category_counts": category_counts,
        "submitted_category_concordance_counts": dict(sorted(concordance_counts.items())),
        "exclusion_counts": dict(sorted(reason_counts.items())),
        "limitations": [
            "Rule transcription is not validated here against the scientific content of the source artifact.",
            "Rows with absent comparators or rule-constrained metadata are excluded rather than imputed.",
            "Rows that fail authenticated source-stage assay-metadata eligibility are excluded before rule matching.",
            "No measurement-unit conversion, method mapping, site mapping, or duplicate-test reconciliation is performed.",
        ],
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_manifest = manifest_path.with_name(manifest_path.name + ".tmp")
    temporary_manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary_manifest.replace(manifest_path)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ast", type=Path, required=True, help="Normalized quantitative AST CSV")
    parser.add_argument(
        "--normalization-manifest",
        type=Path,
        required=True,
        help="Authenticated manifest produced with the normalized AST input",
    )
    parser.add_argument(
        "--allow-unverified-normalization",
        action="store_true",
        help="Allow a synthetic SMOKE_TEST_ONLY normalization manifest; never use for research data.",
    )
    parser.add_argument("--rules", type=Path, required=True, help="Machine-readable breakpoint rules CSV")
    parser.add_argument(
        "--breakpoint-artifact",
        type=Path,
        required=True,
        help="Original frozen breakpoint table or document referenced by every rule",
    )
    parser.add_argument(
        "--context-metadata",
        type=Path,
        help="Exact isolate metadata supplying clinical indication before rule selection",
    )
    parser.add_argument("--output", type=Path, default=Path("data/interim/ast_interpreted.csv"))
    parser.add_argument(
        "--exclusions", type=Path, default=Path("data/interim/breakpoint_exclusions.csv")
    )
    parser.add_argument(
        "--manifest", type=Path, default=Path("data/interim/breakpoint_manifest.json")
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        manifest = run_interpretation(
            ast_path=args.ast,
            normalization_manifest_path=args.normalization_manifest,
            rules_path=args.rules,
            artifact_path=args.breakpoint_artifact,
            output_path=args.output,
            exclusions_path=args.exclusions,
            manifest_path=args.manifest,
            context_path=args.context_metadata,
            allow_unverified_normalization=args.allow_unverified_normalization,
        )
    except BreakpointError as exc:
        raise SystemExit(f"Breakpoint interpretation failed: {exc}") from exc
    counts = manifest["outputs"]
    print(
        "Breakpoint interpretation complete: "
        f"{counts['interpreted']['rows']} interpreted, "
        f"{counts['exclusions']['rows']} excluded. Manifest: {args.manifest}"
    )


if __name__ == "__main__":
    main()
