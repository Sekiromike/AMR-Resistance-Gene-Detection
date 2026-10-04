"""Construct the modeling cohort from interpreted AST and frozen isolate metadata.

The interpreted AST table is the sole source of phenotype rows. Isolate-level
metadata may enrich those rows, but it may not replace phenotype provenance or
repair identity mismatches. Duplicate isolate-antibiotic measurements are
excluded as unresolved groups pending a separately prespecified reconciliation
policy.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


# Fields every endpoint needs: identity, source provenance, and the
# quantitative measurement itself.
AST_BASE_REQUIRED = (
    "source_dataset",
    "source_release",
    "source_fingerprint",
    "source_row_number",
    "source_ast_record_id",
    "isolate_id",
    "biosample_accession",
    "assembly_accession",
    "scientific_name",
    "antibiotic",
    "ast_measurement_type",
    "measurement_sign",
    "raw_ast_value",
    "ast_value",
    "ast_unit",
    "collection_date",
    "country",
    "bioproject_accession",
)

# Fields required only to interpret a measurement into a clinical category.
# docs/ENDPOINT_AMENDMENT.md records why no public source supplies these.
AST_CATEGORICAL_REQUIRED = (
    "submitted_ast_category",
    "ast_method",
    "ast_testing_date",
    "clinical_indication",
    "ast_category",
    "breakpoint_rule_id",
    "breakpoint_rule_row_sha256",
    "breakpoint_artifact_sha256",
    "breakpoint_standard",
    "breakpoint_version",
)

AST_REQUIRED = AST_BASE_REQUIRED + AST_CATEGORICAL_REQUIRED

METADATA_BASE_REQUIRED = (
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
)

# Absent from every public source audited; see docs/ENDPOINT_AMENDMENT.md.
METADATA_CATEGORICAL_REQUIRED = (
    "site",
    "surveillance_network",
    "clinical_indication",
)

METADATA_REQUIRED = METADATA_BASE_REQUIRED + METADATA_CATEGORICAL_REQUIRED

ENDPOINT_PROFILES = {
    "mic_regression": (AST_BASE_REQUIRED, METADATA_BASE_REQUIRED),
    "categorical_sir": (AST_REQUIRED, METADATA_REQUIRED),
}
# Fail closed: an unspecified profile keeps the strict historical contract.
DEFAULT_ENDPOINT_PROFILE = "categorical_sir"

METADATA_ADDED = tuple(field for field in METADATA_REQUIRED if field not in AST_REQUIRED)
INTERPRETER_OWNED = {
    "ast_category",
    "breakpoint_rule_id",
    "breakpoint_rule_row_sha256",
    "breakpoint_artifact_sha256",
    "breakpoint_standard",
    "breakpoint_version",
}
IDENTITY_FIELDS = ("biosample_accession", "assembly_accession")
# A genome assembled here from raw reads has no registered assembly; its
# provenance is the run accession plus genome SHA-256 (validated downstream).
# Mirrors BLANK_ALLOWED_FOR_RAW_READS in validate_research_cohort.py.
BLANK_ALLOWED_FOR_RAW_READS = frozenset({"assembly_accession"})
# Absent in source, kept for the MIC endpoint and never inferred; the analyses
# that need them (specimen-type subgroups) exclude those rows. Mirrors the
# validator's ABSENT_ALLOWED_FOR_MIC for the metadata fields.
ABSENT_ALLOWED_FOR_MIC = frozenset({"specimen_source"})
CONTEXT_FIELDS = ("clinical_indication",)
MIC_CONTEXT_FIELDS: tuple[str, ...] = ()
EXCLUSION_FIELDS = ("cohort_exclusion_reason", "cohort_exclusion_detail")
GENOME_REQUIRED = (
    "isolate_id",
    "biosample_accession",
    "assembly_accession",
    "eligible_for_modeling",
    "overall_status",
    "exclusion_reasons",
)


class CohortConstructionError(ValueError):
    """Raised when the source tables violate a fail-closed cohort contract."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _blank(value: object) -> bool:
    return value is None or str(value).strip().lower() in {"", "nan", "none", "null"}


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise CohortConstructionError(f"Input does not exist: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise CohortConstructionError(f"CSV has no header: {path}")
        fields = list(reader.fieldnames)
        if len(fields) != len(set(fields)):
            raise CohortConstructionError(f"CSV has duplicate column names: {path}")
        return fields, [dict(row) for row in reader]


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise CohortConstructionError(f"Input does not exist: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise CohortConstructionError(f"TSV has no header: {path}")
        fields = list(reader.fieldnames)
        if len(fields) != len(set(fields)):
            raise CohortConstructionError(f"TSV has duplicate column names: {path}")
        return fields, [dict(row) for row in reader]


def _require_columns(fields: list[str], required: tuple[str, ...], label: str) -> None:
    missing = sorted(set(required).difference(fields))
    if missing:
        raise CohortConstructionError(f"{label} is missing required columns: {missing}")


def _write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _unique_isolate_index(
    rows: list[dict[str, str]], label: str
) -> dict[str, dict[str, str]]:
    by_isolate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_isolate[str(row.get("isolate_id", "")).strip()].append(row)
    missing_ids = len(by_isolate.get("", []))
    if missing_ids:
        raise CohortConstructionError(f"{label} has {missing_ids} rows with blank isolate_id.")
    duplicate_ids = sorted(key for key, values in by_isolate.items() if len(values) != 1)
    if duplicate_ids:
        raise CohortConstructionError(
            f"{label} isolate_id is not unique: {duplicate_ids[:20]}"
        )

    for identity in IDENTITY_FIELDS:
        reverse: dict[str, set[str]] = defaultdict(set)
        for isolate_id, values in by_isolate.items():
            value = str(values[0].get(identity, "")).strip()
            if value:
                reverse[value].add(isolate_id)
        ambiguous = sorted(value for value, isolates in reverse.items() if len(isolates) > 1)
        if ambiguous:
            raise CohortConstructionError(
                f"{label} {identity} maps to multiple isolates: {ambiguous[:20]}"
            )
    return {key: values[0] for key, values in by_isolate.items()}


# Amendment 005, rule R: repeated MIC measurements for one isolate and drug
# are reconciled by interval intersection -- the most specific statement
# consistent with every submitted measurement. Nothing is chosen.
RECONCILIATION_FIELDS = ("reconciled_from", "reconciliation")
_LOWER_SIGNS = {">=": True, ">": False}   # value -> inclusive?
_UPPER_SIGNS = {"<=": True, "<": False}
_POINT_SIGNS = {"=", "=="}


def _mic_interval(row: dict[str, str]) -> tuple[tuple[float, bool, str] | None, tuple[float, bool, str] | None]:
    """Return (lower, upper) bounds as (value, inclusive, original text)."""
    sign = str(row.get("measurement_sign", "")).strip()
    text = str(row.get("ast_value", "")).strip()
    try:
        value = float(text)
    except ValueError as error:
        raise CohortConstructionError(f"Non-numeric MIC value in repeat group: {text!r}") from error
    if sign in _POINT_SIGNS:
        return (value, True, text), (value, True, text)
    if sign in _LOWER_SIGNS:
        return (value, _LOWER_SIGNS[sign], text), None
    if sign in _UPPER_SIGNS:
        return None, (value, _UPPER_SIGNS[sign], text)
    raise CohortConstructionError(f"Unsupported MIC comparator in repeat group: {sign!r}")


def reconcile_mic_repeats(rows: list[dict[str, str]]) -> tuple[dict[str, str] | None, str]:
    """Intersect the MIC intervals of one isolate-drug group.

    Returns (reconciled_row, "") on success or (None, exclusion_reason).
    """
    units = {str(row.get("ast_unit", "")).strip() for row in rows}
    if len(units) != 1:
        return None, "repeat_unit_conflict"
    lower: tuple[float, bool, str] | None = None
    upper: tuple[float, bool, str] | None = None
    for row in rows:
        row_lower, row_upper = _mic_interval(row)
        if row_lower is not None:
            # Tighter lower bound: larger value; at a tie, exclusive beats inclusive.
            if lower is None or row_lower[0] > lower[0] or (
                row_lower[0] == lower[0] and not row_lower[1]
            ):
                lower = row_lower
        if row_upper is not None:
            if upper is None or row_upper[0] < upper[0] or (
                row_upper[0] == upper[0] and not row_upper[1]
            ):
                upper = row_upper
    if lower is not None and upper is not None:
        if lower[0] > upper[0] or (lower[0] == upper[0] and not (lower[1] and upper[1])):
            return None, "discordant_repeat_measurements"
        if lower[0] == upper[0]:
            sign, text = "=", lower[2]
        else:
            return None, "repeat_not_representable"
    elif lower is not None:
        sign, text = (">=" if lower[1] else ">"), lower[2]
    elif upper is not None:
        sign, text = ("<=" if upper[1] else "<"), upper[2]
    else:  # pragma: no cover - every row contributes a bound
        return None, "repeat_not_representable"

    def row_number(row: dict[str, str]) -> int:
        try:
            return int(str(row.get("source_row_number", "")).strip())
        except ValueError:
            return 0

    carrier = dict(min(rows, key=row_number))
    carrier["measurement_sign"] = sign
    carrier["ast_value"] = text
    carrier["raw_ast_value"] = ";".join(
        f"{str(r.get('measurement_sign', '')).strip()}{str(r.get('ast_value', '')).strip()}"
        for r in sorted(rows, key=row_number)
    )
    carrier["reconciled_from"] = ";".join(
        sorted(str(r.get("source_ast_record_id", "")).strip() for r in rows)
    )
    carrier["reconciliation"] = "interval_intersection"
    return carrier, ""


def _is_mic_endpoint_row(row: dict[str, str]) -> bool:
    if str(row.get("ast_measurement_type", "")).strip().upper() != "MIC":
        return False
    eligible = str(row.get("mic_eligible", "")).strip().lower()
    return eligible in {"", "true", "1"}


# docs/MIC_EVALUATION_SPEC.md: a value more than 0.1 log2 from a doubling
# dilution is not on the MIC scale (must equal evaluate_mic_predictions.SNAP_TOLERANCE).
DOUBLING_SCALE_TOLERANCE = 0.1


def _on_doubling_scale(value: str) -> bool:
    try:
        number = float(str(value).strip())
    except ValueError:
        return False
    if not (math.isfinite(number) and number > 0):
        return False
    log2 = math.log2(number)
    return abs(log2 - round(log2)) <= DOUBLING_SCALE_TOLERANCE


def prepare_mic_rows(
    ast_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Filter to MIC-endpoint rows (R1) and reconcile repeats (R2-R4)."""
    kept: list[dict[str, str]] = []
    excluded: list[dict[str, str]] = []
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for source in ast_rows:
        row = dict(source)
        if not _is_mic_endpoint_row(row):
            row["cohort_exclusion_reason"] = "not_mic_endpoint_measurement"
            row["cohort_exclusion_detail"] = "Only MIC-eligible rows enter the MIC endpoint."
            excluded.append(row)
            continue
        if not _on_doubling_scale(row.get("ast_value", "")):
            row["cohort_exclusion_reason"] = "mic_value_off_doubling_scale"
            row["cohort_exclusion_detail"] = (
                f"{row.get('ast_value', '')} mg/L is more than {DOUBLING_SCALE_TOLERANCE} log2 "
                "from a doubling dilution (docs/MIC_EVALUATION_SPEC.md)."
            )
            excluded.append(row)
            continue
        key = (str(row.get("isolate_id", "")).strip(), str(row.get("antibiotic", "")).strip())
        groups[key].append(row)
    for key in sorted(groups):
        rows = groups[key]
        if len(rows) == 1:
            row = rows[0]
            row["reconciled_from"] = str(row.get("source_ast_record_id", "")).strip()
            row["reconciliation"] = "single"
            kept.append(row)
            continue
        reconciled, reason = reconcile_mic_repeats(rows)
        if reconciled is None:
            for row in rows:
                row["cohort_exclusion_reason"] = reason
                row["cohort_exclusion_detail"] = (
                    f"{len(rows)} repeated MIC rows could not be reconciled by interval intersection."
                )
                excluded.append(row)
        else:
            kept.append(reconciled)
    return kept, excluded


def construct_rows(
    ast_rows: list[dict[str, str]],
    metadata_rows: list[dict[str, str]],
    genome_rows: list[dict[str, str]] | None = None,
    *,
    profile: str = DEFAULT_ENDPOINT_PROFILE,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    if profile not in ENDPOINT_PROFILES:
        raise CohortConstructionError(
            f"Unknown endpoint profile {profile!r}; expected one of {sorted(ENDPOINT_PROFILES)}"
        )
    _, metadata_required = ENDPOINT_PROFILES[profile]
    context_fields = CONTEXT_FIELDS if profile == "categorical_sir" else MIC_CONTEXT_FIELDS
    metadata = _unique_isolate_index(metadata_rows, "Metadata")
    genomes = _unique_isolate_index(genome_rows, "Genome manifest") if genome_rows is not None else None
    accepted: list[dict[str, str]] = []
    excluded: list[dict[str, str]] = []
    ast_required, _ = ENDPOINT_PROFILES[profile]
    metadata_added = tuple(field for field in metadata_required if field not in ast_required)
    if profile == "mic_regression":
        ast_rows, pre_excluded = prepare_mic_rows(ast_rows)
        excluded.extend(pre_excluded)

    key_counts = Counter(
        (str(row.get("isolate_id", "")).strip(), str(row.get("antibiotic", "")).strip())
        for row in ast_rows
    )
    for source in ast_rows:
        row = dict(source)
        isolate_id = str(row.get("isolate_id", "")).strip()
        antibiotic = str(row.get("antibiotic", "")).strip()
        reason = ""
        detail = ""
        if key_counts[(isolate_id, antibiotic)] != 1:
            reason = "duplicate_isolate_antibiotic_unreconciled"
            detail = "Every source record in the duplicate group was excluded."
        elif isolate_id not in metadata:
            reason = "isolate_metadata_not_found"
            detail = "No exact isolate_id row exists in the frozen metadata table."
        else:
            enrichment = metadata[isolate_id]
            mismatches = [
                field
                for field in IDENTITY_FIELDS + context_fields
                if str(row.get(field, "")).strip() != str(enrichment.get(field, "")).strip()
            ]
            if mismatches:
                reason = "identity_mismatch"
                detail = "AST and metadata disagree for: " + ",".join(mismatches)
            else:
                raw_reads = str(enrichment.get("genome_source", "")).strip() == "raw_reads"
                missing = [
                    field
                    for field in metadata_required
                    if _blank(enrichment.get(field, ""))
                    and not (raw_reads and field in BLANK_ALLOWED_FOR_RAW_READS)
                    and not (profile == "mic_regression" and field in ABSENT_ALLOWED_FOR_MIC)
                ]
                if missing:
                    reason = "incomplete_isolate_metadata"
                    detail = "Required metadata are absent: " + ",".join(sorted(missing))
                elif genomes is not None:
                    genome = genomes.get(isolate_id)
                    if genome is None:
                        reason = "genome_analysis_not_found"
                        detail = "No executed genome-manifest row exists for this isolate."
                    else:
                        genome_mismatches = [
                            field
                            for field in IDENTITY_FIELDS
                            if str(row.get(field, "")).strip()
                            != str(genome.get(field, "")).strip()
                        ]
                        if genome_mismatches:
                            reason = "genome_identity_mismatch"
                            detail = "AST and genome manifest disagree for: " + ",".join(
                                genome_mismatches
                            )
                        elif str(genome.get("eligible_for_modeling", "")).strip().lower() != "true":
                            reason = "genome_not_eligible"
                            detail = str(genome.get("exclusion_reasons", "")).strip()

        if reason:
            row["cohort_exclusion_reason"] = reason
            row["cohort_exclusion_detail"] = detail
            excluded.append(row)
            continue

        enrichment = metadata[isolate_id]
        for field in metadata_added:
            row[field] = str(enrichment[field]).strip()
        if genomes is not None:
            row["genome_qc_status"] = str(genomes[isolate_id]["overall_status"]).strip()
        if profile == "categorical_sir":
            row["breakpoint_rule_sha256"] = row["breakpoint_rule_row_sha256"]
        accepted.append(row)
    return accepted, excluded


def run_construction(
    interpreted_ast: Path | list[Path],
    isolate_metadata: Path,
    output: Path,
    exclusions: Path,
    manifest: Path,
    genome_manifest: Path | None = None,
    *,
    profile: str = DEFAULT_ENDPOINT_PROFILE,
) -> dict[str, object]:
    # Several AST sources (development and external) may be combined into one
    # cohort table; each keeps its own provenance and hash in the manifest.
    ast_paths = list(interpreted_ast) if isinstance(interpreted_ast, (list, tuple)) else [interpreted_ast]
    if not ast_paths:
        raise CohortConstructionError("At least one AST input is required.")
    input_paths = {path.resolve() for path in ast_paths} | {isolate_metadata.resolve()}
    if len(input_paths) != len(ast_paths) + 1:
        raise CohortConstructionError("AST inputs and isolate metadata must be distinct files.")
    if genome_manifest is not None:
        input_paths.add(genome_manifest.resolve())
    output_paths = {output.resolve(), exclusions.resolve(), manifest.resolve()}
    if len(output_paths) != 3 or input_paths.intersection(output_paths):
        raise CohortConstructionError("Inputs and the three output paths must be distinct.")

    ast_fields: list[str] = []
    ast_rows: list[dict[str, str]] = []
    ast_sources: list[dict[str, object]] = []
    per_source_fields: list[tuple[Path, list[str]]] = []
    for path in ast_paths:
        fields, rows = _read_csv(path)
        per_source_fields.append((path, fields))
        ast_fields += [field for field in fields if field not in ast_fields]
        ast_rows.extend(rows)
        ast_sources.append({"path": str(path), "sha256": sha256_file(path), "rows": len(rows)})
    metadata_fields, metadata_rows = _read_csv(isolate_metadata)
    genome_rows = None
    genome_source = None
    if genome_manifest is not None:
        genome_fields, genome_rows = _read_tsv(genome_manifest)
        _require_columns(genome_fields, GENOME_REQUIRED, "Genome manifest")
        genome_source = {
            "path": str(genome_manifest),
            "sha256": sha256_file(genome_manifest),
            "rows": len(genome_rows),
        }
    if profile not in ENDPOINT_PROFILES:
        raise CohortConstructionError(
            f"Unknown endpoint profile {profile!r}; expected one of {sorted(ENDPOINT_PROFILES)}"
        )
    ast_required, metadata_required = ENDPOINT_PROFILES[profile]
    label = "Interpreted AST" if profile == "categorical_sir" else "Normalized AST"
    for path, fields in per_source_fields:
        _require_columns(fields, ast_required, f"{label} ({path.name})")
    _require_columns(metadata_fields, metadata_required, "Isolate metadata")
    metadata_added = tuple(field for field in metadata_required if field not in ast_required)
    collisions = sorted(INTERPRETER_OWNED.intersection(metadata_added))
    if collisions:
        raise CohortConstructionError(f"Metadata may not own interpreted AST fields: {collisions}")

    accepted, excluded = construct_rows(ast_rows, metadata_rows, genome_rows, profile=profile)
    output_fields = ast_fields + [field for field in metadata_added if field not in ast_fields]
    if profile == "categorical_sir":
        output_fields += ["breakpoint_rule_sha256"]
    else:
        output_fields += [field for field in RECONCILIATION_FIELDS if field not in output_fields]
    exclusion_fields = ast_fields + list(EXCLUSION_FIELDS)
    _write_csv(output, output_fields, accepted)
    _write_csv(exclusions, exclusion_fields, excluded)

    exclusion_counts = Counter(row["cohort_exclusion_reason"] for row in excluded)
    payload: dict[str, object] = {
        "schema_version": "1.0.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_status": "SOURCE_LINKED_COHORT_CONSTRUCTION",
        "endpoint_profile": profile,
        "policy": {
            "phenotype_source": "interpreted AST rows only",
            "identity_join": "exact isolate_id, BioSample, and assembly accession",
            "duplicates": (
                "MIC: interval-intersection reconciliation (amendment 005)"
                if profile == "mic_regression"
                else "exclude every unresolved isolate-antibiotic duplicate"
            ),
            "missing_metadata": "exclude; no imputation",
        },
        "sources": {
            "interpreted_ast": ast_sources[0],
            "ast_sources": ast_sources,
            "isolate_metadata": {
                "path": str(isolate_metadata),
                "sha256": sha256_file(isolate_metadata),
                "rows": len(metadata_rows),
            },
            "executed_genome_manifest": genome_source,
        },
        "outputs": {
            "cohort": {"path": str(output), "sha256": sha256_file(output), "rows": len(accepted)},
            "exclusions": {
                "path": str(exclusions),
                "sha256": sha256_file(exclusions),
                "rows": len(excluded),
            },
        },
        "exclusion_counts": dict(sorted(exclusion_counts.items())),
        "limitations": [
            "The isolate metadata table must itself be produced and frozen by genome and epidemiology workflows.",
            "No duplicate AST reconciliation, metadata imputation, or identity repair is performed.",
        ],
    }
    manifest.parent.mkdir(parents=True, exist_ok=True)
    temporary = manifest.with_name(manifest.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(manifest)
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--interpreted-ast",
        type=Path,
        action="append",
        required=True,
        help="AST table; repeat to combine development and external sources.",
    )
    parser.add_argument("--isolate-metadata", type=Path, required=True)
    parser.add_argument("--genome-manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--exclusions", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument(
        "--endpoint-profile",
        choices=sorted(ENDPOINT_PROFILES),
        default=DEFAULT_ENDPOINT_PROFILE,
        help="Endpoint contract to build against. See docs/ENDPOINT_AMENDMENT.md.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        payload = run_construction(
            args.interpreted_ast,
            args.isolate_metadata,
            args.output,
            args.exclusions,
            args.manifest,
            args.genome_manifest,
            profile=args.endpoint_profile,
        )
    except CohortConstructionError as exc:
        raise SystemExit(f"Cohort construction failed: {exc}") from exc
    print(
        f"Cohort construction complete: {payload['outputs']['cohort']['rows']} included, "
        f"{payload['outputs']['exclusions']['rows']} excluded."
    )


if __name__ == "__main__":
    main()
