"""Normalize an NCBI Pathogen Detection AST export without inventing metadata.

The AST Browser and BigQuery exports are submitter-provided source data. This
script preserves the submitted phenotype, quantitative measurement, source
record identity, and raw contextual fields. It emits structural exclusions and
separately marks records that lack source-stage assay metadata needed before
breakpoint interpretation. It does not assign a modern breakpoint version or
claim that the submitted interpretation is correct.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


NORMALIZATION_SCHEMA_VERSION = "2.0.0"
DEFAULT_ACQUISITION_MANIFEST = "ncbi_acquisition_manifest.json"
DEFAULT_ACQUISITION_AMENDMENT = "ncbi_acquisition_manifest_amendment.json"
NULL_LIKE_STRINGS = frozenset({"", "nan", "none", "null", "<na>"})
MIC_AST_METHODS = frozenset(
    {
        "agar dilution",
        "broth macrodilution",
        "broth microdilution",
        "gradient diffusion",
    }
)
ZONE_AST_METHODS = frozenset({"disk diffusion"})
DISK_CONTENT_PATTERN = re.compile(
    r"^\d+(?:\.\d+)?(?:\s*/\s*\d+(?:\.\d+)?)?\s*(?:ug|mcg|µg|μg)$",
    flags=re.IGNORECASE,
)

ALIASES = {
    "source_native_record_id": ("id", "source_ast_record_id"),
    "source_record_checksum": ("checksum", "source_record_checksum"),
    "isolate_id": ("isolate_id", "isolate", "target_acc", "target_accession"),
    "biosample_accession": ("biosample_accession", "biosample_acc", "biosample"),
    "assembly_accession": ("assembly_accession", "assembly_acc", "asm_acc", "assembly"),
    "taxgroup_name": ("taxgroup_name",),
    "scientific_name": ("scientific_name", "organism", "organism_name"),
    "antibiotic": ("antibiotic", "antimicrobial", "drug"),
    "submitted_ast_category": (
        "submitted_ast_category",
        "resistance_phenotype",
        "phenotype",
        "ast_phenotype",
    ),
    "mic": ("mic", "mic_mg_l", "mic_mg_per_l"),
    "mic_secondary": ("mic_secondary", "secondary_mic"),
    "disk_diffusion": (
        "disk_diffusion",
        "disk_diffusion_mm",
        "zone_diameter",
        "zone_diameter_mm",
    ),
    "disk_diffusion_secondary": (
        "disk_diffusion_secondary",
        "secondary_disk_diffusion",
        "secondary_zone_diameter",
    ),
    "measurement_sign": ("measurement_sign", "sign", "comparator"),
    "ast_method": (
        "ast_method",
        "laboratory_typing_method",
        "laboratory_method",
        "method",
    ),
    "submitted_standard": ("submitted_standard", "testing_standard", "standard"),
    "bioproject_accession": ("bioproject_accession", "bioproject_acc", "bioproject"),
    "collection_date": ("collection_date", "collection_dt"),
    # geo_loc_name is kept verbatim. Country is derived from its INSDC prefix
    # only when that prefix is in the frozen INSDC vocabulary (amendment 006).
    "geo_loc_name": ("geo_loc_name", "geographic_location"),
    "country": ("country",),
    "epi_type": ("epi_type",),
    "isolation_type": ("isolation_type",),
    "isolation_source": ("isolation_source",),
    "host": ("host",),
    "ast_record_creation_date": ("ast_record_creation_date", "creation_date"),
    "target_creation_date": ("target_creation_date",),
    "ast_platform": ("ast_platform", "platform", "laboratory_typing_platform"),
    "ast_vendor": ("ast_vendor", "vendor"),
    "ast_reagent": ("ast_reagent", "reagent", "testing_reagent"),
    "disk_content": ("disk_content", "disk_content_ug", "disk_potency"),
    "ast_testing_date": ("ast_testing_date", "testing_date", "test_date"),
}

CATEGORY_MAP = {
    "susceptible": "S",
    "s": "S",
    "intermediate": "I",
    "i": "I",
    "resistant": "R",
    "r": "R",
    "nonsusceptible": "NS",
    "non_susceptible": "NS",
    "non-susceptible": "NS",
    "ns": "NS",
    "susceptible-dose dependent": "SDD",
    "susceptible_dose_dependent": "SDD",
    "sdd": "SDD",
}


class NormalizationError(ValueError):
    """Raised when source identity or provenance cannot be normalized safely."""


def canonical_column(name: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower())).strip("_")


def canonicalize_frame(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    canonical = frame.copy()
    canonical.columns = [canonical_column(column) for column in canonical.columns]
    duplicates = canonical.columns[canonical.columns.duplicated()].tolist()
    if duplicates:
        raise NormalizationError(
            f"{label} has duplicate columns after canonicalization: {sorted(set(duplicates))}"
        )
    return canonical.reset_index(drop=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_table(path: Path) -> pd.DataFrame:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(8192)
    try:
        delimiter = csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
    except csv.Error:
        delimiter = "\t" if "\t" in sample else ","
    frame = pd.read_csv(path, sep=delimiter, dtype="string", keep_default_na=True)
    return canonicalize_frame(frame, str(path))


def blank_mask(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.strip()
    return text.isna() | text.str.lower().isin(NULL_LIKE_STRINGS)


def clean_string_series(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.strip()
    return text.mask(blank_mask(text), pd.NA)


def first_column(frame: pd.DataFrame, logical_name: str) -> pd.Series:
    result = pd.Series(pd.NA, index=frame.index, dtype="string")
    observed_aliases: list[str] = []
    for candidate in ALIASES[logical_name]:
        if candidate not in frame.columns:
            continue
        observed_aliases.append(candidate)
        values = clean_string_series(frame[candidate])
        overlap = ~blank_mask(result) & ~blank_mask(values)
        conflict = overlap & result.map(normalized_metadata_value).ne(
            values.map(normalized_metadata_value)
        )
        if conflict.any():
            raise NormalizationError(
                f"Conflicting aliases for {logical_name} ({observed_aliases}) "
                f"in {int(conflict.sum())} rows"
            )
        result = coalesce_columns(result, values)
    return result


def coalesce_columns(left: pd.Series, right: pd.Series) -> pd.Series:
    return clean_string_series(left).mask(blank_mask(left), clean_string_series(right))


def normalize_category(value: object) -> str | pd.NA:
    if pd.isna(value):
        return pd.NA
    key = str(value).strip().lower()
    if key in NULL_LIKE_STRINGS:
        return pd.NA
    return CATEGORY_MAP.get(key, key.upper())


def normalized_metadata_value(value: object) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().casefold().split())


def valid_ast_method(value: object, measurement_type: object) -> bool:
    method = normalized_metadata_value(value)
    kind = "" if pd.isna(measurement_type) else str(measurement_type).strip().upper()
    if kind == "MIC":
        return method in MIC_AST_METHODS
    if kind == "ZONE":
        return method in ZONE_AST_METHODS
    return False


def valid_disk_content(value: object) -> bool:
    if pd.isna(value):
        return False
    text = str(value).strip()
    if not DISK_CONTENT_PATTERN.fullmatch(text):
        return False
    numeric = re.findall(r"\d+(?:\.\d+)?", text)
    return bool(numeric) and all(float(part) > 0 for part in numeric)


def valid_iso_date_or_datetime(value: object) -> bool:
    if pd.isna(value):
        return False
    text = str(value).strip()
    try:
        if re.fullmatch(r"\d{4}-\d{2}-\d{2}", text):
            datetime.fromisoformat(text)
        else:
            datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return False
    return bool(re.match(r"^\d{4}-\d{2}-\d{2}(?:T|$)", text))


def canonical_raw_record_id(row: pd.Series, columns: Iterable[str]) -> str:
    fields = []
    for column in sorted(columns):
        value = row.get(column, pd.NA)
        fields.append([column, None if pd.isna(value) else str(value)])
    material = json.dumps(fields, ensure_ascii=False, separators=(",", ":"))
    return "ncbi-ast-sha256-" + hashlib.sha256(material.encode("utf-8")).hexdigest()


def _assert_unique_isolate_identifiers(frame: pd.DataFrame) -> None:
    for field in ("isolate_id", "biosample_accession", "assembly_accession"):
        present = ~blank_mask(frame[field])
        duplicated = present & frame[field].duplicated(keep=False)
        if duplicated.any():
            examples = sorted(frame.loc[duplicated, field].dropna().astype(str).unique())[:5]
            raise NormalizationError(
                f"Isolate metadata is ambiguous and not one-to-one: "
                f"duplicate {field} values {examples}"
            )


def _assert_consistent_identity_mapping(frame: pd.DataFrame, label: str) -> None:
    pairs = (
        ("isolate_id", "biosample_accession"),
        ("isolate_id", "assembly_accession"),
        ("biosample_accession", "isolate_id"),
    )
    for left, right in pairs:
        present = ~blank_mask(frame[left]) & ~blank_mask(frame[right])
        if not present.any():
            continue
        counts = frame.loc[present].groupby(left, dropna=False)[right].nunique(dropna=True)
        conflicting = counts[counts > 1]
        if not conflicting.empty:
            examples = [str(value) for value in conflicting.index[:5]]
            raise NormalizationError(
                f"{label} has conflicting {left}-to-{right} identity mappings for {examples}"
            )


def _merge_isolate_metadata(out: pd.DataFrame, isolates: pd.DataFrame) -> pd.DataFrame:
    isolate_fields = (
        "isolate_id",
        "biosample_accession",
        "assembly_accession",
        "scientific_name",
        "bioproject_accession",
        "collection_date",
        "geo_loc_name",
        "country",
        "epi_type",
        "isolation_type",
        "isolation_source",
        "host",
        "target_creation_date",
    )
    isolate_lookup = pd.DataFrame(index=isolates.index)
    for field in isolate_fields:
        isolate_lookup[field] = first_column(isolates, field)

    _assert_unique_isolate_identifiers(isolate_lookup)
    _assert_consistent_identity_mapping(isolate_lookup, "Isolate metadata")
    candidates = [
        field
        for field in ("isolate_id", "biosample_accession")
        if (~blank_mask(out[field])).any() and (~blank_mask(isolate_lookup[field])).any()
    ]
    if not candidates:
        raise NormalizationError(
            "AST and isolate metadata have no populated isolate_id or biosample_accession join key"
        )
    join_key = candidates[0]
    if blank_mask(isolate_lookup[join_key]).any():
        raise NormalizationError(f"Isolate metadata contains blank {join_key} values")

    renamed = isolate_lookup.rename(
        columns={field: f"__isolate_{field}" for field in isolate_fields if field != join_key}
    )
    merged = out[[join_key]].merge(
        renamed,
        on=join_key,
        how="left",
        sort=False,
        validate="many_to_one",
        indicator="__isolate_merge",
    )
    merged.index = out.index
    if len(merged) != len(out):
        raise NormalizationError("Isolate metadata merge changed the AST row count")
    unmatched = merged["__isolate_merge"].ne("both")
    if unmatched.any():
        raise NormalizationError(
            f"Isolate metadata has no exact {join_key} match for {int(unmatched.sum())} AST rows"
        )

    for field in isolate_fields:
        if field == join_key:
            continue
        right = merged[f"__isolate_{field}"]
        if field in {
            "isolate_id",
            "biosample_accession",
            "assembly_accession",
            "scientific_name",
            "bioproject_accession",
            "collection_date",
            "geo_loc_name",
            "country",
        }:
            conflict = (
                ~blank_mask(out[field])
                & ~blank_mask(right)
                & out[field].map(normalized_metadata_value).ne(
                    right.map(normalized_metadata_value)
                )
            )
            if conflict.any():
                raise NormalizationError(
                    f"AST and isolate metadata disagree on {field} for {int(conflict.sum())} rows"
                )
        out[field] = coalesce_columns(out[field], right)
    return out


def _append_reason(reasons: pd.Series, mask: pd.Series, code: str) -> pd.Series:
    addition = pd.Series(code, index=reasons.index, dtype="string")
    return reasons.mask(mask & reasons.eq(""), addition).mask(
        mask & reasons.ne(""), reasons + ";" + addition
    )


def load_country_vocabulary(path: Path) -> set[str]:
    """Load the frozen INSDC /geo_loc_name vocabulary (current and historical)."""
    frame = pd.read_csv(path, sep="\t", dtype=str)
    if list(frame.columns) != ["name", "status"]:
        raise NormalizationError(f"Country vocabulary has unexpected columns: {list(frame.columns)}")
    if not frame["status"].isin(["current", "historical"]).all():
        raise NormalizationError("Country vocabulary status must be current or historical")
    names = set(frame["name"].str.strip())
    if len(names) != len(frame):
        raise NormalizationError("Country vocabulary contains duplicate names")
    return names


def derive_country(out: pd.DataFrame, vocabulary: set[str] | None) -> pd.DataFrame:
    """Amendment 006: fill country from the INSDC prefix of geo_loc_name.

    INSDC defines geo_loc_name as "<name>[:<region>][, <locality>]" where <name>
    is drawn from a controlled vocabulary, so the prefix is a structured field,
    not free text. An explicit country always wins. A prefix outside the frozen
    vocabulary -- including missingness tokens such as "not collected" -- leaves
    country absent rather than guessing.
    """
    explicit = ~blank_mask(out["country"])
    source = pd.Series("absent", index=out.index, dtype="string")
    source = source.mask(explicit, "explicit")
    if vocabulary is not None:
        prefix = out["geo_loc_name"].astype("string").str.split(":").str[0].str.strip()
        derived = ~explicit & prefix.isin(vocabulary).fillna(False).astype(bool)
        out.loc[derived, "country"] = prefix[derived]
        source = source.mask(derived, "insdc_geo_loc_name_prefix")
    out["country_source"] = source
    return out


def normalize_ast(
    ast: pd.DataFrame,
    isolates: pd.DataFrame | None = None,
    source_fingerprint: str = "UNVERSIONED_TEST_INPUT",
    source_release: str = "UNSPECIFIED",
    country_vocabulary: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    ast = canonicalize_frame(ast, "AST input")
    isolates = canonicalize_frame(isolates, "isolate input") if isolates is not None else None

    out = pd.DataFrame(index=ast.index)
    for field in ALIASES:
        if field not in {"mic", "mic_secondary", "disk_diffusion", "disk_diffusion_secondary"}:
            out[field] = first_column(ast, field)

    _assert_consistent_identity_mapping(out, "AST input")
    if isolates is not None:
        out = _merge_isolate_metadata(out, isolates)

    mic_raw = first_column(ast, "mic")
    disk_raw = first_column(ast, "disk_diffusion")
    mic_secondary = first_column(ast, "mic_secondary")
    disk_secondary = first_column(ast, "disk_diffusion_secondary")
    mic = pd.to_numeric(mic_raw, errors="coerce")
    disk = pd.to_numeric(disk_raw, errors="coerce")
    has_mic_raw = ~blank_mask(mic_raw)
    has_disk_raw = ~blank_mask(disk_raw)
    has_mic_secondary = ~blank_mask(mic_secondary)
    has_disk_secondary = ~blank_mask(disk_secondary)
    has_mic = mic.notna()
    has_disk = disk.notna()

    out["ast_measurement_type"] = pd.Series(pd.NA, index=out.index, dtype="string")
    out.loc[has_mic & ~has_disk, "ast_measurement_type"] = "MIC"
    out.loc[has_disk & ~has_mic, "ast_measurement_type"] = "ZONE"
    out["ast_value"] = mic.where(has_mic, disk)
    out["raw_ast_value"] = mic_raw.where(has_mic_raw, disk_raw)
    out["mic_secondary"] = mic_secondary
    out["disk_diffusion_secondary"] = disk_secondary
    out["ast_unit"] = pd.Series(pd.NA, index=out.index, dtype="string")
    out.loc[has_mic & ~has_disk, "ast_unit"] = "mg/L"
    out.loc[has_disk & ~has_mic, "ast_unit"] = "mm"
    out["raw_antibiotic"] = out["antibiotic"]
    out["raw_submitted_ast_category"] = out["submitted_ast_category"]
    out["raw_measurement_sign"] = out["measurement_sign"]
    out["measurement_sign"] = out["measurement_sign"].replace({"==": "="})
    out["submitted_ast_category"] = out["submitted_ast_category"].map(normalize_category)
    out["antibiotic"] = clean_string_series(out["antibiotic"]).str.lower()
    out["source_fingerprint"] = source_fingerprint
    out["source_row_number"] = pd.Series(range(2, len(out) + 2), index=out.index, dtype="Int64")

    fallback_ids = ast.apply(
        lambda row: canonical_raw_record_id(row, ast.columns), axis=1
    ).astype("string")
    native_id = clean_string_series(out["source_native_record_id"])
    native_present = ~blank_mask(native_id)
    out["source_ast_record_id"] = native_id.where(native_present, fallback_ids)
    out["source_ast_record_id_origin"] = pd.Series(
        np.where(native_present, "native_id", "canonical_raw_row_sha256"),
        index=out.index,
        dtype="string",
    )

    reason = pd.Series(pd.NA, index=out.index, dtype="string")
    scientific_name = out["scientific_name"].fillna("").str.lower()
    wrong_species = ~scientific_name.str.match(r"^escherichia coli(?:\s|$)")
    reason = reason.mask(wrong_species, "not_escherichia_coli")
    reason = reason.mask(reason.isna() & has_mic_raw & has_disk_raw, "multiple_measurement_types")
    has_secondary = has_mic_secondary | has_disk_secondary
    reason = reason.mask(
        reason.isna() & has_secondary,
        "secondary_measurement_requires_reconciliation",
    )
    nonnumeric = (has_mic_raw & ~has_mic) | (has_disk_raw & ~has_disk)
    reason = reason.mask(
        reason.isna() & nonnumeric,
        "unsupported_nonnumeric_or_combination_measurement",
    )
    reason = reason.mask(reason.isna() & ~has_mic_raw & ~has_disk_raw, "missing_quantitative_ast")
    numeric_values = pd.to_numeric(out["ast_value"], errors="coerce")
    finite_positive = numeric_values.gt(0) & np.isfinite(numeric_values)
    reason = reason.mask(
        reason.isna() & numeric_values.notna() & ~finite_positive,
        "nonpositive_or_nonfinite_measurement",
    )
    missing_sign = blank_mask(out["measurement_sign"])
    reason = reason.mask(reason.isna() & missing_sign, "missing_measurement_sign")
    invalid_sign = ~out["measurement_sign"].isin({"=", "<", "<=", ">", ">="})
    reason = reason.mask(reason.isna() & ~missing_sign & invalid_sign, "invalid_measurement_sign")
    reason = reason.mask(reason.isna() & blank_mask(out["antibiotic"]), "missing_antibiotic")
    reason = reason.mask(reason.isna() & blank_mask(out["isolate_id"]), "missing_isolate_id")
    reason = reason.mask(
        reason.isna() & blank_mask(out["biosample_accession"]),
        "missing_biosample_accession",
    )
    duplicate_record_id = out["source_ast_record_id"].duplicated(keep=False)
    reason = reason.mask(
        reason.isna() & duplicate_record_id,
        "duplicate_source_ast_record_id",
    )

    breakpoint_reasons = pd.Series("", index=out.index, dtype="string")
    breakpoint_reasons = _append_reason(
        breakpoint_reasons,
        blank_mask(out["ast_method"]),
        "missing_ast_method",
    )
    method_present = ~blank_mask(out["ast_method"])
    valid_methods = pd.Series(
        [
            valid_ast_method(method, measurement_type)
            for method, measurement_type in zip(
                out["ast_method"], out["ast_measurement_type"], strict=True
            )
        ],
        index=out.index,
        dtype="boolean",
    )
    breakpoint_reasons = _append_reason(
        breakpoint_reasons,
        method_present & ~valid_methods,
        "unsupported_ast_method_or_measurement_mismatch",
    )
    zone_rows = out["ast_measurement_type"].eq("ZONE")
    breakpoint_reasons = _append_reason(
        breakpoint_reasons,
        zone_rows & blank_mask(out["disk_content"]),
        "missing_disk_content",
    )
    disk_present = ~blank_mask(out["disk_content"])
    valid_disks = out["disk_content"].map(valid_disk_content).astype("boolean")
    breakpoint_reasons = _append_reason(
        breakpoint_reasons,
        zone_rows & disk_present & ~valid_disks,
        "invalid_disk_content",
    )
    breakpoint_reasons = _append_reason(
        breakpoint_reasons,
        blank_mask(out["ast_testing_date"]),
        "missing_ast_testing_date",
    )
    testing_date_present = ~blank_mask(out["ast_testing_date"])
    valid_testing_dates = out["ast_testing_date"].map(valid_iso_date_or_datetime).astype("boolean")
    breakpoint_reasons = _append_reason(
        breakpoint_reasons,
        testing_date_present & ~valid_testing_dates,
        "invalid_ast_testing_date",
    )
    out["breakpoint_eligible"] = breakpoint_reasons.eq("")
    out["breakpoint_ineligibility_reasons"] = breakpoint_reasons.mask(
        breakpoint_reasons.eq(""), pd.NA
    )

    # Eligibility for the interval-censored MIC endpoint (docs/ENDPOINT_AMENDMENT.md).
    # This deliberately does not consult AST method, testing date, or disk
    # potency: no measurement is interpreted under that endpoint, so those
    # fields are not part of its reference standard. They remain recorded as
    # absent above and still gate the categorical secondary endpoint.
    mic_reasons = pd.Series("", index=out.index, dtype="string")
    mic_reasons = _append_reason(
        mic_reasons,
        ~out["ast_measurement_type"].eq("MIC"),
        "not_a_mic_measurement",
    )
    mic_reasons = _append_reason(
        mic_reasons,
        blank_mask(out["ast_value"].astype("string")),
        "missing_mic_value",
    )
    mic_reasons = _append_reason(
        mic_reasons,
        blank_mask(out["measurement_sign"].astype("string")),
        "missing_measurement_sign",
    )
    mic_reasons = _append_reason(
        mic_reasons,
        ~out["ast_unit"].astype("string").str.strip().eq("mg/L"),
        "unexpected_mic_unit",
    )
    out["mic_eligible"] = mic_reasons.eq("")
    out["mic_ineligibility_reasons"] = mic_reasons.mask(mic_reasons.eq(""), pd.NA)

    out = derive_country(out, country_vocabulary)
    out.insert(0, "source_dataset", "NCBI Pathogen Detection AST")
    out.insert(1, "source_release", source_release)
    keep_columns = [
        "source_dataset",
        "source_release",
        "source_fingerprint",
        "source_row_number",
        "source_ast_record_id",
        "source_ast_record_id_origin",
        "source_record_checksum",
        "isolate_id",
        "biosample_accession",
        "assembly_accession",
        "taxgroup_name",
        "scientific_name",
        "raw_antibiotic",
        "antibiotic",
        "raw_submitted_ast_category",
        "submitted_ast_category",
        "ast_measurement_type",
        "raw_measurement_sign",
        "measurement_sign",
        "raw_ast_value",
        "ast_value",
        "ast_unit",
        "mic_secondary",
        "disk_diffusion_secondary",
        "ast_method",
        "submitted_standard",
        "ast_platform",
        "ast_vendor",
        "ast_reagent",
        "disk_content",
        "ast_testing_date",
        "breakpoint_eligible",
        "breakpoint_ineligibility_reasons",
        "mic_eligible",
        "mic_ineligibility_reasons",
        "collection_date",
        "ast_record_creation_date",
        "target_creation_date",
        "geo_loc_name",
        "country",
        "country_source",
        "bioproject_accession",
        "epi_type",
        "isolation_type",
        "isolation_source",
        "host",
    ]
    included = out.loc[reason.isna(), keep_columns].reset_index(drop=True)
    excluded = out.loc[reason.notna(), keep_columns].copy()
    excluded.insert(0, "exclusion_reason", reason[reason.notna()].values)
    return included, excluded.reset_index(drop=True)


def _verify_hash(path: Path, expected: str, label: str) -> str:
    if not path.is_file():
        raise NormalizationError(f"Acquisition manifest references missing {label}: {path}")
    actual = sha256_file(path)
    if actual != expected:
        raise NormalizationError(
            f"Acquisition manifest hash mismatch for {label}: expected {expected}, observed {actual}"
        )
    return actual


def verify_acquisition_provenance(
    manifest_path: Path,
    ast_path: Path,
    isolates_path: Path | None,
    ast_rows: int,
    isolate_rows: int | None,
    amendment_path: Path | None = None,
) -> dict[str, Any]:
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise NormalizationError(f"Cannot read acquisition manifest {manifest_path}: {exc}") from exc
    if manifest.get("status") != "complete":
        raise NormalizationError("Acquisition manifest status is not complete")

    retrieval = manifest.get("retrieval") or {}
    outputs = retrieval.get("outputs") or {}
    verified: dict[str, str] = {}
    for logical_name, path, rows in (
        ("ast", ast_path, ast_rows),
        ("isolates", isolates_path, isolate_rows),
    ):
        if path is None:
            continue
        record = outputs.get(logical_name) or {}
        expected_hash = record.get("sha256")
        expected_rows = record.get("rows")
        if not expected_hash or expected_rows is None:
            raise NormalizationError(
                f"Acquisition manifest lacks hash or row count for {logical_name} output"
            )
        verified[f"output_{logical_name}"] = _verify_hash(
            path, str(expected_hash), f"{logical_name} output"
        )
        if int(expected_rows) != int(rows):
            raise NormalizationError(
                f"Acquisition manifest row-count mismatch for {logical_name}: "
                f"expected {expected_rows}, observed {rows}"
            )

    base_dir = manifest_path.parent
    for logical_name in ("ast", "isolates"):
        record = (manifest.get("queries") or {}).get(logical_name) or {}
        if not record.get("file") or not record.get("sha256"):
            raise NormalizationError(
                f"Acquisition manifest lacks file or hash for {logical_name} query"
            )
        artifact = base_dir / str(record["file"])
        verified[f"query_{logical_name}"] = _verify_hash(
            artifact, str(record["sha256"]), f"{logical_name} query"
        )
    metadata_record = retrieval.get("source_table_metadata") or {}
    if not metadata_record.get("file") or not metadata_record.get("sha256"):
        raise NormalizationError(
            "Acquisition manifest lacks source-table-metadata file or hash"
        )
    artifact = base_dir / str(metadata_record["file"])
    verified["source_table_metadata"] = _verify_hash(
        artifact,
        str(metadata_record["sha256"]),
        "source table metadata",
    )

    manifest_hash = sha256_file(manifest_path)
    result: dict[str, Any] = {
        "path": str(manifest_path),
        "sha256": manifest_hash,
        "verified_artifacts": verified,
        "amendment": None,
    }
    amendment: dict[str, Any] | None = None
    if amendment_path is not None:
        try:
            amendment = json.loads(amendment_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise NormalizationError(
                f"Cannot read acquisition manifest amendment {amendment_path}: {exc}"
            ) from exc
        if amendment.get("amendment_status") != "complete":
            raise NormalizationError("Acquisition amendment status is not complete")
        original_hash = ((amendment.get("original_manifest") or {}).get("sha256"))
        if original_hash != manifest_hash:
            raise NormalizationError(
                "Acquisition amendment does not reference the supplied acquisition manifest hash"
            )
        unchanged = amendment.get("unchanged_artifacts") or {}
        for amendment_key, verified_key in (
            ("ast_sha256", "output_ast"),
            ("isolates_sha256", "output_isolates"),
            ("source_table_metadata_sha256", "source_table_metadata"),
        ):
            if amendment_key not in unchanged:
                raise NormalizationError(
                    f"Acquisition amendment lacks required unchanged hash {amendment_key}"
                )
            if unchanged[amendment_key] != verified.get(verified_key):
                raise NormalizationError(
                    f"Acquisition amendment {amendment_key} disagrees with verified provenance"
                )
        result["amendment"] = {
            "path": str(amendment_path),
            "sha256": sha256_file(amendment_path),
        }

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
            required_corrections["/queries/parameters/antibiotics_csv"] = expected_parameter
    executable = str(((retrieval.get("client") or {}).get("executable")) or "")
    portable_executable = executable.replace("\\", "/").rsplit("/", 1)[-1]
    if executable and executable != portable_executable:
        required_corrections["/retrieval/client/executable"] = portable_executable
    if required_corrections:
        if amendment is None:
            raise NormalizationError(
                "Acquisition manifest has provenance serialization errors but no verified amendment"
            )
        corrections = {
            str(item.get("json_pointer")): str(item.get("corrected_value"))
            for item in amendment.get("corrections", [])
            if isinstance(item, dict)
        }
        for pointer, expected in required_corrections.items():
            if corrections.get(pointer) != expected:
                raise NormalizationError(
                    f"Acquisition amendment does not supply the required correction for {pointer}"
                )
    return result


def _temporary_path(destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=destination.parent,
        delete=False,
    )
    handle.close()
    return Path(handle.name)


def write_csv_atomic(frame: pd.DataFrame, destination: Path) -> None:
    temporary = _temporary_path(destination)
    try:
        frame.to_csv(temporary, index=False, encoding="utf-8", lineterminator="\n")
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def write_json_atomic(payload: dict[str, Any], destination: Path) -> None:
    temporary = _temporary_path(destination)
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _missing_counts(frame: pd.DataFrame, fields: Iterable[str]) -> dict[str, int]:
    return {field: int(blank_mask(frame[field]).sum()) for field in fields}


def _reason_counts(series: pd.Series) -> dict[str, int]:
    counts: dict[str, int] = {}
    for value in series.dropna().astype(str):
        for reason in value.split(";"):
            counts[reason] = counts.get(reason, 0) + 1
    return dict(sorted(counts.items()))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ast", type=Path, required=True, help="NCBI AST Browser or BigQuery CSV/TSV export")
    parser.add_argument("--isolates", type=Path, help="Optional NCBI Isolates Browser export for assembly metadata")
    parser.add_argument("--output", type=Path, default=Path("data/interim/ncbi_ast_normalized.csv"))
    parser.add_argument("--exclusions", type=Path, default=Path("data/interim/ncbi_ast_exclusions.csv"))
    parser.add_argument("--manifest", type=Path, default=Path("data/interim/ncbi_ast_manifest.json"))
    parser.add_argument(
        "--country-vocabulary",
        type=Path,
        default=Path("config/insdc_geo_loc_name_vocabulary.tsv"),
        help="Frozen INSDC /geo_loc_name vocabulary used to derive country (amendment 006).",
    )
    parser.add_argument(
        "--acquisition-manifest",
        type=Path,
        help="Completed acquisition manifest to authenticate source inputs (auto-detected beside --ast).",
    )
    parser.add_argument(
        "--acquisition-amendment",
        type=Path,
        help="Optional append-only correction that references the acquisition manifest SHA-256.",
    )
    parser.add_argument(
        "--source-release",
        help="Source snapshot/release identifier; defaults to the source file SHA-256",
    )
    parser.add_argument(
        "--allow-unverified-input",
        action="store_true",
        help=(
            "Permit inputs without an acquisition manifest only for explicitly synthetic "
            "software tests; such output is marked SMOKE_TEST_ONLY_UNVERIFIED_INPUT."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ast = read_table(args.ast)
    isolates = read_table(args.isolates) if args.isolates else None
    source_hash = sha256_file(args.ast)
    source_release = args.source_release or source_hash

    acquisition_manifest = args.acquisition_manifest
    inferred_manifest = args.ast.parent / DEFAULT_ACQUISITION_MANIFEST
    if acquisition_manifest is None and inferred_manifest.is_file():
        acquisition_manifest = inferred_manifest
    acquisition_amendment = args.acquisition_amendment
    inferred_amendment = args.ast.parent / DEFAULT_ACQUISITION_AMENDMENT
    if acquisition_amendment is None and inferred_amendment.is_file():
        acquisition_amendment = inferred_amendment
    if acquisition_amendment is not None and acquisition_manifest is None:
        raise NormalizationError("An acquisition amendment requires an acquisition manifest")
    if acquisition_manifest is None and not args.allow_unverified_input:
        raise NormalizationError(
            "No completed acquisition manifest was supplied or found beside the AST input. "
            "Use --allow-unverified-input only for synthetic software tests."
        )
    provenance = (
        verify_acquisition_provenance(
            acquisition_manifest,
            args.ast,
            args.isolates,
            len(ast),
            len(isolates) if isolates is not None else None,
            acquisition_amendment,
        )
        if acquisition_manifest is not None
        else None
    )

    country_vocabulary = load_country_vocabulary(args.country_vocabulary)
    included, excluded = normalize_ast(
        ast,
        isolates,
        source_fingerprint=source_hash,
        source_release=source_release,
        country_vocabulary=country_vocabulary,
    )
    write_csv_atomic(included, args.output)
    write_csv_atomic(excluded, args.exclusions)

    duplicate_groups = (
        included.loc[
            included.duplicated(["isolate_id", "antibiotic"], keep=False),
            ["isolate_id", "antibiotic"],
        ]
        .drop_duplicates()
        .shape[0]
    )
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "mapping_version": NORMALIZATION_SCHEMA_VERSION,
        "normalization_schema_version": NORMALIZATION_SCHEMA_VERSION,
        "normalizer": {
            "script": "scripts/normalize_ncbi_ast.py",
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "python_version": sys.version.split()[0],
            "pandas_version": pd.__version__,
            "numpy_version": np.__version__,
            "arguments": {
                "ast": str(args.ast),
                "isolates": str(args.isolates) if args.isolates else None,
                "source_release": source_release,
                "acquisition_manifest": (
                    str(acquisition_manifest) if acquisition_manifest else None
                ),
                "acquisition_amendment": (
                    str(acquisition_amendment) if acquisition_amendment else None
                ),
                "allow_unverified_input": bool(args.allow_unverified_input),
            },
        },
        "source": {
            "path": str(args.ast),
            "sha256": source_hash,
            "bytes": args.ast.stat().st_size,
            "rows": int(len(ast)),
            "columns": list(ast.columns),
            "release": source_release,
        },
        "isolates_source": (
            {
                "path": str(args.isolates),
                "sha256": sha256_file(args.isolates),
                "bytes": args.isolates.stat().st_size,
                "rows": int(len(isolates)),
                "columns": list(isolates.columns),
            }
            if args.isolates and isolates is not None
            else None
        ),
        "acquisition_provenance": provenance,
        "input_schemas": {
            "ast": list(ast.columns),
            "isolates": list(isolates.columns) if isolates is not None else None,
        },
        "output_schema": list(included.columns),
        "outputs": {
            "normalized": {
                "path": str(args.output),
                "rows": int(len(included)),
                "bytes": args.output.stat().st_size,
                "sha256": sha256_file(args.output),
            },
            "exclusions": {
                "path": str(args.exclusions),
                "rows": int(len(excluded)),
                "bytes": args.exclusions.stat().st_size,
                "sha256": sha256_file(args.exclusions),
            },
        },
        "n_source_rows": int(len(ast)),
        "n_normalized_rows": int(len(included)),
        "n_excluded_rows": int(len(excluded)),
        "exclusion_counts": excluded["exclusion_reason"].value_counts().sort_index().to_dict(),
        "country_vocabulary": {
            "path": str(args.country_vocabulary),
            "sha256": sha256_file(args.country_vocabulary),
            "amendment": "docs/COUNTRY_AND_POPULATION_AMENDMENT.md",
        },
        "country_source_counts": {
            str(key): int(value)
            for key, value in included["country_source"].value_counts().sort_index().items()
        },
        "n_breakpoint_eligible_rows": int(included["breakpoint_eligible"].sum()),
        "breakpoint_ineligibility_counts": _reason_counts(
            included["breakpoint_ineligibility_reasons"]
        ),
        "n_mic_eligible_rows": int(included["mic_eligible"].sum()),
        "mic_ineligibility_counts": _reason_counts(included["mic_ineligibility_reasons"]),
        "mic_eligibility_scope": (
            "Quantitative-measurement readiness for the interval-censored MIC endpoint only: an "
            "MIC measurement type with a value, an explicit comparator, and mg/L units. "
            "Eligibility here does not establish assay method, AST testing date, clinical "
            "indication, laboratory site, cohort completeness, or analysis readiness, and it "
            "confers no clinical category. See docs/ENDPOINT_AMENDMENT.md."
        ),
        "breakpoint_eligibility_scope": (
            "Source-stage assay metadata only: an explicitly supported method consistent with "
            "measurement type and a valid ISO testing date for every record, plus positive disk "
            "potency with an explicit mass unit for zone records. Eligibility here does not establish a clinical "
            "indication, breakpoint authority/version, cohort completeness, or analysis readiness."
        ),
        "field_missing_counts": _missing_counts(
            included,
            (
                "source_record_checksum",
                "assembly_accession",
                "ast_method",
                "disk_content",
                "ast_testing_date",
                "collection_date",
                "geo_loc_name",
                "country",
                "bioproject_accession",
            ),
        ),
        "n_repeated_isolate_antibiotic_groups": int(duplicate_groups),
        "source_record_id_origin_counts": (
            included["source_ast_record_id_origin"].value_counts().sort_index().to_dict()
        ),
        "scientific_status": (
            "SOURCE_NORMALIZATION_ONLY"
            if provenance is not None
            else "SMOKE_TEST_ONLY_UNVERIFIED_INPUT"
        ),
        "limitations": [
            "NCBI AST values and metadata are submitter-provided and not method-verified by NCBI.",
            "No breakpoint version is inferred from a testing-standard name.",
            "The submitted phenotype is preserved separately from any later independent interpretation.",
            "Country is derived from the INSDC geo_loc_name prefix only when the prefix is in the frozen INSDC vocabulary; geo_loc_name itself is preserved verbatim (amendment 006).",
            "Repeated isolate-antibiotic records require a prespecified reconciliation rule.",
            "The supported AST-method vocabulary is deliberately conservative; new source terms require an audited mapping-version update.",
        ],
    }
    write_json_atomic(manifest, args.manifest)
    print(
        f"Normalized {len(included)} AST rows; structurally excluded {len(excluded)}; "
        f"source-stage breakpoint eligible {manifest['n_breakpoint_eligible_rows']}; "
        f"MIC-endpoint eligible {manifest['n_mic_eligible_rows']}. "
        f"Manifest: {args.manifest}"
    )


if __name__ == "__main__":
    main()
