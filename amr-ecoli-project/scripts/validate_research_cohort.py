"""Validate a curated isolate-level AMR modeling cohort.

The input is the one-row-per-isolate-per-antibiotic modeling table described in
docs/STUDY_PROTOCOL.md. This validator intentionally rejects the aggregate NCBI
AST_phenotypes table because it lacks quantitative AST provenance.

Run from the repository root:
    python scripts/validate_research_cohort.py --cohort data/curated/cohort.csv
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import date
from pathlib import Path

import pandas as pd

try:
    from scripts.apply_breakpoints import (
        ADDED_OUTPUT_COLUMNS,
        BreakpointError,
        interpret_rows,
        load_breakpoint_rules,
    )
except ModuleNotFoundError:  # Direct execution: python scripts/validate_research_cohort.py
    from apply_breakpoints import (  # type: ignore[no-redef]
        ADDED_OUTPUT_COLUMNS,
        BreakpointError,
        interpret_rows,
        load_breakpoint_rules,
    )


# Columns every endpoint profile requires. These carry isolate identity, genome
# provenance, the quantitative measurement itself, and the evaluation structure
# that keeps near neighbours and duplicates from crossing a split boundary.
BASE_REQUIRED_COLUMNS = [
    "isolate_id",
    "biosample_accession",
    "assembly_accession",
    "scientific_name",
    "species_method",
    "species_result",
    "genome_qc_status",
    "antibiotic",
    "source_dataset",
    "source_release",
    "source_fingerprint",
    "source_row_number",
    "source_ast_record_id",
    "ast_measurement_type",
    "measurement_sign",
    "raw_ast_value",
    "ast_value",
    "ast_unit",
    "collection_date",
    "country",
    "bioproject_accession",
    "specimen_source",
    "deduplication_group",
    "lineage_group",
    "genomic_cluster",
    "evaluation_split",
    "intended_use_population",
    "genome_source",
    "genome_source_accession",
    "genome_sha256",
]

# Columns required only to interpret a measurement into a clinical category.
# docs/ENDPOINT_AMENDMENT.md records why no public source supplies these and why
# the categorical endpoint is therefore a blocked secondary rather than primary.
CATEGORICAL_ONLY_COLUMNS = [
    "ast_method",
    "ast_testing_date",
    "submitted_ast_category",
    "ast_category",
    "breakpoint_standard",
    "breakpoint_version",
    "breakpoint_rule_id",
    "breakpoint_rule_sha256",
    "breakpoint_artifact_sha256",
    "site",
    "clinical_indication",
    "surveillance_network",
]

MIC_REGRESSION_COLUMNS = list(BASE_REQUIRED_COLUMNS)
CATEGORICAL_SIR_COLUMNS = BASE_REQUIRED_COLUMNS + CATEGORICAL_ONLY_COLUMNS

ENDPOINT_PROFILES = {
    "mic_regression": MIC_REGRESSION_COLUMNS,
    "categorical_sir": CATEGORICAL_SIR_COLUMNS,
}

# Amendment 006: country, collection date and specimen source stay absent when
# the source does not supply them (never proxied). Under the MIC endpoint such
# rows remain usable for lineage-grouped evaluation; the analyses that need the
# field (geographic and forward-time holdouts, specimen-type subgroups) exclude
# them. Their counts are reported as warnings, not failures.
ABSENT_ALLOWED_FOR_MIC = frozenset({"collection_date", "country", "specimen_source"})

# INSDC collection dates carry the submitter's precision: YYYY, YYYY-MM or
# YYYY-MM-DD. Missing month or day is never filled in.
PARTIAL_DATE_RE = re.compile(r"^(\d{4})(?:-(\d{2})(?:-(\d{2}))?)?$")
MIN_COLLECTION_YEAR = 1900


def valid_partial_date(value: object, today: date | None = None) -> bool:
    match = PARTIAL_DATE_RE.match(str(value).strip())
    if not match:
        return False
    year, month, day = (int(part) if part else None for part in match.groups())
    today = today or date.today()
    if year < MIN_COLLECTION_YEAR or year > today.year:
        return False
    if month is not None and not 1 <= month <= 12:
        return False
    if day is not None:
        try:
            if date(year, month, day) > today:
                return False
        except ValueError:
            return False
    elif month is not None and (year, month) > (today.year, today.month):
        return False
    return True


# Fail closed: an unspecified profile keeps the strictest historical contract.
DEFAULT_ENDPOINT_PROFILE = "categorical_sir"

# Retained so existing callers and tests that import the strict contract keep
# the behaviour they had before docs/ENDPOINT_AMENDMENT.md.
REQUIRED_COLUMNS = CATEGORICAL_SIR_COLUMNS

ALLOWED_CATEGORIES = {"S", "I", "R"}
ALLOWED_MEASUREMENTS = {"MIC", "ZONE"}
ALLOWED_STANDARDS = {"CLSI", "EUCAST"}
ALLOWED_SPLITS = {"development", "external"}
ALLOWED_POPULATIONS = {"human_clinical", "non_human_or_environmental", "undetermined"}
ALLOWED_SIGNS = {"=", "==", "<", "<=", ">", ">="}
ASSEMBLY_RE = re.compile(r"^GC[AF]_\d+\.\d+$")
RUN_RE = re.compile(r"^[SED]RR\d+$")
GENOME_SOURCES = {"registered_assembly", "raw_reads"}
# Genomes assembled here from raw reads have no NCBI assembly accession; their
# provenance is the run accession plus the checksum of the assembly analysed.
BLANK_ALLOWED_FOR_RAW_READS = {"assembly_accession"}
BIOSAMPLE_RE = re.compile(r"^SAM[NED][A-Z]?\d+$")


def resolve_endpoint_profile(study_config: dict | None, override: str | None = None) -> str:
    """Pick the endpoint profile, preferring an explicit override.

    An absent or unrecognised value falls back to the strict categorical
    contract so a malformed config cannot silently relax the gate.
    """
    if override:
        candidate = str(override).strip()
    elif study_config:
        candidate = str(study_config.get("primary_endpoint", "")).strip()
    else:
        candidate = ""
    if candidate in ENDPOINT_PROFILES:
        return candidate
    return DEFAULT_ENDPOINT_PROFILE


def configured_antibiotic_endpoints(study_config: dict) -> dict:
    """Read categorical endpoints from the secondary block, then the legacy key."""
    secondary = study_config.get("secondary_endpoint", {})
    if isinstance(secondary, dict) and secondary.get("antibiotic_endpoints"):
        return secondary["antibiotic_endpoints"]
    return study_config.get("antibiotic_endpoints", {})


def blank_mask(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.strip().str.lower()
    return series.isna() | text.isin({"", "nan", "none", "null"})


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_breakpoint_artifact(
    df: pd.DataFrame,
    study_config: dict,
    config_root: Path,
    profile: str = DEFAULT_ENDPOINT_PROFILE,
) -> list[str]:
    failures: list[str] = []
    if profile == "mic_regression":
        # No measurement is interpreted under this endpoint, so a breakpoint
        # artifact is not part of the reference standard. The artifact stays
        # frozen in config for the blocked categorical secondary.
        return failures
    ast_config = study_config.get("ast", {})
    standard = str(ast_config.get("interpretation_standard", "")).strip()
    version = str(ast_config.get("breakpoint_version", "")).strip()
    configured_hash = str(ast_config.get("breakpoint_table_sha256", "")).strip().lower()
    artifact_value = str(ast_config.get("breakpoint_table", "")).strip()
    if not standard or not version or not artifact_value:
        return ["Study config must define AST standard, breakpoint version, and breakpoint table."]
    if configured_hash.startswith("required") or not re.fullmatch(r"[0-9a-f]{64}", configured_hash):
        failures.append("Study config does not contain a frozen breakpoint-table SHA-256.")
    artifact = Path(artifact_value)
    if not artifact.is_absolute():
        artifact = config_root / artifact
    if not artifact.exists():
        failures.append(f"Configured breakpoint table does not exist: {artifact}")
    elif re.fullmatch(r"[0-9a-f]{64}", configured_hash):
        observed_hash = sha256_file(artifact)
        if observed_hash != configured_hash:
            failures.append(
                f"Breakpoint table SHA-256 mismatch: configured={configured_hash}, observed={observed_hash}"
            )
    if "breakpoint_standard" in df.columns:
        observed_standards = set(df["breakpoint_standard"].astype(str).str.upper())
        if observed_standards != {standard.upper()}:
            failures.append(
                f"Cohort breakpoint standard {sorted(observed_standards)} does not match config {standard}."
            )
    if "breakpoint_version" in df.columns:
        observed_versions = set(df["breakpoint_version"].astype(str))
        if observed_versions != {version}:
            failures.append(
                f"Cohort breakpoint version {sorted(observed_versions)} does not match config {version}."
            )
    if "breakpoint_artifact_sha256" in df.columns and re.fullmatch(r"[0-9a-f]{64}", configured_hash):
        row_hashes = set(df["breakpoint_artifact_sha256"].astype(str).str.lower())
        if row_hashes != {configured_hash}:
            failures.append("Cohort rows do not reference the configured breakpoint artifact SHA-256.")
    configured_drugs = set(map(str, study_config.get("primary_antibiotics", [])))
    observed_drugs = set(df["antibiotic"].astype(str)) if "antibiotic" in df.columns else set()
    missing_drugs = sorted(configured_drugs.difference(observed_drugs))
    if missing_drugs:
        failures.append(f"Configured primary antibiotics absent from cohort: {missing_drugs}")
    for drug, endpoint in configured_antibiotic_endpoints(study_config).items():
        expected_indication = str(endpoint.get("clinical_indication", "")).strip()
        if not expected_indication or drug not in observed_drugs:
            continue
        observed_indications = set(
            df.loc[df["antibiotic"].astype(str).eq(str(drug)), "clinical_indication"]
            .astype(str)
            .str.strip()
        )
        if observed_indications != {expected_indication}:
            failures.append(
                f"Cohort clinical indication for {drug} {sorted(observed_indications)} "
                f"does not match configured endpoint {expected_indication}."
            )
    return failures


def validate_breakpoint_derivation(
    df: pd.DataFrame, rules_path: Path, artifact_path: Path
) -> list[str]:
    """Re-derive every category and compare rule-level provenance exactly."""
    failures: list[str] = []
    try:
        rules = load_breakpoint_rules(rules_path, artifact_path)
    except BreakpointError as exc:
        return [f"Breakpoint rules failed validation: {exc}"]

    source_rows = []
    for row in df.to_dict(orient="records"):
        source_row = {
            key: value for key, value in row.items() if key not in ADDED_OUTPUT_COLUMNS
        }
        # Curated cohorts retain the assay evidence but not the normalizer's
        # derived eligibility columns. Reconstruct the candidate flag here;
        # interpret_rows independently revalidates method/date/disk metadata and
        # rejects any inconsistent flag before applying a rule.
        source_row.setdefault("disk_content", "")
        source_row["breakpoint_eligible"] = True
        source_row["breakpoint_ineligibility_reasons"] = ""
        source_rows.append(source_row)
    interpreted, excluded = interpret_rows(source_rows, rules)
    if excluded:
        details = [
            f"{row.get('source_ast_record_id', '<missing>')}:{row['exclusion_reason']}"
            for row in excluded[:20]
        ]
        failures.append(
            "Cohort rows cannot be independently re-interpreted by the frozen rules: "
            + str(details)
        )

    expected_by_id = {str(row["source_ast_record_id"]): row for row in interpreted}
    mismatches: list[str] = []
    comparisons = {
        "ast_category": "ast_category",
        "breakpoint_rule_id": "breakpoint_rule_id",
        "breakpoint_rule_sha256": "breakpoint_rule_row_sha256",
        "breakpoint_artifact_sha256": "breakpoint_artifact_sha256",
        "breakpoint_standard": "breakpoint_standard",
        "breakpoint_version": "breakpoint_version",
    }
    for observed in df.to_dict(orient="records"):
        record_id = str(observed.get("source_ast_record_id", ""))
        expected = expected_by_id.get(record_id)
        if expected is None:
            continue
        for observed_field, expected_field in comparisons.items():
            if str(observed.get(observed_field, "")).strip() != str(
                expected.get(expected_field, "")
            ).strip():
                mismatches.append(f"{record_id}:{observed_field}")
    if mismatches:
        failures.append(
            "Cohort breakpoint-derived values disagree with independent re-interpretation: "
            + str(sorted(mismatches)[:20])
        )
    return failures


def validate_cohort(
    df: pd.DataFrame,
    min_per_category: int = 30,
    profile: str = DEFAULT_ENDPOINT_PROFILE,
) -> dict:
    failures: list[str] = []
    warnings: list[str] = []
    if profile not in ENDPOINT_PROFILES:
        return {
            "status": "FAIL",
            "n_rows": int(len(df)),
            "endpoint_profile": str(profile),
            "failures": [
                f"Unknown endpoint profile {profile!r}; expected one of "
                f"{sorted(ENDPOINT_PROFILES)}."
            ],
            "warnings": warnings,
        }
    required_columns = ENDPOINT_PROFILES[profile]
    interprets_categories = profile == "categorical_sir"
    missing_columns = sorted(set(required_columns).difference(df.columns))
    if missing_columns:
        failures.append(f"Missing required columns: {missing_columns}")
        return {
            "status": "FAIL",
            "n_rows": int(len(df)),
            "endpoint_profile": profile,
            "failures": failures,
            "warnings": warnings,
        }
    if df.empty:
        failures.append("Cohort table has no rows.")

    # Required columns are all present here (checked above).
    raw_reads = df["genome_source"].astype(str).eq("raw_reads")
    blank_counts = {}
    for column in required_columns:
        blank = blank_mask(df[column])
        if column in BLANK_ALLOWED_FOR_RAW_READS:
            blank = blank & ~raw_reads
        blank_counts[column] = int(blank.sum())
    absent_by_source = {}
    if not interprets_categories:
        absent_by_source = {
            key: blank_counts.pop(key) for key in sorted(ABSENT_ALLOWED_FOR_MIC) if blank_counts.get(key)
        }
    if absent_by_source:
        warnings.append(
            "Absent in source, kept for lineage-grouped evaluation and excluded only from "
            f"analyses that need the field: {absent_by_source}"
        )
    nonzero_blanks = {key: value for key, value in blank_counts.items() if value}
    if nonzero_blanks:
        failures.append(f"Required values are blank: {nonzero_blanks}")

    registered = df["genome_source"].astype(str).eq("registered_assembly")
    invalid_accessions = registered & ~df["assembly_accession"].astype(str).str.match(ASSEMBLY_RE)
    if invalid_accessions.any():
        failures.append(f"Invalid assembly accessions: {int(invalid_accessions.sum())}")
    invalid_sources = ~df["genome_source"].astype(str).isin(GENOME_SOURCES)
    if invalid_sources.any():
        failures.append(f"Invalid genome sources: {sorted(map(str, df.loc[invalid_sources, 'genome_source'].unique()))}")
    source_accession = df["genome_source_accession"].astype(str)
    registered_mismatch = registered & source_accession.ne(df["assembly_accession"].astype(str))
    if registered_mismatch.any():
        failures.append(f"Registered genomes whose source accession differs from the assembly accession: {int(registered_mismatch.sum())}")
    invalid_runs = raw_reads & ~source_accession.str.match(RUN_RE)
    if invalid_runs.any():
        failures.append(f"Raw-read genomes without a valid run accession: {int(invalid_runs.sum())}")
    invalid_genome_hashes = ~df["genome_sha256"].astype(str).str.lower().str.fullmatch(r"[0-9a-f]{64}")
    if invalid_genome_hashes.any():
        failures.append(f"Invalid genome SHA-256 values: {int(invalid_genome_hashes.sum())}")

    invalid_biosamples = ~df["biosample_accession"].astype(str).str.match(BIOSAMPLE_RE)
    if invalid_biosamples.any():
        failures.append(f"Invalid BioSample accessions: {int(invalid_biosamples.sum())}")

    normalized_category = None
    if interprets_categories:
        normalized_category = df["ast_category"].astype(str).str.upper()
        invalid_categories = ~df["ast_category"].astype(str).isin(ALLOWED_CATEGORIES)
        if invalid_categories.any():
            values = sorted(map(str, df.loc[invalid_categories, "ast_category"].unique()))
            failures.append(f"Invalid AST categories: {values}")

    normalized_measurement = df["ast_measurement_type"].astype(str).str.upper()
    invalid_measurements = ~normalized_measurement.isin(ALLOWED_MEASUREMENTS)
    if invalid_measurements.any():
        values = sorted(map(str, df.loc[invalid_measurements, "ast_measurement_type"].unique()))
        failures.append(
            f"Invalid AST measurement types: {values}"
        )

    invalid_signs = ~df["measurement_sign"].astype(str).str.strip().isin(ALLOWED_SIGNS)
    if invalid_signs.any():
        values = sorted(map(str, df.loc[invalid_signs, "measurement_sign"].unique()))
        failures.append(f"Invalid AST measurement comparators: {values}")

    normalized_standard = None
    if interprets_categories:
        normalized_standard = df["breakpoint_standard"].astype(str).str.upper()
        invalid_standards = ~normalized_standard.isin(ALLOWED_STANDARDS)
        if invalid_standards.any():
            values = sorted(map(str, df.loc[invalid_standards, "breakpoint_standard"].unique()))
            failures.append(
                f"Invalid breakpoint standards: {values}"
            )

    invalid_populations = ~df["intended_use_population"].astype(str).isin(ALLOWED_POPULATIONS)
    if invalid_populations.any():
        values = sorted(map(str, df.loc[invalid_populations, "intended_use_population"].unique()))
        failures.append(f"Invalid intended-use populations: {values}")

    invalid_splits = ~df["evaluation_split"].astype(str).isin(ALLOWED_SPLITS)
    if invalid_splits.any():
        values = sorted(map(str, df.loc[invalid_splits, "evaluation_split"].unique()))
        failures.append(f"Invalid evaluation splits: {values}")

    numeric_ast = pd.to_numeric(df["ast_value"], errors="coerce")
    if numeric_ast.isna().any():
        failures.append(f"Non-numeric AST measurements: {int(numeric_ast.isna().sum())}")
    if (numeric_ast <= 0).any():
        failures.append(f"Non-positive AST measurements: {int((numeric_ast <= 0).sum())}")

    expected_units = normalized_measurement.map({"MIC": "mg/L", "ZONE": "mm"})
    invalid_units = df["ast_unit"].astype(str).ne(expected_units)
    if invalid_units.any():
        failures.append(f"AST measurement type/unit mismatches: {int(invalid_units.sum())}")

    zone_rows = normalized_measurement.eq("ZONE")
    if zone_rows.any():
        if "disk_content" not in df.columns:
            failures.append("Disk-zone rows require a disk_content column.")
        elif blank_mask(df.loc[zone_rows, "disk_content"]).any():
            failures.append(
                f"Disk-zone rows with blank disk content: {int(blank_mask(df.loc[zone_rows, 'disk_content']).sum())}"
            )

    present_dates = df["collection_date"][~blank_mask(df["collection_date"])]
    invalid_dates = ~present_dates.map(valid_partial_date).astype(bool)
    if invalid_dates.any():
        examples = sorted(map(str, present_dates[invalid_dates].unique()))[:5]
        failures.append(f"Invalid collection dates: {int(invalid_dates.sum())} (e.g. {examples})")
    if interprets_categories:
        parsed_testing_dates = pd.to_datetime(df["ast_testing_date"], errors="coerce", utc=True)
        if parsed_testing_dates.isna().any():
            failures.append(f"Invalid AST testing dates: {int(parsed_testing_dates.isna().sum())}")

    source_rows = pd.to_numeric(df["source_row_number"], errors="coerce")
    if source_rows.isna().any() or (source_rows < 1).any():
        failures.append("source_row_number must contain positive integers.")

    hash_columns = ["source_fingerprint"]
    if interprets_categories:
        hash_columns += ["breakpoint_rule_sha256", "breakpoint_artifact_sha256"]
    for hash_column in hash_columns:
        invalid_hashes = ~df[hash_column].astype(str).str.lower().str.fullmatch(r"[0-9a-f]{64}")
        if invalid_hashes.any():
            failures.append(f"Invalid SHA-256 values in {hash_column}: {int(invalid_hashes.sum())}")

    key = ["isolate_id", "antibiotic"]
    duplicate_keys = df.duplicated(key, keep=False)
    if duplicate_keys.any():
        failures.append(
            f"Duplicate isolate-antibiotic rows require source-level reconciliation: {int(duplicate_keys.sum())}"
        )

    duplicate_source_ids = df["source_ast_record_id"].astype(str).duplicated(keep=False)
    if duplicate_source_ids.any():
        failures.append(f"Duplicate source AST record identifiers: {int(duplicate_source_ids.sum())}")

    for identity_column in ("biosample_accession", "assembly_accession"):
        counts_by_isolate = df.groupby("isolate_id")[identity_column].nunique()
        if (counts_by_isolate > 1).any():
            failures.append(
                f"Isolates map to multiple {identity_column} values: {int((counts_by_isolate > 1).sum())}"
            )
        isolates_by_identity = df.groupby(identity_column)["isolate_id"].nunique()
        if (isolates_by_identity > 1).any():
            failures.append(
                f"{identity_column} values map to multiple isolates: {int((isolates_by_identity > 1).sum())}"
            )

    cluster_splits = df.groupby("genomic_cluster")["evaluation_split"].nunique()
    crossing_clusters = cluster_splits[cluster_splits > 1]
    if not crossing_clusters.empty:
        failures.append(
            "Near-neighbor genomic clusters cross development/external boundaries: "
            + str(sorted(map(str, crossing_clusters.index))[:20])
        )

    dedup_splits = df.groupby("deduplication_group")["evaluation_split"].nunique()
    if (dedup_splits > 1).any():
        failures.append(
            "Technical/patient/outbreak deduplication groups cross evaluation boundaries: "
            + str(sorted(map(str, dedup_splits[dedup_splits > 1].index))[:20])
        )

    provenance_columns = [
        column
        for column in ("site", "bioproject_accession", "surveillance_network")
        if column in df.columns
    ]
    unavailable_provenance = [
        column
        for column in ("site", "surveillance_network")
        if column not in provenance_columns
    ]
    if unavailable_provenance:
        # Absent source metadata stays absent. Record the reduced separation
        # guarantee instead of substituting a proxy for it.
        warnings.append(
            "External/development separation could not be checked on "
            f"{unavailable_provenance}; these fields are absent from the cohort "
            "and must not be proxied by country, BioProject, or collection date."
        )
    for provenance_column in provenance_columns:
        development_values = set(
            df.loc[df["evaluation_split"].eq("development"), provenance_column].astype(str)
        )
        external_values = set(
            df.loc[df["evaluation_split"].eq("external"), provenance_column].astype(str)
        )
        overlap = sorted(development_values.intersection(external_values))
        if overlap:
            failures.append(
                f"External cohort overlaps development {provenance_column}: {overlap[:20]}"
            )

    if interprets_categories:
        submitted = df["submitted_ast_category"].astype(str).str.upper()
        comparable = submitted.isin(ALLOWED_CATEGORIES)
        discordant = comparable & submitted.ne(normalized_category)
        if discordant.any():
            warnings.append(
                f"Submitted and independently derived AST categories disagree for {int(discordant.sum())} rows."
            )

        standard_keys = normalized_standard + " " + df["breakpoint_version"].astype(str).str.strip()
        standard_counts = pd.DataFrame(
            {"antibiotic": df["antibiotic"].astype(str), "standard": standard_keys}
        ).groupby("antibiotic")["standard"].nunique()
        mixed = standard_counts[standard_counts > 1]
        if not mixed.empty:
            failures.append(f"Mixed breakpoint standards/versions within antibiotics: {mixed.to_dict()}")

    split_sets = df.groupby("antibiotic")["evaluation_split"].agg(lambda values: set(map(str, values)))
    missing_external = [drug for drug, splits in split_sets.items() if splits != ALLOWED_SPLITS]
    if missing_external:
        failures.append(f"Antibiotics without both development and external cohorts: {missing_external}")

    censoring = None
    if interprets_categories:
        category_sets = df.groupby(["antibiotic", "evaluation_split"])["ast_category"].agg(set)
        missing_binary_support = [
            f"{drug}/{split}"
            for (drug, split), categories in category_sets.items()
            if not {"S", "R"}.issubset(categories)
        ]
        if missing_binary_support:
            failures.append(
                "Drug/partition groups without both S and R support: " + str(missing_binary_support)
            )

        counts = (
            df.assign(ast_category=normalized_category)
            .groupby(["antibiotic", "evaluation_split", "ast_category"])
            .size()
            .rename("n")
            .reset_index()
        )
    else:
        # Under an interval-censored MIC endpoint there is no category to count.
        # Report measurement support and the censoring profile instead, because
        # a drug whose measurements are almost entirely censored at one panel
        # limit carries little quantitative signal.
        counts = (
            df.groupby(["antibiotic", "evaluation_split"])
            .size()
            .rename("n")
            .reset_index()
        )
        censored_mask = df["measurement_sign"].astype(str).str.strip().isin({"<", "<=", ">", ">="})
        censoring = (
            df.assign(censored=censored_mask)
            .groupby(["antibiotic", "evaluation_split"])["censored"]
            .agg(["size", "sum"])
            .rename(columns={"size": "n", "sum": "n_censored"})
            .reset_index()
        )
        censoring["fraction_censored"] = (
            censoring["n_censored"].astype(float) / censoring["n"].astype(float)
        ).round(4)
        fully_censored = censoring[censoring["fraction_censored"] >= 1.0]
        if not fully_censored.empty:
            failures.append(
                "Drug/partition groups with no exact MIC measurement: "
                + str(fully_censored[["antibiotic", "evaluation_split"]].to_dict(orient="records"))
            )
        censoring = censoring.to_dict(orient="records")
        # The external set is sealed until models and thresholds are frozen.
        # Its censoring fraction is close to a resistance prevalence (a '>'
        # MIC is usually above the breakpoint), so it is checked above but
        # never reported.
        for record in censoring:
            if str(record["evaluation_split"]) == "external":
                record["n_censored"] = "sealed"
                record["fraction_censored"] = "sealed"

    sparse = counts[counts["n"] < min_per_category]
    if not sparse.empty:
        warnings.append(
            "Counts below the screening threshold; perform a formal precision/power calculation: "
            + str(sparse.to_dict(orient="records"))
        )

    report = {
        "status": "PASS" if not failures else "FAIL",
        "n_rows": int(len(df)),
        "n_isolates": int(df["isolate_id"].nunique()),
        "endpoint_profile": profile,
        "antibiotics": sorted(map(str, df["antibiotic"].unique())),
        "counts": counts.to_dict(orient="records"),
        "failures": failures,
        "warnings": warnings,
    }
    if censoring is not None:
        report["censoring"] = censoring
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--report", type=Path, default=Path("reports/cohort_validation.json"))
    parser.add_argument("--study-config", type=Path, help="Freeze and verify the breakpoint artifact")
    parser.add_argument("--breakpoint-rules", type=Path, help="Machine-readable frozen breakpoint rules")
    parser.add_argument("--breakpoint-artifact", type=Path, help="Original frozen breakpoint artifact")
    parser.add_argument("--min-per-category", type=int, default=30)
    parser.add_argument(
        "--endpoint-profile",
        choices=sorted(ENDPOINT_PROFILES),
        help=(
            "Endpoint contract to enforce. Defaults to the study config's "
            "primary_endpoint, then to the strict categorical contract. "
            "See docs/ENDPOINT_AMENDMENT.md."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cohort = pd.read_csv(args.cohort, low_memory=False)
    study_config = (
        json.loads(args.study_config.read_text(encoding="utf-8")) if args.study_config else None
    )
    profile = resolve_endpoint_profile(study_config, args.endpoint_profile)
    report = validate_cohort(
        cohort, min_per_category=args.min_per_category, profile=profile
    )
    if study_config is not None:
        breakpoint_failures = validate_breakpoint_artifact(
            cohort,
            study_config,
            Path.cwd(),
            profile=profile,
        )
        report["failures"].extend(breakpoint_failures)
        report["study_config"] = str(args.study_config)
        if profile == "categorical_sir":
            ast_config = study_config.get("ast", {})
            configured_rules = ast_config.get("breakpoint_rules")
            configured_artifact = ast_config.get("breakpoint_table")
            rules_path = args.breakpoint_rules or (
                Path(configured_rules) if configured_rules else None
            )
            artifact_path = args.breakpoint_artifact or (
                Path(configured_artifact) if configured_artifact else None
            )
            if rules_path is None or artifact_path is None:
                report["failures"].append(
                    "Study config validation requires machine-readable breakpoint rules and the source artifact."
                )
            else:
                report["failures"].extend(
                    validate_breakpoint_derivation(cohort, rules_path, artifact_path)
                )
        report["status"] = "PASS" if not report["failures"] else "FAIL"
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Cohort validation: {report['status']} | rows={report['n_rows']}")
    for warning in report["warnings"]:
        print(f"WARNING: {warning}")
    for failure in report["failures"]:
        print(f"FAILURE: {failure}")
    print(f"Report: {args.report}")
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
