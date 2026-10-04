"""Normalize the frozen JARBS external MIC table into the cohort AST schema.

The external cohort stays sealed: this script reshapes rows so they can enter
the same cohort table as development data, and reports only structural counts
(rows, missingness, eligibility). It never reports MIC values, distributions,
or resistance prevalence, and nothing here is used to choose features,
models, thresholds or epochs.

Study-level facts are recorded as such, not as per-isolate measurements:
MICs were centrally remeasured by broth microdilution on MicroScan WalkAway
panels; the per-isolate panel and the AST testing date are not published and
stay blank. Country comes from the same INSDC geo_loc_name rule as development
(amendment 006).

Run from the repository root:
    python scripts/normalize_jarbs_ast.py --jarbs ... --output ... --manifest ...
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

try:
    from scripts.normalize_ncbi_ast import derive_country, load_country_vocabulary
except ModuleNotFoundError:  # Direct execution
    from normalize_ncbi_ast import derive_country, load_country_vocabulary  # type: ignore[no-redef]

SOURCE_DATASET = "JARBS-GNR Supplementary Data 6"
BIOPROJECT = "PRJDB10842"
DRUGS = ("ceftriaxone", "ciprofloxacin", "gentamicin")
MIC_PATTERN = re.compile(r"^(<=|>|<|>=)?\s*(\d+(?:\.\d+)?)$")
OUTPUT_FIELDS = [
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
    "raw_measurement_sign",
    "measurement_sign",
    "raw_ast_value",
    "ast_value",
    "ast_unit",
    "ast_method",
    "ast_method_basis",
    "ast_platform",
    "ast_testing_date",
    "mic_eligible",
    "mic_ineligibility_reasons",
    "collection_date",
    "geo_loc_name",
    "country",
    "country_source",
    "bioproject_accession",
    "isolation_source",
]


class JarbsNormalizationError(ValueError):
    """Raised when the frozen table does not match the expected structure."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_mic(raw: str) -> tuple[str, str, str]:
    """Return (sign, value, reason). An empty reason means eligible."""
    text = raw.strip()
    match = MIC_PATTERN.match(text)
    if not match:
        return "", "", "unparseable_mic"
    sign = match.group(1) or "="
    value = match.group(2)
    if float(value) <= 0:
        return sign, value, "non_positive_mic"
    return sign, value, ""


def normalize_jarbs(
    jarbs_path: Path,
    vocabulary_path: Path,
    source_release: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    table = pd.read_csv(jarbs_path, dtype=str, keep_default_na=False)
    required = {"isolate_id", "genome_species", "biosample_accession", "collection_date_raw",
                "geo_loc_name_raw", "isolation_source_raw"} | {f"{d}_mic_raw" for d in DRUGS}
    missing = sorted(required.difference(table.columns))
    if missing:
        raise JarbsNormalizationError(f"JARBS table is missing columns: {missing}")
    if table["isolate_id"].duplicated().any():
        raise JarbsNormalizationError("JARBS table has duplicate isolate identifiers")
    species = {value.strip() for value in table["genome_species"]}
    if not species <= {"Escherichia coli", "E. coli"}:
        raise JarbsNormalizationError(f"Unexpected genome_species values: {sorted(species)}")

    fingerprint = sha256_file(jarbs_path)
    rows: list[dict[str, str]] = []
    not_tested = {drug: 0 for drug in DRUGS}
    for index, record in table.iterrows():
        date = record["collection_date_raw"].strip()
        date = "" if date.casefold() in {"missing", ""} else date
        for drug in DRUGS:
            raw = record[f"{drug}_mic_raw"].strip()
            if not raw:
                not_tested[drug] += 1
                continue
            sign, value, reason = parse_mic(raw)
            rows.append({
                "source_dataset": SOURCE_DATASET,
                "source_release": source_release,
                "source_fingerprint": fingerprint,
                "source_row_number": str(int(index) + 2),  # header is line 1
                "source_ast_record_id": f"JARBS:{record['isolate_id']}:{drug}",
                "isolate_id": record["isolate_id"],
                "biosample_accession": record["biosample_accession"].strip(),
                "assembly_accession": "",  # assembled here from raw reads
                "scientific_name": "Escherichia coli",
                "antibiotic": drug,
                "ast_measurement_type": "MIC",
                "raw_measurement_sign": sign,
                "measurement_sign": sign,
                "raw_ast_value": raw,
                "ast_value": value,
                "ast_unit": "mg/L",
                "ast_method": "broth microdilution",
                "ast_method_basis": "study-level documentation (central remeasurement, MicroScan WalkAway panels)",
                "ast_platform": "MicroScan WalkAway",
                "ast_testing_date": "",  # not published at isolate level
                "mic_eligible": "true" if not reason else "false",
                "mic_ineligibility_reasons": reason,
                "collection_date": date,
                "geo_loc_name": record["geo_loc_name_raw"].strip(),
                "country": "",
                "bioproject_accession": BIOPROJECT,
                "isolation_source": record["isolation_source_raw"].strip(),
            })
    frame = pd.DataFrame(rows, columns=[f for f in OUTPUT_FIELDS if f != "country_source"])
    frame = derive_country(frame, load_country_vocabulary(vocabulary_path))
    frame = frame[OUTPUT_FIELDS]

    counts = {
        "isolates": int(len(table)),
        "rows": int(len(frame)),
        "rows_by_drug": {drug: int((frame["antibiotic"] == drug).sum()) for drug in DRUGS},
        "not_tested_by_drug": not_tested,
        "mic_eligible_rows": int((frame["mic_eligible"] == "true").sum()),
        "ineligibility_reasons": {
            str(k): int(v) for k, v in frame.loc[frame["mic_ineligibility_reasons"] != "",
                                                  "mic_ineligibility_reasons"].value_counts().items()
        },
        "country_source": {str(k): int(v) for k, v in frame["country_source"].value_counts().items()},
        "collection_date_absent_rows": int((frame["collection_date"] == "").sum()),
    }
    return frame, counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--jarbs", type=Path, required=True)
    parser.add_argument("--country-vocabulary", type=Path, default=Path("config/insdc_geo_loc_name_vocabulary.tsv"))
    parser.add_argument("--source-release", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    try:
        frame, counts = normalize_jarbs(args.jarbs, args.country_vocabulary, args.source_release)
    except JarbsNormalizationError as exc:
        raise SystemExit(f"JARBS normalization failed: {exc}") from exc
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, index=False, lineterminator="\n")
    payload = {
        "schema_version": "1.0.0",
        "operation": "normalize_jarbs_ast",
        "generated_at_utc": utc_now(),
        "inputs": {
            "jarbs": {"path": str(args.jarbs), "sha256": sha256_file(args.jarbs)},
            "country_vocabulary": {"path": str(args.country_vocabulary), "sha256": sha256_file(args.country_vocabulary)},
        },
        "counts": counts,
        "outputs": {"normalized": {"path": str(args.output), "sha256": sha256_file(args.output)}},
        "scientific_boundary": {
            "external_cohort_sealed": True,
            "values_reported": "structural counts only; no MIC values, distributions or prevalence",
        },
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(counts, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
