"""Derive per-isolate specimen source and intended-use population.

Amendment 006. The primary analysis population is human clinical isolates; the
all-source cohort is a prespecified sensitivity analysis. Classification uses
only source metadata recorded before any phenotype is modelled, and never
guesses: an isolate whose evidence is incomplete or conflicting is
`undetermined`, not assigned.

Development (NCBI Pathogen Detection), per isolate:

| epi_type      | host         | population                   |
|---------------|--------------|------------------------------|
| clinical      | Homo sapiens | human_clinical               |
| clinical      | absent       | undetermined                 |
| clinical      | non-human    | non_human_or_environmental   |
| env./other    | Homo sapiens | undetermined (conflict)      |
| env./other    | other/absent | non_human_or_environmental   |
| absent        | any          | undetermined                 |

NCBI's epi_type "clinical" alone is not sufficient: in the frozen snapshot 259
clinical isolates have animal hosts (veterinary clinical), so a human host is
required.

External (JARBS): human_clinical by documented study design -- hospital
clinical-isolate surveillance -- recorded as the basis, not inferred.

Run from the repository root:
    python scripts/build_source_attributes.py --development-ast ... --jarbs ... --output ...
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

OUTPUT_FIELDS = [
    "isolate_id",
    "cohort",
    "specimen_source",
    "intended_use_population",
    "population_basis",
]
POPULATIONS = ("human_clinical", "non_human_or_environmental", "undetermined")
HUMAN_HOSTS = frozenset({"homo sapiens"})
NULL_LIKE = frozenset({"", "nan", "none", "null", "<na>", "missing", "not collected",
                       "not applicable", "not provided", "unknown"})
JARBS_BASIS = (
    "study design: JARBS-GNR hospital clinical-isolate surveillance "
    "(Supplementary Data 6); not inferred from specimen"
)


class AttributeError_(ValueError):
    """Raised when source attributes cannot be derived consistently."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _clean(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)) or value is pd.NA:
        return ""
    text = " ".join(str(value).split())
    return "" if text.casefold() in NULL_LIKE else text


def classify_development(epi_type: str, host: str) -> tuple[str, str]:
    epi = epi_type.casefold()
    human = host.casefold() in HUMAN_HOSTS
    if not epi:
        return "undetermined", "epi_type absent"
    if epi == "clinical":
        if human:
            return "human_clinical", "epi_type=clinical and host=Homo sapiens"
        if not host:
            return "undetermined", "epi_type=clinical but host absent"
        return "non_human_or_environmental", f"epi_type=clinical with non-human host ({host})"
    if human:
        return "undetermined", f"host=Homo sapiens conflicts with epi_type={epi_type}"
    return "non_human_or_environmental", f"epi_type={epi_type}"


def build_source_attributes(
    development_ast_path: Path,
    jarbs_path: Path,
    output_path: Path,
    manifest_path: Path,
) -> dict[str, Any]:
    development = pd.read_csv(development_ast_path, dtype=str, keep_default_na=False)
    required = {"isolate_id", "epi_type", "host", "isolation_source"}
    missing = sorted(required.difference(development.columns))
    if missing:
        raise AttributeError_(f"Development AST is missing columns: {missing}")
    rows: list[dict[str, str]] = []
    for isolate_id, group in development.groupby("isolate_id", sort=True):
        values = {}
        for field in ("epi_type", "host", "isolation_source"):
            distinct = {_clean(v) for v in group[field]}
            if len(distinct) > 1:
                raise AttributeError_(f"{isolate_id} has conflicting {field} values: {sorted(distinct)}")
            values[field] = distinct.pop()
        population, basis = classify_development(values["epi_type"], values["host"])
        rows.append(
            {
                "isolate_id": isolate_id,
                "cohort": "development",
                "specimen_source": values["isolation_source"],
                "intended_use_population": population,
                "population_basis": basis,
            }
        )

    jarbs = pd.read_csv(jarbs_path, dtype=str, keep_default_na=False)
    for field in ("isolate_id", "isolation_source_raw"):
        if field not in jarbs.columns:
            raise AttributeError_(f"JARBS table is missing column: {field}")
    if jarbs["isolate_id"].duplicated().any():
        raise AttributeError_("JARBS table has duplicate isolate identifiers")
    for record in jarbs.sort_values("isolate_id").itertuples(index=False):
        rows.append(
            {
                "isolate_id": record.isolate_id,
                "cohort": "external",
                "specimen_source": _clean(record.isolation_source_raw),
                "intended_use_population": "human_clinical",
                "population_basis": JARBS_BASIS,
            }
        )

    identifiers = [row["isolate_id"] for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise AttributeError_("An isolate identifier appears in both cohorts")

    frame = pd.DataFrame(rows, columns=OUTPUT_FIELDS)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, sep="\t", index=False, lineterminator="\n")

    counts = {
        cohort: dict(Counter(group["intended_use_population"]))
        for cohort, group in frame.groupby("cohort")
    }
    specimen_absent = {
        cohort: int((group["specimen_source"] == "").sum())
        for cohort, group in frame.groupby("cohort")
    }
    payload = {
        "schema_version": "1.0.0",
        "operation": "build_source_attributes",
        "amendment": "docs/COUNTRY_AND_POPULATION_AMENDMENT.md",
        "generated_at_utc": utc_now(),
        "inputs": {
            "development_ast": {"path": str(development_ast_path), "sha256": sha256_file(development_ast_path)},
            "jarbs": {"path": str(jarbs_path), "sha256": sha256_file(jarbs_path)},
        },
        "population_counts": counts,
        "specimen_source_absent": specimen_absent,
        "outputs": {"attributes": {"path": str(output_path), "sha256": sha256_file(output_path)}},
        "scientific_boundary": {
            "phenotypes_read": False,
            "primary_population": "human_clinical",
            "sensitivity_population": "all_sources",
        },
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--development-ast", type=Path, required=True)
    parser.add_argument("--jarbs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        payload = build_source_attributes(args.development_ast, args.jarbs, args.output, args.manifest)
    except AttributeError_ as exc:
        raise SystemExit(f"Source attribute derivation failed: {exc}") from exc
    print(json.dumps({"population_counts": payload["population_counts"],
                      "specimen_source_absent": payload["specimen_source_absent"]}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
