"""Join external reference MICs to LOCKED predictions (seal step 3).

docs/EXTERNAL_EVALUATION_PLAN.md. Refuses unless the prediction file's
SHA-256 matches the entry recorded in the committed lock file
(config/external_prediction_lock.tsv) before any external MIC was read, and
unless every external prediction row has exactly one reference row.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

import pandas as pd

REFERENCE_COLUMNS = ["measurement_sign", "ast_value"]


class LockError(ValueError):
    """Raised when predictions are not the locked ones or rows do not match."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def locked_hashes(lock: Path) -> dict[str, str]:
    with lock.open("r", encoding="utf-8", newline="") as handle:
        return {row["file"]: row["sha256"] for row in csv.DictReader(handle, delimiter="\t")}


def attach(predictions: Path, cohort: Path, lock: Path) -> pd.DataFrame:
    recorded = locked_hashes(lock).get(predictions.name)
    if recorded is None:
        raise LockError(f"{predictions.name} is not in the lock file")
    if sha256_file(predictions) != recorded:
        raise LockError(f"{predictions.name} differs from its locked SHA-256")
    frame = pd.read_csv(predictions, dtype={"isolate_id": str, "lineage_group": str})
    if set(REFERENCE_COLUMNS) & set(frame.columns):
        raise LockError("Locked predictions must not carry reference columns")
    references = pd.read_csv(cohort, dtype=str, keep_default_na=False, low_memory=False,
                             usecols=["isolate_id", "antibiotic", "evaluation_split", *REFERENCE_COLUMNS])
    references = references[references["evaluation_split"].eq("external")].drop(columns="evaluation_split")
    merged = frame.merge(references, on=["isolate_id", "antibiotic"], how="left", validate="one_to_one",
                         indicator=True)
    if not merged["_merge"].eq("both").all():
        raise LockError(f"{int((merged['_merge'] != 'both').sum())} predictions have no external reference")
    if len(references) != len(frame):
        raise LockError(f"{len(references)} external reference rows but {len(frame)} predictions")
    merged = merged.drop(columns="_merge")
    merged["ast_value"] = merged["ast_value"].astype(float)
    return merged


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--lock", type=Path, default=Path("config/external_prediction_lock.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        merged = attach(args.predictions, args.cohort, args.lock)
    except LockError as exc:
        raise SystemExit(f"Refusing to attach references: {exc}") from exc
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
