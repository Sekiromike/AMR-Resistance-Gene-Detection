"""No-genome baseline for the MIC endpoint (docs/MIC_EVALUATION_SPEC.md).

Per antibiotic and per development fold, predict for every held-out row the
single dilution that maximises essential agreement on the training folds
(ties: higher exact-only agreement, then the lower dilution). This is the
strongest constant prediction under the primary metric, so a model credited
against it has learned something from the genome.

The originally specified intercept-only interval-censored normal was replaced
on 2026-09-27 (decision owner) after its first run: development MICs are
bimodal and heavily censored at both panel limits, and the single normal fitted
degenerately (e.g. ceftriaxone mean 2^-24 mg/L, sd 30 log2). No genomic model
had been fitted and no external MIC read.

Reference MICs are read with the evaluator's rules: a reported dilution d
means the latent log2 MIC is in (d-1, d]; comparators open the interval.

Development only. The external set is never read here.

Run from the repository root:
    python scripts/run_mic_baseline.py --cohort ... --splits ... \\
        --population human_clinical --output ... --manifest ...
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from scripts.evaluate_mic_predictions import reference_bounds
except ModuleNotFoundError:  # Direct execution
    from evaluate_mic_predictions import reference_bounds  # type: ignore[no-redef]

POPULATIONS = ("human_clinical", "all_sources")


class BaselineError(ValueError):
    """Raised when the cohort or splits violate the baseline contract."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def best_constant_dilution(lo: np.ndarray, hi: np.ndarray) -> tuple[int, float, float]:
    """Dilution maximising training EA; returns (dilution, EA, exact-only EA)."""
    finite = np.concatenate([lo[np.isfinite(lo)], hi[np.isfinite(hi)]])
    if finite.size == 0:
        raise BaselineError("Cannot fit a baseline: every training interval is unbounded")
    exact = lo == hi
    best: tuple[float, float, int] | None = None
    for candidate in range(int(finite.min()) - 1, int(finite.max()) + 2):
        error = np.maximum(np.maximum(candidate - hi, lo - candidate), 0.0)
        agree = error <= 1
        ea = float(agree.mean())
        ea_exact = float(agree[exact].mean()) if exact.any() else 0.0
        key = (ea, ea_exact, -candidate)
        if best is None or key > best:
            best = key
    assert best is not None
    return -best[2], best[0], best[1]


def run_baseline(cohort: pd.DataFrame, splits: pd.DataFrame, population: str) -> tuple[pd.DataFrame, list[dict]]:
    if population not in POPULATIONS:
        raise BaselineError(f"population must be one of {POPULATIONS}")
    development = cohort.loc[cohort["evaluation_split"].astype(str).eq("development")].copy()
    if population == "human_clinical":
        development = development.loc[development["intended_use_population"].astype(str).eq("human_clinical")]
    folds = splits.loc[splits["evaluation_split"].astype(str).eq("development"), ["isolate_id", "cv_fold"]]
    development = development.merge(folds, on="isolate_id", how="left", validate="many_to_one")
    if development["cv_fold"].isna().any():
        raise BaselineError("Every development isolate needs a cv_fold")
    development["cv_fold"] = development["cv_fold"].astype(int)
    bounds = [reference_bounds(str(s), float(v))
              for s, v in zip(development["measurement_sign"], development["ast_value"])]
    development["_lo"] = [lo for lo, _ in bounds]
    development["_hi"] = [hi for _, hi in bounds]

    predictions: list[pd.DataFrame] = []
    fits: list[dict[str, Any]] = []
    for antibiotic, drug in development.groupby("antibiotic", sort=True):
        for fold in sorted(drug["cv_fold"].unique()):
            train = drug.loc[drug["cv_fold"] != fold]
            test = drug.loc[drug["cv_fold"] == fold]
            dilution, ea, ea_exact = best_constant_dilution(train["_lo"].to_numpy(float),
                                                            train["_hi"].to_numpy(float))
            fits.append({"antibiotic": str(antibiotic), "cv_fold": int(fold), "n_train": int(len(train)),
                         "n_test": int(len(test)), "dilution_log2": dilution, "mic_mg_per_l": 2.0 ** dilution,
                         "training_essential_agreement": ea, "training_essential_agreement_exact_only": ea_exact})
            predictions.append(test.assign(predicted_log2_mic=float(dilution)))
    columns = ["isolate_id", "antibiotic", "evaluation_split", "lineage_group", "intended_use_population",
               "cv_fold", "measurement_sign", "ast_value", "predicted_log2_mic"]
    output = pd.concat(predictions, ignore_index=True)[columns]
    return output.sort_values(["antibiotic", "isolate_id"]).reset_index(drop=True), fits


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--population", choices=POPULATIONS, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cohort = pd.read_csv(args.cohort, dtype={"isolate_id": str, "lineage_group": str}, low_memory=False)
    splits = pd.read_csv(args.splits, dtype={"isolate_id": str})
    try:
        output, fits = run_baseline(cohort, splits, args.population)
    except BaselineError as exc:
        raise SystemExit(f"MIC baseline failed: {exc}") from exc
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    payload = {
        "schema_version": "1.0.0",
        "operation": "run_mic_baseline",
        "model": "best constant dilution per antibiotic under training-fold essential agreement",
        "generated_at_utc": utc_now(),
        "population": args.population,
        "inputs": {"cohort": {"path": str(args.cohort), "sha256": sha256_file(args.cohort)},
                   "splits": {"path": str(args.splits), "sha256": sha256_file(args.splits)}},
        "fits": fits,
        "outputs": {"predictions": {"path": str(args.output), "sha256": sha256_file(args.output),
                                    "rows": int(len(output))}},
        "scientific_boundary": {"external_rows_read": False, "genome_features_used": False},
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
