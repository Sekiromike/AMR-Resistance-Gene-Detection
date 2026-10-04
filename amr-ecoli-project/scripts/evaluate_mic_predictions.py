"""Evaluate MIC predictions against interval-censored reference MICs.

Primary endpoint (docs/ENDPOINT_AMENDMENT.md): log2 MIC, interval-censored.
Scoring rules (docs/MIC_EVALUATION_SPEC.md):

Reference MIC as a set of doubling dilutions, on the integer log2 scale
(d = log2 mg/L; conventional labels 0.12, 0.06, ... snap to the nearest
dilution within 0.1 log2; anything further off the doubling scale is refused):

    "= y"  -> {y}          "<= y" -> {.., y}      "< y"  -> {.., y-1}
    ">= y" -> {y, ..}      "> y"  -> {y+1, ..}

A point prediction is rounded to the nearest dilution p. Its error is the
number of dilution steps from p to the nearest member of the reference set
(0 when p is consistent with it). Then:

- essential agreement (EA): error <= 1;
- exact agreement: error == 0;
- over-call / under-call by >= 2 dilutions: p above / below the set by >= 2;
- bias: mean signed error p - y on exact references only.

Each metric is reported over all rows and over exact-reference rows only.
Addendum A: every result is also broken down by reference class (left-
censored, exact, right-censored); for right-censored rows (clearly resistant)
essential agreement and the under-call rate carry whole-lineage bootstrap
intervals, the under-call being the MIC analogue of a very major error.

Interval log-likelihood, when a Gaussian predictive sd is supplied: a reported
dilution d means the latent log2 MIC lies in (d-1, d], so the reference set
[lo, hi] is the latent interval (lo-1, hi].

Confidence intervals resample lineage_group (config evaluation.bootstrap_unit)
with replacement. Each drug is also reported for the dominant lineage (config
evaluation.dominant_lineage, amendment 007) and for all other lineages; the
dominant stratum is a single lineage, so its interval resamples isolates and
is labelled as such. External rows are refused unless --locked-external is given,
so the sealed set cannot be scored by accident.

Run from the repository root:
    python scripts/evaluate_mic_predictions.py --predictions ... --output ...
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import norm

SNAP_TOLERANCE = 0.1  # log2 units
BOOTSTRAP_UNIT = "lineage_group"
REQUIRED_COLUMNS = (
    "isolate_id",
    "antibiotic",
    "evaluation_split",
    "lineage_group",
    "measurement_sign",
    "ast_value",
    "predicted_log2_mic",
)


class MicEvaluationError(ValueError):
    """Raised when predictions or references violate the evaluation contract."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dilution_index(value: float) -> int:
    """Return the integer log2 dilution for a positive mg/L value, or raise."""
    if not (isinstance(value, (int, float)) and math.isfinite(value) and value > 0):
        raise MicEvaluationError(f"MIC value must be positive and finite, got {value!r}")
    log2 = math.log2(value)
    nearest = round(log2)
    if abs(log2 - nearest) > SNAP_TOLERANCE:
        raise MicEvaluationError(f"MIC value {value} is not on the doubling-dilution scale")
    return int(nearest)


def reference_bounds(sign: str, value: float) -> tuple[float, float]:
    """Inclusive dilution bounds (lo, hi) of the reference set; +-inf when open."""
    d = dilution_index(value)
    sign = sign.strip()
    if sign in {"=", "=="}:
        return float(d), float(d)
    if sign == "<=":
        return -math.inf, float(d)
    if sign == "<":
        return -math.inf, float(d - 1)
    if sign == ">=":
        return float(d), math.inf
    if sign == ">":
        return float(d + 1), math.inf
    raise MicEvaluationError(f"Unknown MIC comparator {sign!r}")


def predicted_dilution(log2_mic: np.ndarray) -> np.ndarray:
    """Round to the nearest dilution; halves round up (deterministic)."""
    return np.floor(np.asarray(log2_mic, dtype=float) + 0.5)


def interval_log_likelihood(lo: np.ndarray, hi: np.ndarray, mu: np.ndarray, sd: np.ndarray) -> np.ndarray:
    """log P(latent in (lo-1, hi]) under N(mu, sd), computed stably in both tails."""
    a = (np.asarray(lo, dtype=float) - 1.0 - mu) / sd
    b = (np.asarray(hi, dtype=float) - mu) / sd
    upper_tail = a > 0  # both bounds in the upper tail: use survival functions
    prob = np.where(upper_tail, norm.sf(a) - norm.sf(b), norm.cdf(b) - norm.cdf(a))
    return np.log(np.clip(prob, 1e-300, None))


def score_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Add reference bounds, rounded prediction and per-row agreement columns."""
    missing = [column for column in REQUIRED_COLUMNS if column not in frame.columns]
    if missing:
        raise MicEvaluationError(f"Predictions lack columns: {missing}")
    if frame["predicted_log2_mic"].isna().any():
        raise MicEvaluationError("Every row needs a prediction; abstention is reported separately.")
    duplicated = frame.duplicated(["isolate_id", "antibiotic"])
    if duplicated.any():
        raise MicEvaluationError(f"{int(duplicated.sum())} duplicate isolate-antibiotic predictions")
    bounds = [reference_bounds(str(s), float(v)) for s, v in zip(frame["measurement_sign"], frame["ast_value"])]
    out = frame.copy()
    out["ref_lo"] = [lo for lo, _ in bounds]
    out["ref_hi"] = [hi for _, hi in bounds]
    out["exact_reference"] = out["ref_lo"] == out["ref_hi"]
    p = predicted_dilution(out["predicted_log2_mic"].to_numpy())
    out["predicted_dilution"] = p
    above = p - out["ref_hi"].to_numpy()
    below = out["ref_lo"].to_numpy() - p
    out["error_dilutions"] = np.maximum(np.maximum(above, below), 0.0)
    out["essential_agreement"] = out["error_dilutions"] <= 1
    out["exact_agreement"] = out["error_dilutions"] == 0
    out["over_call_2plus"] = above >= 2
    out["under_call_2plus"] = below >= 2
    # A point-prediction model (the baselines) may carry an entirely empty sd
    # column: that means "no predictive distribution", not an error. A
    # partially empty column is still refused.
    if "predicted_log2_sd" in out.columns and out["predicted_log2_sd"].isna().all():
        out = out.drop(columns="predicted_log2_sd")
    if "predicted_log2_sd" in out.columns:
        sd = out["predicted_log2_sd"].to_numpy(dtype=float)
        if not np.all(np.isfinite(sd) & (sd > 0)):
            raise MicEvaluationError("predicted_log2_sd must be positive and finite")
        out["log_likelihood"] = interval_log_likelihood(
            out["ref_lo"].to_numpy(), out["ref_hi"].to_numpy(),
            out["predicted_log2_mic"].to_numpy(dtype=float), sd)
    return out


def _rate(mask: pd.Series) -> float:
    return float(mask.mean()) if len(mask) else float("nan")


def summarize(scored: pd.DataFrame) -> dict[str, Any]:
    exact = scored[scored["exact_reference"]]
    summary: dict[str, Any] = {
        "n": int(len(scored)),
        "n_exact_reference": int(len(exact)),
        "n_censored_reference": int(len(scored) - len(exact)),
        "essential_agreement": _rate(scored["essential_agreement"]),
        "essential_agreement_exact_only": _rate(exact["essential_agreement"]),
        "exact_agreement": _rate(scored["exact_agreement"]),
        "exact_agreement_exact_only": _rate(exact["exact_agreement"]),
        "over_call_2plus": _rate(scored["over_call_2plus"]),
        "under_call_2plus": _rate(scored["under_call_2plus"]),
        "bias_exact_only": (
            float((exact["predicted_dilution"] - exact["ref_lo"]).mean()) if len(exact) else float("nan")
        ),
    }
    if "log_likelihood" in scored.columns:
        summary["mean_interval_log_likelihood"] = float(scored["log_likelihood"].mean())
    return summary


def reference_class(scored: pd.DataFrame) -> pd.Series:
    cls = pd.Series("exact", index=scored.index)
    cls[np.isneginf(scored["ref_lo"].to_numpy(float))] = "left_censored"
    cls[np.isposinf(scored["ref_hi"].to_numpy(float))] = "right_censored"
    return cls


def proportion_bootstrap(indicator: pd.Series, groups: pd.Series, iterations: int, seed: int) -> list[float]:
    """95% interval of a proportion, resampling whole groups."""
    per = pd.DataFrame({"y": indicator.astype(float), "g": groups.astype(str)}).groupby("g")["y"].agg(["sum", "count"])
    if len(per) < 2:
        return [float("nan"), float("nan")]
    sums, counts = per["sum"].to_numpy(), per["count"].to_numpy()
    draws = np.random.default_rng(seed).integers(0, len(per), size=(iterations, len(per)))
    values = sums[draws].sum(axis=1) / counts[draws].sum(axis=1)
    return [float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))]


def by_reference_class(scored: pd.DataFrame, iterations: int, seed: int) -> dict[str, Any]:
    cls = reference_class(scored)
    out: dict[str, Any] = {}
    for name in ("left_censored", "exact", "right_censored"):
        part = scored[cls == name]
        entry: dict[str, Any] = {"n": int(len(part))}
        if len(part):
            entry.update(essential_agreement=_rate(part["essential_agreement"]),
                         under_call_2plus=_rate(part["under_call_2plus"]),
                         over_call_2plus=_rate(part["over_call_2plus"]))
            if name == "right_censored":
                entry["essential_agreement_95ci"] = proportion_bootstrap(
                    part["essential_agreement"], part[BOOTSTRAP_UNIT], iterations, seed)
                entry["under_call_2plus_95ci"] = proportion_bootstrap(
                    part["under_call_2plus"], part[BOOTSTRAP_UNIT], iterations, seed)
        out[name] = entry
    return out


def lineage_bootstrap(
    scored: pd.DataFrame, iterations: int, seed: int, unit: str = BOOTSTRAP_UNIT
) -> dict[str, list[float]]:
    """95% percentile intervals for EA, resampling whole lineages."""
    if scored[unit].astype(str).str.strip().eq("").any():
        raise MicEvaluationError(f"Every row needs a {unit} for the bootstrap")
    per_group = scored.assign(
        _ea=scored["essential_agreement"].astype(float),
        _exact=scored["exact_reference"].astype(float),
        _ea_exact=(scored["essential_agreement"] & scored["exact_reference"]).astype(float),
    ).groupby(scored[unit].astype(str))[["_ea", "_exact", "_ea_exact"]].agg(["sum", "count"])
    n = per_group[("_ea", "count")].to_numpy()
    ea = per_group[("_ea", "sum")].to_numpy()
    n_exact = per_group[("_exact", "sum")].to_numpy()
    ea_exact = per_group[("_ea_exact", "sum")].to_numpy()
    groups = len(n)
    if groups < 2:
        raise MicEvaluationError(f"Bootstrap needs at least two {unit} values")
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, groups, size=(iterations, groups))
    ea_all = ea[draws].sum(axis=1) / n[draws].sum(axis=1)
    exact_den = n_exact[draws].sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        ea_ex = np.where(exact_den > 0, ea_exact[draws].sum(axis=1) / exact_den, np.nan)

    def interval(values: np.ndarray) -> list[float]:
        values = values[np.isfinite(values)]
        if len(values) < iterations * 0.5:
            return [float("nan"), float("nan")]
        return [float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))]

    return {
        "essential_agreement": interval(ea_all),
        "essential_agreement_exact_only": interval(ea_ex),
        "groups": groups,
    }


def stratified(scored: pd.DataFrame, dominant: str, iterations: int, seed: int) -> dict[str, Any]:
    in_dominant = scored["lineage_group"].astype(str).eq(dominant)
    strata: dict[str, Any] = {}
    for label, subset, unit in ((dominant, scored[in_dominant], "isolate_id"),
                                ("other_lineages", scored[~in_dominant], BOOTSTRAP_UNIT)):
        if subset.empty:
            strata[label] = {"n": 0}
            continue
        result = summarize(subset)
        if subset[unit].nunique() >= 2:
            result["bootstrap_95ci"] = lineage_bootstrap(subset, iterations, seed, unit=unit)
            result["bootstrap_unit"] = "isolate (single lineage)" if unit == "isolate_id" else unit
        strata[label] = result
    return strata


def evaluate(
    frame: pd.DataFrame, *, iterations: int, seed: int, locked_external: bool = False,
    dominant_lineage: str | None = None,
) -> dict[str, Any]:
    splits = set(frame["evaluation_split"].astype(str))
    if "external" in splits and not locked_external:
        raise MicEvaluationError(
            "External rows present; the external set is scored once, after models and "
            "thresholds are frozen, with --locked-external."
        )
    if locked_external and splits != {"external"}:
        raise MicEvaluationError("A locked external evaluation must contain external rows only")
    scored = score_rows(frame)
    per_drug: dict[str, Any] = {}
    for antibiotic, group in scored.groupby("antibiotic", sort=True):
        result = summarize(group)
        result["bootstrap_95ci"] = lineage_bootstrap(group, iterations, seed)
        result["by_reference_class"] = by_reference_class(group, iterations, seed)
        if dominant_lineage:
            result["strata"] = stratified(group, dominant_lineage, iterations, seed)
        per_drug[str(antibiotic)] = result
    return {"per_antibiotic": per_drug, "scored": scored}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--study-config", type=Path, default=Path("config/study.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scored-rows", type=Path, help="Optional per-row scoring table")
    parser.add_argument("--locked-external", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = json.loads(args.study_config.read_text(encoding="utf-8"))
    evaluation = config["evaluation"]
    if evaluation.get("bootstrap_unit") != BOOTSTRAP_UNIT:
        raise SystemExit(f"config evaluation.bootstrap_unit must be {BOOTSTRAP_UNIT}")
    frame = pd.read_csv(args.predictions, dtype={"isolate_id": str, "lineage_group": str})
    try:
        result = evaluate(frame, iterations=int(evaluation["cluster_bootstrap_iterations"]),
                          seed=int(evaluation["seed"]), locked_external=args.locked_external,
                          dominant_lineage=evaluation.get("dominant_lineage"))
    except MicEvaluationError as exc:
        raise SystemExit(f"MIC evaluation failed: {exc}") from exc
    if args.scored_rows:
        result["scored"].to_csv(args.scored_rows, index=False)
    payload = {
        "schema_version": "1.0.0",
        "operation": "evaluate_mic_predictions",
        "generated_at_utc": utc_now(),
        "specification": "docs/MIC_EVALUATION_SPEC.md",
        "inputs": {"predictions": {"path": str(args.predictions), "sha256": sha256_file(args.predictions)},
                   "study_config": {"path": str(args.study_config), "sha256": sha256_file(args.study_config)}},
        "locked_external": args.locked_external,
        "bootstrap": {"unit": BOOTSTRAP_UNIT, "iterations": int(evaluation["cluster_bootstrap_iterations"]),
                      "seed": int(evaluation["seed"])},
        "per_antibiotic": result["per_antibiotic"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
