"""Paired model comparison and external precision planning (amendment 007).

1. Paired comparison. For each antibiotic, a model is compared with each
   baseline on identical rows: the difference in essential agreement, with a
   95% interval from resampling whole lineage groups (the same draws for both
   models, so the interval is for the paired difference).

2. Source robustness (MIC_EVALUATION_SPEC addendum B). Each paired
   difference is recomputed with each major source (a BioProject with at least
   MAJOR_SOURCE_SHARE of the drug's rows) excluded in turn.

3. Precision planning. The within-lineage intraclass correlation (ANOVA
   estimator) of the model's essential-agreement errors in development is the
   planning ICC. With the external lineage-group sizes (genotype only; no
   external MIC value or comparator is read) the external 95% half-width at the
   planning EA is 1.96 * sqrt(p (1 - p) DEFF / n), DEFF = 1 + (m_w - 1) ICC,
   m_w the size-weighted mean group size. A drug above the target half-width
   is labelled underpowered.

Run from the repository root:
    python scripts/compare_mic_models.py --model amrfinder=... \\
        --baseline constant=... --baseline neighbor=... \\
        --cohort ... --output ...
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

try:
    from scripts.evaluate_mic_predictions import score_rows
except ModuleNotFoundError:  # Direct execution
    from evaluate_mic_predictions import score_rows  # type: ignore[no-redef]

KEY = ["isolate_id", "antibiotic"]
MAJOR_SOURCE_SHARE = 0.10


class ComparisonError(ValueError):
    """Raised when models cannot be compared on identical rows."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def paired_difference(model: pd.DataFrame, baseline: pd.DataFrame, iterations: int, seed: int) -> dict[str, Any]:
    """EA(model) - EA(baseline) on identical rows with a paired lineage bootstrap."""
    left = score_rows(model)[KEY + ["lineage_group", "essential_agreement"]]
    right = score_rows(baseline)[KEY + ["essential_agreement"]]
    merged = left.merge(right, on=KEY, how="outer", suffixes=("_model", "_baseline"), indicator=True)
    if not merged["_merge"].eq("both").all():
        raise ComparisonError("Model and baseline were not scored on identical rows")
    merged["diff"] = merged["essential_agreement_model"].astype(float) - merged["essential_agreement_baseline"].astype(float)
    per_group = merged.groupby(merged["lineage_group"].astype(str))["diff"].agg(["sum", "count"])
    sums, counts = per_group["sum"].to_numpy(), per_group["count"].to_numpy()
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(sums), size=(iterations, len(sums)))
    boot = sums[draws].sum(axis=1) / counts[draws].sum(axis=1)
    return {
        "n": int(len(merged)),
        "ea_model": float(merged["essential_agreement_model"].mean()),
        "ea_baseline": float(merged["essential_agreement_baseline"].mean()),
        "difference": float(merged["diff"].mean()),
        "difference_95ci": [float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))],
        "lineage_groups": int(len(sums)),
    }


def source_robustness(model: pd.DataFrame, baseline: pd.DataFrame, sources: pd.Series,
                      iterations: int, seed: int) -> dict[str, Any]:
    """Paired difference with each major source excluded in turn."""
    source = model["isolate_id"].map(sources).fillna("unknown")
    shares = source.value_counts(normalize=True)
    out: dict[str, Any] = {}
    for project, share in shares[shares >= MAJOR_SOURCE_SHARE].items():
        keep = set(model.loc[source != project, "isolate_id"])
        if model.loc[model["isolate_id"].isin(keep), "lineage_group"].nunique() < 2:
            # e.g. the single-laboratory external set: nothing left to compare.
            out[str(project)] = {"share_excluded": float(share), "skipped": "too few rows remain"}
            continue
        entry = paired_difference(model[model["isolate_id"].isin(keep)],
                                  baseline[baseline["isolate_id"].isin(keep)], iterations, seed)
        out[str(project)] = {"share_excluded": float(share), **entry}
    return out


def anova_icc(values: np.ndarray, groups: np.ndarray) -> float:
    frame = pd.DataFrame({"y": values.astype(float), "g": groups})
    stats = frame.groupby("g")["y"].agg(["mean", "count"])
    n, k = len(frame), len(stats)
    if k < 2 or n <= k:
        raise ComparisonError("ICC needs at least two groups and more rows than groups")
    grand = frame["y"].mean()
    ssb = float(((stats["mean"] - grand) ** 2 * stats["count"]).sum())
    ssw = float(((frame["y"] - frame.groupby("g")["y"].transform("mean")) ** 2).sum())
    msb, msw = ssb / (k - 1), ssw / (n - k)
    m0 = (n - float((stats["count"] ** 2).sum()) / n) / (k - 1)
    denominator = msb + (m0 - 1) * msw
    return 0.0 if denominator <= 0 else max(0.0, (msb - msw) / denominator)


def weighted_group_size(groups: pd.Series) -> float:
    sizes = groups.value_counts().to_numpy(dtype=float)
    return float((sizes ** 2).sum() / sizes.sum())


def plan_precision(model: pd.DataFrame, external_groups: pd.Series, planning_ea: float,
                   target_half_width: float) -> dict[str, Any]:
    scored = score_rows(model)
    error = (~scored["essential_agreement"]).to_numpy(dtype=float)
    icc = anova_icc(error, scored["lineage_group"].astype(str).to_numpy())
    m_w = weighted_group_size(external_groups)
    deff = 1 + (m_w - 1) * icc
    n = int(len(external_groups))
    half_width = 1.96 * math.sqrt(planning_ea * (1 - planning_ea) * deff / n)
    return {
        "development_error_icc": icc,
        "external_n": n,
        "external_lineage_groups": int(external_groups.nunique()),
        "external_weighted_group_size": m_w,
        "design_effect": deff,
        "effective_n": n / deff,
        "planning_ea": planning_ea,
        "half_width": half_width,
        "target_half_width": target_half_width,
        "status": "adequate" if half_width <= target_half_width else "underpowered",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True, help="name=predictions.csv")
    parser.add_argument("--baseline", action="append", required=True, help="name=predictions.csv; repeatable")
    parser.add_argument("--cohort", type=Path, required=True, help="for external lineage-group sizes only")
    parser.add_argument("--study-config", type=Path, default=Path("config/study.json"))
    parser.add_argument("--planning-ea", type=float, default=0.90)
    parser.add_argument("--target-half-width", type=float, default=0.05)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    evaluation = json.loads(args.study_config.read_text(encoding="utf-8"))["evaluation"]
    iterations, seed = int(evaluation["cluster_bootstrap_iterations"]), int(evaluation["seed"])

    def load(spec: str) -> tuple[str, Path, pd.DataFrame]:
        name, _, path = spec.partition("=")
        return name, Path(path), pd.read_csv(path, dtype={"isolate_id": str, "lineage_group": str})

    model_name, model_path, model = load(args.model)
    baselines = [load(spec) for spec in args.baseline]
    # Genotype columns only: no external MIC value or comparator is read.
    cohort = pd.read_csv(args.cohort, usecols=["isolate_id", "antibiotic", "evaluation_split", "lineage_group",
                                               "bioproject_accession"], dtype=str, keep_default_na=False)
    external = cohort[cohort["evaluation_split"].eq("external")]
    sources = cohort.drop_duplicates("isolate_id").set_index("isolate_id")["bioproject_accession"]

    result: dict[str, Any] = {}
    try:
        for antibiotic in sorted(model["antibiotic"].unique()):
            m = model[model["antibiotic"].eq(antibiotic)]
            entry: dict[str, Any] = {"comparisons": {}}
            entry["source_robustness"] = {}
            for name, _, frame in baselines:
                b = frame[frame["antibiotic"].eq(antibiotic)]
                entry["comparisons"][name] = paired_difference(m, b, iterations, seed)
                entry["source_robustness"][name] = source_robustness(m, b, sources, iterations, seed)
            groups = external.loc[external["antibiotic"].eq(antibiotic), "lineage_group"]
            entry["external_precision"] = plan_precision(m, groups, args.planning_ea, args.target_half_width)
            result[str(antibiotic)] = entry
    except ComparisonError as exc:
        raise SystemExit(f"Model comparison failed: {exc}") from exc
    payload = {
        "schema_version": "1.0.0",
        "operation": "compare_mic_models",
        "generated_at_utc": utc_now(),
        "model": model_name,
        "inputs": {name: {"path": str(path), "sha256": sha256_file(path)}
                   for name, path, _ in [(model_name, model_path, None), *baselines]}
                  | {"cohort": {"path": str(args.cohort), "sha256": sha256_file(args.cohort),
                                "columns_read": ["isolate_id", "antibiotic", "evaluation_split", "lineage_group",
                                                 "bioproject_accession"]}},
        "source_robustness": {"major_source_share": MAJOR_SOURCE_SHARE,
                              "specification": "docs/MIC_EVALUATION_SPEC.md addendum B"},
        "bootstrap": {"unit": "lineage_group", "iterations": iterations, "seed": seed, "paired": True},
        "per_antibiotic": result,
        "scientific_boundary": {"external_mic_read": False},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
