"""Lock thresholds on validation predictions and evaluate an untouched external set."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    matthews_corrcoef,
    roc_auc_score,
)


REQUIRED_COLUMNS = {
    "isolate_id",
    "antibiotic",
    "model",
    "partition",
    "y_true",
    "y_score",
    "genomic_cluster",
    "lineage_group",
}

# Amendment 005: resample by lineage. Single-linkage ANI components chain --
# one external component holds 1,314 of 3,159 isolates -- so resampling them
# would make every interval degenerate. ANI components remain for leakage
# control only.
BOOTSTRAP_UNIT = "lineage_group"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def confusion_metrics(y_true: np.ndarray, y_score: np.ndarray, threshold: float) -> dict[str, float]:
    y_pred = (y_score >= threshold).astype(int)
    tn = int(np.sum((y_true == 0) & (y_pred == 0)))
    fp = int(np.sum((y_true == 0) & (y_pred == 1)))
    fn = int(np.sum((y_true == 1) & (y_pred == 0)))
    tp = int(np.sum((y_true == 1) & (y_pred == 1)))

    def ratio(numerator: float, denominator: float) -> float:
        return float(numerator / denominator) if denominator else float("nan")

    metrics = {
        "n": float(len(y_true)),
        "n_susceptible": float(np.sum(y_true == 0)),
        "n_resistant": float(np.sum(y_true == 1)),
        "tn": float(tn),
        "fp": float(fp),
        "fn": float(fn),
        "tp": float(tp),
        "sensitivity": ratio(tp, tp + fn),
        "specificity": ratio(tn, tn + fp),
        "false_susceptible_rate": ratio(fn, tp + fn),
        "false_resistant_rate": ratio(fp, tn + fp),
        "ppv": ratio(tp, tp + fp),
        "npv": ratio(tn, tn + fn),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "mcc": float(matthews_corrcoef(y_true, y_pred)),
        "brier": float(brier_score_loss(y_true, y_score)),
    }
    if len(np.unique(y_true)) == 2:
        metrics["auroc"] = float(roc_auc_score(y_true, y_score))
        metrics["auprc"] = float(average_precision_score(y_true, y_score))
    else:
        metrics["auroc"] = float("nan")
        metrics["auprc"] = float("nan")
    return metrics


def select_threshold(y_true: np.ndarray, y_score: np.ndarray, min_sensitivity: float) -> tuple[float, str]:
    if len(np.unique(y_true)) < 2:
        raise ValueError("Validation data must contain both susceptible and resistant isolates.")
    candidates = np.unique(np.r_[0.0, y_score, 1.0])
    rows = []
    for threshold in candidates:
        metrics = confusion_metrics(y_true, y_score, float(threshold))
        rows.append((float(threshold), metrics["sensitivity"], metrics["specificity"]))
    feasible = [row for row in rows if row[1] >= min_sensitivity]
    if feasible:
        threshold, _, _ = max(feasible, key=lambda row: (row[2], row[1], row[0]))
        return threshold, f"max_specificity_subject_to_sensitivity>={min_sensitivity:.3f}"
    threshold, _, _ = max(rows, key=lambda row: (row[1], row[2], row[0]))
    return threshold, "fallback_max_sensitivity_then_specificity"


def cluster_bootstrap(
    frame: pd.DataFrame,
    threshold: float,
    iterations: int,
    seed: int,
    unit: str = BOOTSTRAP_UNIT,
) -> dict[str, tuple[float, float]]:
    groups = sorted(frame[unit].astype(str).unique())
    if len(groups) < 2:
        raise ValueError(f"Cluster bootstrap requires at least two {unit} values.")
    grouped = {group: frame.loc[frame[unit].astype(str).eq(group)] for group in groups}
    rng = np.random.default_rng(seed)
    samples: dict[str, list[float]] = {}
    accepted = 0
    for _ in range(iterations):
        chosen = rng.choice(groups, size=len(groups), replace=True)
        sample = pd.concat([grouped[group] for group in chosen], ignore_index=True)
        y_true = sample["y_true"].to_numpy(dtype=int)
        if len(np.unique(y_true)) < 2:
            continue
        metrics = confusion_metrics(y_true, sample["y_score"].to_numpy(dtype=float), threshold)
        for name, value in metrics.items():
            if np.isfinite(value):
                samples.setdefault(name, []).append(float(value))
        accepted += 1
    if accepted < max(100, int(iterations * 0.5)):
        raise ValueError(f"Only {accepted}/{iterations} bootstrap replicates contained both classes.")
    return {
        name: (float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975)))
        for name, values in samples.items()
    }


def evaluate(
    predictions: pd.DataFrame,
    min_sensitivity: float = 0.90,
    bootstrap_iterations: int = 1000,
    seed: int = 20260815,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing = sorted(REQUIRED_COLUMNS.difference(predictions.columns))
    if missing:
        raise ValueError(f"Missing prediction columns: {missing}")
    if predictions.duplicated(["isolate_id", "antibiotic", "model", "partition"]).any():
        raise ValueError("Predictions contain duplicate isolate/drug/model/partition rows.")
    if not set(predictions["partition"].astype(str)).issubset({"validation", "external"}):
        raise ValueError("partition must contain only validation or external.")
    scores = pd.to_numeric(predictions["y_score"], errors="coerce")
    labels = pd.to_numeric(predictions["y_true"], errors="coerce")
    if scores.isna().any() or (~scores.between(0, 1)).any():
        raise ValueError("y_score must be numeric probabilities in [0, 1].")
    if labels.isna().any() or not labels.isin([0, 1]).all():
        raise ValueError("y_true must contain binary labels 0/1.")
    predictions = predictions.assign(y_score=scores.astype(float), y_true=labels.astype(int))

    for (antibiotic, partition), comparison in predictions.groupby(
        ["antibiotic", "partition"], sort=False
    ):
        expected = None
        for model, model_rows in comparison.groupby("model", sort=False):
            signature = set(
                zip(
                    model_rows["isolate_id"].astype(str),
                    model_rows["y_true"].astype(int),
                    model_rows["genomic_cluster"].astype(str),
                    model_rows["lineage_group"].astype(str),
                )
            )
            if expected is None:
                expected = signature
            elif signature != expected:
                raise ValueError(
                    f"Models do not use identical isolates/labels/clusters for {antibiotic}/{partition}: {model}"
                )

    threshold_rows = []
    metric_rows = []
    for (antibiotic, model), task in predictions.groupby(["antibiotic", "model"], sort=True):
        validation = task.loc[task["partition"].eq("validation")]
        external = task.loc[task["partition"].eq("external")]
        if validation.empty or external.empty:
            raise ValueError(f"{antibiotic}/{model} requires both validation and external predictions.")
        threshold, objective = select_threshold(
            validation["y_true"].to_numpy(), validation["y_score"].to_numpy(), min_sensitivity
        )
        threshold_rows.append(
            {
                "antibiotic": antibiotic,
                "model": model,
                "threshold": threshold,
                "selection_partition": "validation",
                "objective": objective,
                "n_validation": len(validation),
            }
        )
        estimates = confusion_metrics(
            external["y_true"].to_numpy(), external["y_score"].to_numpy(), threshold
        )
        intervals = cluster_bootstrap(external, threshold, bootstrap_iterations, seed)
        for metric, estimate in estimates.items():
            lower, upper = intervals.get(metric, (float("nan"), float("nan")))
            metric_rows.append(
                {
                    "antibiotic": antibiotic,
                    "model": model,
                    "partition": "external",
                    "metric": metric,
                    "estimate": estimate,
                    "ci_95_lower": lower,
                    "ci_95_upper": upper,
                    "threshold": threshold,
                    "bootstrap_unit": BOOTSTRAP_UNIT,
                    "bootstrap_iterations": bootstrap_iterations,
                }
            )
    return pd.DataFrame(threshold_rows), pd.DataFrame(metric_rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--thresholds", type=Path, default=Path("results/locked_thresholds.csv"))
    parser.add_argument("--metrics", type=Path, default=Path("results/external_metrics.csv"))
    parser.add_argument("--manifest", type=Path, default=Path("results/evaluation_manifest.json"))
    parser.add_argument("--min-sensitivity", type=float, default=0.90)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260815)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    predictions = pd.read_csv(args.predictions)
    thresholds, metrics = evaluate(
        predictions,
        min_sensitivity=args.min_sensitivity,
        bootstrap_iterations=args.bootstrap_iterations,
        seed=args.seed,
    )
    args.thresholds.parent.mkdir(parents=True, exist_ok=True)
    args.metrics.parent.mkdir(parents=True, exist_ok=True)
    thresholds.to_csv(args.thresholds, index=False)
    metrics.to_csv(args.metrics, index=False)
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "predictions_sha256": sha256_file(args.predictions),
        "thresholds_sha256": sha256_file(args.thresholds),
        "metrics_sha256": sha256_file(args.metrics),
        "threshold_source": "validation_only",
        "evaluation_partition": "external_only",
        "bootstrap_unit": BOOTSTRAP_UNIT,
        "bootstrap_iterations": args.bootstrap_iterations,
        "seed": args.seed,
        "min_sensitivity": args.min_sensitivity,
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Evaluated {len(thresholds)} drug/model tasks on the locked external partition.")


if __name__ == "__main__":
    main()
