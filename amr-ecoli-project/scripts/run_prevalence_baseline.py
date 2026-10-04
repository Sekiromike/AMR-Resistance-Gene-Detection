"""Generate out-of-fold and external predictions for a prevalence-only baseline."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def make_predictions(
    cohort: pd.DataFrame,
    splits: pd.DataFrame,
    intermediate_policy: str = "error",
) -> pd.DataFrame:
    required_cohort = {"isolate_id", "antibiotic", "ast_category", "evaluation_split", "genomic_cluster", "lineage_group"}
    required_splits = {"isolate_id", "evaluation_split", "genomic_cluster", "lineage_group", "cv_fold"}
    missing = sorted(required_cohort.difference(cohort.columns) | required_splits.difference(splits.columns))
    if missing:
        raise ValueError(f"Missing prevalence-baseline columns: {missing}")
    data = cohort.merge(
        splits[["isolate_id", "evaluation_split", "genomic_cluster", "lineage_group", "cv_fold"]],
        on=["isolate_id", "evaluation_split", "genomic_cluster", "lineage_group"],
        how="left",
        validate="many_to_one",
    )
    observed_categories = set(data["ast_category"].astype(str))
    invalid_categories = observed_categories.difference({"S", "I", "R"})
    if invalid_categories:
        raise ValueError(f"Non-canonical AST categories: {sorted(invalid_categories)}")
    if "I" in observed_categories and intermediate_policy == "error":
        raise ValueError("Cohort contains EUCAST I rows; select an explicit intermediate policy.")
    if intermediate_policy not in {"error", "exclude_with_flow_report"}:
        raise ValueError(f"Unsupported intermediate policy: {intermediate_policy}")
    data = data.loc[data["ast_category"].isin(["S", "R"])].copy()
    data["y_true"] = data["ast_category"].map({"S": 0, "R": 1}).astype(int)
    rows = []
    for antibiotic, task in data.groupby("antibiotic", sort=True):
        development = task.loc[task["evaluation_split"].eq("development")]
        external = task.loc[task["evaluation_split"].eq("external")]
        if development.empty or external.empty:
            raise ValueError(f"{antibiotic} requires both development and external rows.")
        if development["cv_fold"].isna().any():
            raise ValueError(f"{antibiotic} has development isolates without a CV fold.")
        for fold in sorted(development["cv_fold"].astype(int).unique()):
            held_out = development.loc[development["cv_fold"].astype(int).eq(fold)]
            training = development.loc[~development["cv_fold"].astype(int).eq(fold)]
            if training.empty:
                raise ValueError(f"{antibiotic} fold {fold} has no training isolates.")
            prevalence = float(training["y_true"].mean())
            for row in held_out.itertuples(index=False):
                rows.append(
                    {
                        "isolate_id": row.isolate_id,
                        "antibiotic": antibiotic,
                        "model": "prevalence-only",
                        "partition": "validation",
                        "y_true": row.y_true,
                        "y_score": prevalence,
                        "genomic_cluster": row.genomic_cluster,
                        "lineage_group": row.lineage_group,
                    }
                )
        development_prevalence = float(development["y_true"].mean())
        for row in external.itertuples(index=False):
            rows.append(
                {
                    "isolate_id": row.isolate_id,
                    "antibiotic": antibiotic,
                    "model": "prevalence-only",
                    "partition": "external",
                    "y_true": row.y_true,
                    "y_score": development_prevalence,
                    "genomic_cluster": row.genomic_cluster,
                        "lineage_group": row.lineage_group,
                }
            )
    predictions = pd.DataFrame(rows)
    if predictions.duplicated(["isolate_id", "antibiotic", "partition"]).any():
        raise ValueError("Baseline prediction construction produced duplicate rows.")
    return predictions.sort_values(["antibiotic", "partition", "isolate_id"]).reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("results/predictions/prevalence.csv"))
    parser.add_argument("--manifest", type=Path, default=Path("results/predictions/prevalence_manifest.json"))
    parser.add_argument(
        "--intermediate-policy",
        choices=["error", "exclude_with_flow_report"],
        default="error",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cohort = pd.read_csv(args.cohort)
    predictions = make_predictions(
        cohort,
        pd.read_csv(args.splits),
        intermediate_policy=args.intermediate_policy,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(args.output, index=False)
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "model": "prevalence-only",
        "cohort_sha256": sha256_file(args.cohort),
        "splits_sha256": sha256_file(args.splits),
        "predictions_sha256": sha256_file(args.output),
        "validation_predictions": "out-of-fold development prevalence",
        "external_predictions": "full-development prevalence",
        "n_predictions": len(predictions),
        "intermediate_policy": args.intermediate_policy,
        "n_intermediate_rows_excluded": int(cohort["ast_category"].eq("I").sum()),
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {len(predictions)} prevalence-only predictions to {args.output}")


if __name__ == "__main__":
    main()
