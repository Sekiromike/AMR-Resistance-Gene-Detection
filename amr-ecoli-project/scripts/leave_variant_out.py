"""Leave-one-variant-out experiment (docs/LEAVE_VARIANT_OUT_PLAN.md).

For a target allele A, every human-clinical development isolate with a
ceftriaxone reference that carries A is held out; models are fitted on the
remaining isolates (settings chosen on the frozen folds restricted to them,
then refitted) and predict the held-out carriers. Development data only.

Subcommands:
    run      one target allele, all learner x feature-set combinations
    readout  pre-stated read-out over all targets (plan section "Read-out")
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from scripts import predict_external as pe
    from scripts import run_mic_models as models
    from scripts.evaluate_mic_predictions import score_rows
except ModuleNotFoundError:  # Direct execution
    import predict_external as pe  # type: ignore[no-redef]
    import run_mic_models as models  # type: ignore[no-redef]
    from evaluate_mic_predictions import score_rows  # type: ignore[no-redef]

DRUG = "ceftriaxone"
POPULATION = "human_clinical"
LEARNERS = ("ridge_censored", "xgboost_aft")
FEATURE_SETS = ("allele", "family", "evo2", "amr+evo2")
PRIMARY_TARGETS = ("blaCTX-M-15", "blaCTX-M-27", "blaCTX-M-14", "blaCTX-M-55")
RESISTANT_SIDE_LOG2 = 1.0
NON_INFERIORITY_MARGIN = -0.05


class LeaveOutError(ValueError):
    """Raised when the experiment contract cannot be honoured."""


def split_on_allele(rows: pd.DataFrame, allele_features: pd.DataFrame, target: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    carriers = set(allele_features.loc[allele_features["element_symbol"].eq(target), "isolate_id"])
    test = rows[rows["isolate_id"].isin(carriers)].reset_index(drop=True)
    train = rows[~rows["isolate_id"].isin(carriers)].reset_index(drop=True)
    if test.empty:
        raise LeaveOutError(f"No carriers of {target}")
    return train, test


def run_target(target: str, development: pd.DataFrame, splits: pd.DataFrame, allele_features: pd.DataFrame,
               family_features: pd.DataFrame, embeddings: dict[str, Any] | None, seed: int,
               feature_sets: tuple[str, ...] = FEATURE_SETS) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    rows = pe.training_rows(development, splits, POPULATION)
    rows = rows[rows["antibiotic"].eq(DRUG)].reset_index(drop=True)
    train, test = split_on_allele(rows, allele_features, target)
    if (allele_features["element_symbol"].eq(target) & allele_features["isolate_id"].isin(set(train["isolate_id"]))).any():
        raise LeaveOutError("A training isolate carries the held-out allele")
    frames, records = [], []
    for feature_set in feature_sets:
        table = family_features if feature_set == "family" else allele_features
        builder_set = "amr" if feature_set in ("allele", "family") else feature_set
        for learner_name in LEARNERS:
            learner = models.LEARNERS[learner_name]
            builder = models.FeatureBuilder(builder_set, table, embeddings)
            mu, sd, chosen, report = pe.fit_learner(train, test["isolate_id"].tolist(), table, learner, builder, seed)
            frames.append(test[["isolate_id", "antibiotic", "lineage_group", "measurement_sign", "ast_value"]].assign(
                evaluation_split="development", predicted_log2_mic=mu, predicted_log2_sd=sd,
                model=f"{learner_name}:{feature_set}", target=target))
            records.append({"target": target, "model": f"{learner_name}:{feature_set}", "n_train": int(len(train)),
                            "n_test": int(len(test)), "selected": chosen, "selection_ea": report.get("selection_ea")})
    return pd.concat(frames, ignore_index=True), records


def resistant_side(frame: pd.DataFrame) -> pd.DataFrame:
    """Right-censored rows with a detected flag. Scored per target: an isolate
    carrying two target alleles appears once per target, and (target, isolate)
    is the unit."""
    parts = []
    for _, part in frame.groupby("target", sort=True):
        scored = score_rows(part.drop(columns=["predicted_log2_sd"]))
        high = scored[np.isposinf(scored["ref_hi"])]
        parts.append(high.assign(detected=(high["predicted_dilution"] > RESISTANT_SIDE_LOG2).astype(float)))
    return pd.concat(parts, ignore_index=True)


def paired_bootstrap(a: pd.DataFrame, b: pd.DataFrame, iterations: int, seed: int) -> dict[str, float]:
    """Difference in detection share (a - b) on identical (target, isolate) rows, resampling lineages."""
    key = ["target", "isolate_id"]
    merged = a[key + ["lineage_group", "detected"]].merge(b[key + ["detected"]], on=key, suffixes=("_a", "_b"),
                                                          validate="one_to_one")
    merged["diff"] = merged["detected_a"] - merged["detected_b"]
    per = merged.groupby("lineage_group")["diff"].agg(["sum", "count"])
    sums, counts = per["sum"].to_numpy(), per["count"].to_numpy()
    draws = np.random.default_rng(seed).integers(0, len(per), size=(iterations, len(per)))
    boot = sums[draws].sum(axis=1) / counts[draws].sum(axis=1)
    return {"difference": float(merged["diff"].mean()), "ci_low": float(np.quantile(boot, 0.025)),
            "ci_high": float(np.quantile(boot, 0.975)), "n": int(len(merged))}


def readout(predictions: pd.DataFrame, iterations: int, seed: int) -> dict[str, Any]:
    out: dict[str, Any] = {"per_target": {}, "pooled_primary": {}, "hypotheses": {}}
    detected = {model: resistant_side(part) for model, part in predictions.groupby("model")}
    for target in sorted(predictions["target"].unique()):
        entry = {}
        for model, part in predictions[predictions["target"].eq(target)].groupby("model"):
            scored = score_rows(part.drop(columns=["predicted_log2_sd"]))
            high = detected[model][detected[model]["target"].eq(target)]
            entry[model] = {"n_test": int(len(part)), "n_right_censored": int(len(high)),
                            "resistant_side_share": float(high["detected"].mean()) if len(high) else None,
                            "essential_agreement": float(scored["essential_agreement"].mean())}
        out["per_target"][target] = entry
    for model, high in detected.items():
        primary = high[high["target"].isin(PRIMARY_TARGETS)]
        out["pooled_primary"][model] = {"n": int(len(primary)), "resistant_side_share": float(primary["detected"].mean())}
    for learner in LEARNERS:
        pick = {fs: detected[f"{learner}:{fs}"] for fs in ("family", "allele", "evo2")}
        pick = {fs: d[d["target"].isin(PRIMARY_TARGETS)] for fs, d in pick.items()}
        h1 = paired_bootstrap(pick["family"], pick["evo2"], iterations, seed)
        h2 = paired_bootstrap(pick["family"], pick["allele"], iterations, seed)
        out["hypotheses"][learner] = {
            "H1_family_minus_evo2": {**h1, "supported": h1["ci_low"] > NON_INFERIORITY_MARGIN},
            "H2_family_minus_allele": {**h2, "supported": h2["ci_low"] > 0.0},
        }
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run")
    run.add_argument("--target", required=True)
    run.add_argument("--cohort", type=Path, required=True)
    run.add_argument("--splits", type=Path, required=True)
    run.add_argument("--features", type=Path, required=True)
    run.add_argument("--family-features", type=Path, required=True)
    run.add_argument("--embeddings", type=Path, required=True)
    run.add_argument("--unit-mapping", type=Path, required=True)
    run.add_argument("--study-config", type=Path, default=Path("config/study.json"))
    run.add_argument("--output", type=Path, required=True)
    ro = sub.add_parser("readout")
    ro.add_argument("--predictions", type=Path, nargs="+", required=True)
    ro.add_argument("--study-config", type=Path, default=Path("config/study.json"))
    ro.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.study_config.read_text(encoding="utf-8"))
    if args.command == "run":
        development, _ = pe.load_cohort_sealed(args.cohort)
        splits = pd.read_csv(args.splits, dtype={"isolate_id": str, "cv_group": str})
        allele = pd.read_csv(args.features, sep="\t", dtype=str, keep_default_na=False)
        family = pd.read_csv(args.family_features, sep="\t", dtype=str, keep_default_na=False)
        from embed_resistance_units import LAYERS  # type: ignore[import-not-found]
        from extract_resistance_loci import PANEL  # type: ignore[import-not-found]
        import embedding_features as emb  # type: ignore[import-not-found]
        embeddings = {"mapping": emb.load_mapping(args.unit_mapping), "panel_loci": list(PANEL),
                      "layers": {layer: emb.load_layer(args.embeddings, layer) for layer in LAYERS}}
        predictions, records = run_target(args.target, development, splits, allele, family, embeddings,
                                          int(config["splitting"]["seed"]))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        predictions.to_csv(args.output, index=False)
        args.output.with_suffix(".json").write_text(json.dumps(records, indent=2, default=str) + "\n", encoding="utf-8")
    else:
        predictions = pd.concat([pd.read_csv(p, dtype={"isolate_id": str, "lineage_group": str}) for p in args.predictions],
                                ignore_index=True)
        evaluation = config["evaluation"]
        result = readout(predictions, int(evaluation["cluster_bootstrap_iterations"]), int(evaluation["seed"]))
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(result["hypotheses"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
