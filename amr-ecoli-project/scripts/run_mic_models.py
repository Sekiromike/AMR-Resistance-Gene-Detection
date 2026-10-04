"""Nested cross-validation for MIC models on AMRFinderPlus determinants.

Learners (each a hyperparameter grid ordered simplest first; inner ties go to
the simpler setting):

- ridge_censored: ridge-penalised interval-censored normal regression
  (scripts/run_censored_regression.py). Round 2 widens the penalty grid
  downward to 0.001 because round 1 (job 64960056) selected its smallest
  value, 0.1, in most folds.
- xgboost_aft: gradient-boosted trees with XGBoost's accelerated-failure-time
  objective, which fits interval-censored targets natively; it can learn
  interactions such as gyrA plus parC mutations that an additive model cannot.
  MIC bounds are passed on the mg/L scale (latent MIC in (2^(lo-1), 2^hi]).

Feature sets (docs/FOUNDATION_MODEL_PROTOCOL.md revision 0.2):
- amr: AMRFinderPlus determinants only;
- amr+evo2: determinants plus Evo 2 resistance-unit features;
- evo2: Evo 2 resistance-unit features only;
- amr+evo2_panel / amr+evo2_acquired: diagnostic ablation adding only the
  chromosomal-panel or only the acquired-unit Evo 2 features, to locate a gain.
Evo 2 features come from scripts/embedding_features.py; their PCA is refitted
inside every inner and outer training set, and the embedding layer is chosen
on inner folds together with the learner's hyperparameters.

Nested CV (docs/EVALUATION_DESIGN_AMENDMENT.md): outer folds from the frozen
split table only report; features are filtered and hyperparameters chosen on
inner grouped folds within each outer training set; the chosen setting is
refitted on the whole outer training set and predicts the held-out fold once.
Development only.

Run from the repository root:
    python scripts/run_mic_models.py --learner xgboost_aft --cohort ... \\
        --splits ... --features ... --population human_clinical --output ... --manifest ...
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
    from scripts import embedding_features as emb
    from scripts import run_censored_regression as ridge
    from scripts.evaluate_mic_predictions import reference_bounds
except ModuleNotFoundError:  # Direct execution
    import embedding_features as emb  # type: ignore[no-redef]
    import run_censored_regression as ridge  # type: ignore[no-redef]
    from evaluate_mic_predictions import reference_bounds  # type: ignore[no-redef]

FEATURE_SETS = ("amr", "amr+evo2", "evo2", "amr+evo2_panel", "amr+evo2_acquired")
EVO2_PARTS = {"amr+evo2": ("panel", "acquired"), "evo2": ("panel", "acquired"),
              "amr+evo2_panel": ("panel",), "amr+evo2_acquired": ("acquired",)}

POPULATIONS = ("human_clinical", "all_sources")
INNER_FOLDS = 4
MIN_FEATURE_COUNT = 5
LN2 = math.log(2.0)


class RidgeCensored:
    name = "ridge_censored"
    # Simplest (strongest penalty) first.
    grid = [{"penalty": p} for p in (1000.0, 100.0, 10.0, 1.0, 0.1, 0.01, 0.001)]

    def __init__(self) -> None:
        self.fits = 0
        self.not_converged: list[str] = []

    def fit(self, X: np.ndarray, lo: np.ndarray, hi: np.ndarray, params: dict[str, Any], seed: int) -> Any:
        theta, converged, message = ridge.fit_checked(X, lo, hi, params["penalty"])
        self.fits += 1
        if not converged:
            self.not_converged.append(message)
        return theta

    def predict(self, model: Any, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        mu, sigma = ridge.predict(model, X)
        return mu, np.full(len(mu), sigma)

    def describe(self, model: Any, vocabulary: list[str]) -> dict[str, Any]:
        effects = sorted(zip(vocabulary, model[2:]), key=lambda item: -abs(item[1]))[:10]
        return {"intercept_log2": float(model[0]), "sigma_log2": float(np.exp(model[1])),
                "top_features": [{"feature": f, "log2_effect": float(c)} for f, c in effects]}

    def take_fit_report(self) -> dict[str, Any]:
        report = {"fits": self.fits, "not_converged": len(self.not_converged),
                  "not_converged_messages": sorted(set(self.not_converged))}
        self.fits, self.not_converged = 0, []
        return report


class XGBoostAFT:
    name = "xgboost_aft"
    # Simplest (shallow, fewer trees) first. AFT scale is in natural-log units.
    grid = [{"max_depth": d, "n_trees": t, "aft_scale": 1.0} for d in (2, 3, 4) for t in (100, 300)]

    def _dmatrix(self, X: np.ndarray, lo: np.ndarray | None = None, hi: np.ndarray | None = None) -> Any:
        import xgboost as xgb
        matrix = xgb.DMatrix(X)
        if lo is not None:
            lower = np.where(np.isneginf(lo), 0.0, np.exp2(lo - 1.0))
            upper = np.where(np.isposinf(hi), np.inf, np.exp2(hi))
            matrix.set_float_info("label_lower_bound", lower)
            matrix.set_float_info("label_upper_bound", upper)
        return matrix

    def fit(self, X: np.ndarray, lo: np.ndarray, hi: np.ndarray, params: dict[str, Any], seed: int) -> Any:
        import xgboost as xgb
        booster = xgb.train(
            {"objective": "survival:aft", "eval_metric": "aft-nloglik",
             "aft_loss_distribution": "normal", "aft_loss_distribution_scale": params["aft_scale"],
             "max_depth": params["max_depth"], "learning_rate": 0.05, "tree_method": "hist",
             "min_child_weight": 5, "subsample": 1.0, "seed": seed, "nthread": 1},
            self._dmatrix(X, lo, hi), num_boost_round=params["n_trees"])
        return {"booster": booster, "scale": params["aft_scale"]}

    def predict(self, model: Any, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        # Prediction is exp(margin) on the mg/L scale; its log2 is the location.
        mic = model["booster"].predict(self._dmatrix(X))
        return np.log2(mic), np.full(len(mic), model["scale"] / LN2)

    def take_fit_report(self) -> dict[str, Any]:
        return {}

    def describe(self, model: Any, vocabulary: list[str]) -> dict[str, Any]:
        gain = model["booster"].get_score(importance_type="total_gain")
        named = sorted(((vocabulary[int(k[1:])], v) for k, v in gain.items()), key=lambda item: -item[1])[:10]
        return {"top_features": [{"feature": f, "total_gain": float(v)} for f, v in named]}


LEARNERS = {learner.name: learner for learner in (RidgeCensored(), XGBoostAFT())}


class FeatureBuilder:
    """Assembles the design matrix for one feature set; Evo 2 parts are fitted in-fold."""

    def __init__(self, feature_set: str, features: pd.DataFrame, embeddings: dict[str, Any] | None):
        if feature_set not in FEATURE_SETS:
            raise ridge.ModelError(f"feature set must be one of {FEATURE_SETS}")
        if feature_set != "amr" and embeddings is None:
            raise ridge.ModelError(f"feature set {feature_set} needs embeddings")
        self.feature_set, self.features, self.embeddings = feature_set, features, embeddings

    def layers(self) -> list[str | None]:
        return [None] if self.feature_set == "amr" else list(self.embeddings["layers"])

    def build(self, fit_ids: list[str], apply_ids: list[list[str]], vocabulary: list[str], layer: str | None,
              seed: int) -> tuple[list[np.ndarray], list[str]]:
        """Fit on fit_ids, return a matrix for each list in apply_ids, plus column names."""
        blocks: list[list[np.ndarray]] = [[] for _ in apply_ids]
        names: list[str] = []
        if self.feature_set != "evo2":
            for i, ids in enumerate(apply_ids):
                blocks[i].append(ridge.design(pd.Series(ids), self.features, vocabulary))
            names += vocabulary
        if self.feature_set in EVO2_PARTS:
            index, vectors = self.embeddings["layers"][layer]
            featurizer = emb.EmbeddingFeaturizer(self.embeddings["mapping"], index, vectors,
                                                 self.embeddings["panel_loci"],
                                                 EVO2_PARTS[self.feature_set]).fit(fit_ids, seed)
            for i, ids in enumerate(apply_ids):
                blocks[i].append(featurizer.transform(ids))
            names += [f"evo2[{j}]" for j in range(featurizer.width())]
        return [np.hstack(b) for b in blocks], names


def run_nested(cohort: pd.DataFrame, splits: pd.DataFrame, features: pd.DataFrame, population: str,
               learner: Any, seed: int, feature_set: str = "amr",
               embeddings: dict[str, Any] | None = None) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    builder = FeatureBuilder(feature_set, features, embeddings)
    if population not in POPULATIONS:
        raise ridge.ModelError(f"population must be one of {POPULATIONS}")
    development = cohort.loc[cohort["evaluation_split"].astype(str).eq("development")].copy()
    if population == "human_clinical":
        development = development.loc[development["intended_use_population"].astype(str).eq("human_clinical")]
    folds = splits.loc[splits["evaluation_split"].astype(str).eq("development"), ["isolate_id", "cv_fold", "cv_group"]]
    development = development.merge(folds, on="isolate_id", how="left", validate="many_to_one")
    if development["cv_fold"].isna().any() or development["cv_group"].isna().any():
        raise ridge.ModelError("Every development isolate needs a cv_fold and cv_group")
    development["cv_fold"] = development["cv_fold"].astype(int)
    bounds = [reference_bounds(str(s), float(v)) for s, v in zip(development["measurement_sign"], development["ast_value"])]
    development["_lo"] = [b[0] for b in bounds]
    development["_hi"] = [b[1] for b in bounds]

    predictions: list[pd.DataFrame] = []
    records: list[dict[str, Any]] = []
    for antibiotic, drug in development.groupby("antibiotic", sort=True):
        for fold in sorted(drug["cv_fold"].unique()):
            train = drug.loc[drug["cv_fold"] != fold].reset_index(drop=True)
            test = drug.loc[drug["cv_fold"] == fold].reset_index(drop=True)
            learner.take_fit_report()  # counts start afresh for this outer fold
            present = features[features["isolate_id"].isin(set(train["isolate_id"]))]
            counts = present.groupby("element_symbol")["isolate_id"].nunique()
            vocabulary = sorted(counts[counts >= MIN_FEATURE_COUNT].index)
            lo, hi = train["_lo"].to_numpy(float), train["_hi"].to_numpy(float)
            train_ids = train["isolate_id"].tolist()

            inner = ridge.inner_folds(train["cv_group"].astype(str), INNER_FOLDS, seed)
            layers = builder.layers()
            hits = np.zeros((len(layers), len(learner.grid)))
            for li, layer in enumerate(layers):
                for k in range(INNER_FOLDS):
                    fit_rows, eval_rows = inner != k, inner == k
                    fit_ids = [i for i, keep in zip(train_ids, fit_rows) if keep]
                    eval_ids = [i for i, keep in zip(train_ids, eval_rows) if keep]
                    (X_fit, X_eval), _ = builder.build(fit_ids, [fit_ids, eval_ids], vocabulary, layer, seed)
                    for gi, params in enumerate(learner.grid):
                        model = learner.fit(X_fit, lo[fit_rows], hi[fit_rows], params, seed)
                        mu, _ = learner.predict(model, X_eval)
                        hits[li, gi] += ridge.essential_agreement(mu, lo[eval_rows], hi[eval_rows]) * eval_rows.sum()
            scores = hits / len(train)
            # Ties: earlier layer, then the simpler setting.
            li_best, gi_best = max(((li, gi) for li in range(len(layers)) for gi in range(len(learner.grid))),
                                   key=lambda t: (scores[t], -t[0], -t[1]))
            best, layer = learner.grid[gi_best], layers[li_best]

            test_ids = test["isolate_id"].tolist()
            (X_train, X_test), names = builder.build(train_ids, [train_ids, test_ids], vocabulary, layer, seed)
            model = learner.fit(X_train, lo, hi, best, seed)
            mu, sd = learner.predict(model, X_test)
            chosen = {**best, "layer": layer}
            predictions.append(test.assign(predicted_log2_mic=mu, predicted_log2_sd=sd,
                                           selected=json.dumps(chosen, sort_keys=True)))
            records.append({
                "antibiotic": str(antibiotic), "cv_fold": int(fold), "n_train": int(len(train)),
                "n_test": int(len(test)), "n_features": len(names), "n_amr_features": len(vocabulary),
                "inner_ea": [{"layer": layers[li], "params": learner.grid[gi], "ea": float(scores[li, gi])}
                             for li in range(len(layers)) for gi in range(len(learner.grid))],
                "selected": chosen, "selected_at_grid_edge": gi_best in (0, len(learner.grid) - 1),
                "optimiser": learner.take_fit_report(),
                **learner.describe(model, names),
            })
    columns = ["isolate_id", "antibiotic", "evaluation_split", "lineage_group", "intended_use_population",
               "cv_fold", "measurement_sign", "ast_value", "predicted_log2_mic", "predicted_log2_sd", "selected"]
    output = pd.concat(predictions, ignore_index=True)[columns]
    return output.sort_values(["antibiotic", "isolate_id"]).reset_index(drop=True), records


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--learner", choices=sorted(LEARNERS), required=True)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--population", choices=POPULATIONS, required=True)
    parser.add_argument("--feature-set", choices=FEATURE_SETS, default="amr")
    parser.add_argument("--embeddings", type=Path, help="embedding shard directory (Evo 2 feature sets)")
    parser.add_argument("--unit-mapping", type=Path, help="unit_mapping.tsv from extract_resistance_loci.py")
    parser.add_argument("--study-config", type=Path, default=Path("config/study.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    learner = LEARNERS[args.learner]
    seed = int(json.loads(args.study_config.read_text(encoding="utf-8"))["splitting"]["seed"])
    cohort = pd.read_csv(args.cohort, dtype={"isolate_id": str, "lineage_group": str}, low_memory=False)
    splits = pd.read_csv(args.splits, dtype={"isolate_id": str, "cv_group": str})
    features = pd.read_csv(args.features, sep="\t", dtype=str, keep_default_na=False)
    embeddings = None
    if args.feature_set != "amr":
        if args.embeddings is None or args.unit_mapping is None:
            raise SystemExit("Evo 2 feature sets need --embeddings and --unit-mapping")
        try:
            from scripts.extract_resistance_loci import PANEL
        except ModuleNotFoundError:
            from extract_resistance_loci import PANEL  # type: ignore[no-redef]
        try:
            from scripts.embed_resistance_units import LAYERS
        except ModuleNotFoundError:
            from embed_resistance_units import LAYERS  # type: ignore[no-redef]
        try:
            embeddings = {"mapping": emb.load_mapping(args.unit_mapping), "panel_loci": list(PANEL),
                          "layers": {layer: emb.load_layer(args.embeddings, layer) for layer in LAYERS}}
        except emb.EmbeddingFeatureError as exc:
            raise SystemExit(f"Embedding features failed: {exc}") from exc
    try:
        output, records = run_nested(cohort, splits, features, args.population, learner, seed,
                                     args.feature_set, embeddings)
    except ridge.ModelError as exc:
        raise SystemExit(f"{learner.name} failed: {exc}") from exc
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    versions: dict[str, str] = {"numpy": np.__version__, "pandas": pd.__version__}
    if learner.name == "xgboost_aft":
        import xgboost
        versions["xgboost"] = xgboost.__version__
    payload = {
        "schema_version": "1.0.0",
        "operation": "run_mic_models",
        "learner": learner.name,
        "feature_set": args.feature_set,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "population": args.population,
        "nested_cv": {"grid": learner.grid, "inner_folds": INNER_FOLDS, "min_feature_count": MIN_FEATURE_COUNT,
                      "selection_metric": "inner essential agreement", "tie_break": "simpler setting", "seed": seed},
        "versions": versions,
        "inputs": {name: {"path": str(path), "sha256": sha256_file(path)} for name, path in
                   (("cohort", args.cohort), ("splits", args.splits), ("features", args.features),
                    ("study_config", args.study_config))},
        "folds": records,
        "outputs": {"predictions": {"path": str(args.output), "sha256": sha256_file(args.output),
                                    "rows": int(len(output))}},
        "scientific_boundary": {"external_rows_read": False},
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
