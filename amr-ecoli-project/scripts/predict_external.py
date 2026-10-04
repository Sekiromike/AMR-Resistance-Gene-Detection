"""Frozen-model predictions for the locked external evaluation (seal step 1).

docs/EXTERNAL_EVALUATION_PLAN.md. Each model is fitted on development only:
hyperparameters (and, for Evo 2 feature sets, the layer) are chosen by grouped
cross-validation over the five frozen development folds, then the model is
refitted on all development rows of the training population and predicts the
external rows.

The seal is enforced here, not assumed. The cohort is streamed and every
external row is reduced to genotype columns as it is read, so no external MIC
value or comparator is parsed, retained or used. Predictions carry no
reference columns; seal step 2 checksums and records them before
scripts/attach_external_references.py joins references for evaluation.

Models: constant, neighbor, and the run_mic_models learners on any feature set.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from scripts import run_censored_regression as ridge
    from scripts import run_mic_models as models
    from scripts.evaluate_mic_predictions import reference_bounds
    from scripts.run_mic_baseline import best_constant_dilution
    from scripts.run_neighbor_baseline import boundary_prediction
except ModuleNotFoundError:  # Direct execution
    import run_censored_regression as ridge  # type: ignore[no-redef]
    import run_mic_models as models  # type: ignore[no-redef]
    from evaluate_mic_predictions import reference_bounds  # type: ignore[no-redef]
    from run_mic_baseline import best_constant_dilution  # type: ignore[no-redef]
    from run_neighbor_baseline import boundary_prediction  # type: ignore[no-redef]

GENOTYPE_COLUMNS = ("isolate_id", "antibiotic", "evaluation_split", "lineage_group", "intended_use_population")
OUTPUT_COLUMNS = ["isolate_id", "antibiotic", "evaluation_split", "lineage_group", "intended_use_population",
                  "predicted_log2_mic", "predicted_log2_sd", "model", "training_population", "selected"]
NEIGHBOR_TOP_K = 200


class ExternalPredictionError(ValueError):
    """Raised when the frozen-model contract cannot be honoured."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_cohort_sealed(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Development rows in full; external rows reduced to genotype columns on read."""
    development, external = [], []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            split = row["evaluation_split"]
            if split == "development":
                development.append(row)
            elif split == "external":
                external.append({column: row[column] for column in GENOTYPE_COLUMNS})
            else:
                raise ExternalPredictionError(f"Unknown evaluation_split {split!r}")
    return pd.DataFrame(development), pd.DataFrame(external, columns=list(GENOTYPE_COLUMNS))


def training_rows(development: pd.DataFrame, splits: pd.DataFrame, population: str) -> pd.DataFrame:
    rows = development
    if population == "human_clinical":
        rows = rows[rows["intended_use_population"].eq("human_clinical")]
    folds = splits.loc[splits["evaluation_split"].eq("development"), ["isolate_id", "cv_fold", "cv_group"]]
    rows = rows.merge(folds, on="isolate_id", how="left", validate="many_to_one")
    if rows["cv_fold"].isna().any():
        raise ExternalPredictionError("Every development isolate needs a frozen cv_fold")
    rows = rows.assign(cv_fold=rows["cv_fold"].astype(int))
    bounds = [reference_bounds(str(s), float(v)) for s, v in zip(rows["measurement_sign"], rows["ast_value"])]
    return rows.assign(_lo=[b[0] for b in bounds], _hi=[b[1] for b in bounds]).reset_index(drop=True)


def fit_learner(train: pd.DataFrame, external_ids: list[str], features: pd.DataFrame, learner: Any,
                builder: Any, seed: int) -> tuple[np.ndarray, np.ndarray, dict[str, Any], dict[str, Any]]:
    """Choose settings on the frozen folds, refit on all rows, predict external_ids."""
    present = features[features["isolate_id"].isin(set(train["isolate_id"]))]
    counts = present.groupby("element_symbol")["isolate_id"].nunique()
    vocabulary = sorted(counts[counts >= models.MIN_FEATURE_COUNT].index)
    lo, hi = train["_lo"].to_numpy(float), train["_hi"].to_numpy(float)
    ids = train["isolate_id"].tolist()
    folds = train["cv_fold"].to_numpy()
    layers = builder.layers()
    hits = np.zeros((len(layers), len(learner.grid)))
    learner.take_fit_report()
    for li, layer in enumerate(layers):
        for k in sorted(set(folds)):
            fit_rows, eval_rows = folds != k, folds == k
            fit_ids = [i for i, keep in zip(ids, fit_rows) if keep]
            eval_ids = [i for i, keep in zip(ids, eval_rows) if keep]
            (X_fit, X_eval), _ = builder.build(fit_ids, [fit_ids, eval_ids], vocabulary, layer, seed)
            for gi, params in enumerate(learner.grid):
                model = learner.fit(X_fit, lo[fit_rows], hi[fit_rows], params, seed)
                mu, _ = learner.predict(model, X_eval)
                hits[li, gi] += ridge.essential_agreement(mu, lo[eval_rows], hi[eval_rows]) * eval_rows.sum()
    scores = hits / len(train)
    li_best, gi_best = max(((li, gi) for li in range(len(layers)) for gi in range(len(learner.grid))),
                           key=lambda t: (scores[t], -t[0], -t[1]))
    chosen = {**learner.grid[gi_best], "layer": layers[li_best]}
    (X_train, X_ext), names = builder.build(ids, [ids, external_ids], vocabulary, layers[li_best], seed)
    model = learner.fit(X_train, lo, hi, learner.grid[gi_best], seed)
    mu, sd = learner.predict(model, X_ext)
    report = {"selection_ea": float(scores[li_best, gi_best]), "n_features": len(names),
              "optimiser": learner.take_fit_report(), **learner.describe(model, names)}
    return mu, sd, chosen, report


def neighbor_predictions(comparisons: Path, external: pd.DataFrame, train: pd.DataFrame) -> pd.Series:
    """Boundary of the highest-ANI training relative with that drug's MIC, per external row."""
    reference = {(r.isolate_id, r.antibiotic): (r.measurement_sign, float(r.ast_value)) for r in train.itertuples()}
    pool = set(train["isolate_id"])
    wanted = set(external["isolate_id"])
    heaps: dict[str, list] = {}
    with comparisons.open("r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[:3] != ["development_sequence_id", "external_sequence_id", "ani_percent"]:
            raise ExternalPredictionError(f"Unexpected comparison header {header}")
        for line in handle:
            dev, ext, ani, af_dev, af_ext = line.rstrip("\n").split("\t")[:5]
            dev, ext = dev.removeprefix("development::"), ext.removeprefix("external::")
            if ext not in wanted or dev not in pool:
                continue
            item = (float(ani), min(float(af_dev), float(af_ext)), "".join(chr(0x10FFFF - ord(c)) for c in dev))
            heap = heaps.setdefault(ext, [])
            if len(heap) < NEIGHBOR_TOP_K:
                heapq.heappush(heap, item)
            elif item > heap[0]:
                heapq.heapreplace(heap, item)
    ranked = {ext: ["".join(chr(0x10FFFF - ord(c)) for c in r) for _, _, r in sorted(heap, reverse=True)]
              for ext, heap in heaps.items()}
    out = []
    for row in external.itertuples():
        relative = next((d for d in ranked.get(row.isolate_id, []) if (d, row.antibiotic) in reference), None)
        if relative is None:
            raise ExternalPredictionError(f"No development relative with {row.antibiotic} for {row.isolate_id}")
        out.append(boundary_prediction(*reference[(relative, row.antibiotic)]))
    return pd.Series(out, index=external.index)


def predict(model_name: str, development: pd.DataFrame, external: pd.DataFrame, splits: pd.DataFrame,
            features: pd.DataFrame, population: str, seed: int, embeddings: dict[str, Any] | None = None,
            comparisons: Path | None = None) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    train_all = training_rows(development, splits, population)
    frames, records = [], []
    for antibiotic in sorted(external["antibiotic"].unique()):
        train = train_all[train_all["antibiotic"].eq(antibiotic)].reset_index(drop=True)
        ext = external[external["antibiotic"].eq(antibiotic)].reset_index(drop=True)
        sd = np.full(len(ext), np.nan)
        if model_name == "constant":
            dilution, ea, _ = best_constant_dilution(train["_lo"].to_numpy(float), train["_hi"].to_numpy(float))
            mu, chosen, report = np.full(len(ext), float(dilution)), {"dilution_log2": dilution}, {"training_ea": ea}
        elif model_name == "neighbor":
            if comparisons is None:
                raise ExternalPredictionError("neighbor needs --comparisons")
            mu, chosen, report = neighbor_predictions(comparisons, ext, train).to_numpy(float), {}, {}
        else:
            learner_name, _, feature_set = model_name.partition(":")
            learner = models.LEARNERS[learner_name]
            builder = models.FeatureBuilder(feature_set or "amr", features, embeddings)
            mu, sd, chosen, report = fit_learner(train, ext["isolate_id"].tolist(), features, learner, builder, seed)
        frames.append(ext.assign(predicted_log2_mic=mu, predicted_log2_sd=sd, model=model_name,
                                 training_population=population, selected=json.dumps(chosen, sort_keys=True)))
        records.append({"antibiotic": antibiotic, "n_train": int(len(train)), "n_external": int(len(ext)),
                        "selected": chosen, **report})
    return pd.concat(frames, ignore_index=True)[OUTPUT_COLUMNS], records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True, help="constant | neighbor | <learner>:<feature set>")
    parser.add_argument("--population", choices=models.POPULATIONS, required=True)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--embeddings", type=Path)
    parser.add_argument("--unit-mapping", type=Path)
    parser.add_argument("--comparisons", type=Path)
    parser.add_argument("--study-config", type=Path, default=Path("config/study.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    seed = int(json.loads(args.study_config.read_text(encoding="utf-8"))["splitting"]["seed"])
    development, external = load_cohort_sealed(args.cohort)
    splits = pd.read_csv(args.splits, dtype={"isolate_id": str, "cv_group": str})
    features = pd.read_csv(args.features, sep="\t", dtype=str, keep_default_na=False)
    embeddings = None
    if "evo2" in args.model:
        try:
            from scripts.embed_resistance_units import LAYERS
            from scripts.extract_resistance_loci import PANEL
            from scripts import embedding_features as emb
        except ModuleNotFoundError:
            from embed_resistance_units import LAYERS  # type: ignore[no-redef]
            from extract_resistance_loci import PANEL  # type: ignore[no-redef]
            import embedding_features as emb  # type: ignore[no-redef]
        embeddings = {"mapping": emb.load_mapping(args.unit_mapping), "panel_loci": list(PANEL),
                      "layers": {layer: emb.load_layer(args.embeddings, layer) for layer in LAYERS}}
    try:
        output, records = predict(args.model, development, external, splits, features, args.population, seed,
                                  embeddings, args.comparisons)
    except (ExternalPredictionError, ridge.ModelError) as exc:
        raise SystemExit(f"External prediction failed: {exc}") from exc
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    payload = {
        "schema_version": "1.0.0",
        "operation": "predict_external",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "model": args.model,
        "training_population": args.population,
        "plan": "docs/EXTERNAL_EVALUATION_PLAN.md",
        "inputs": {name: {"path": str(path), "sha256": sha256_file(path)} for name, path in
                   (("cohort", args.cohort), ("splits", args.splits), ("features", args.features))},
        "per_antibiotic": records,
        "outputs": {"predictions": {"path": str(args.output), "sha256": sha256_file(args.output),
                                    "rows": int(len(output))}},
        "scientific_boundary": {"external_reference_columns_read": False,
                                "external_columns_retained": list(GENOTYPE_COLUMNS)},
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
