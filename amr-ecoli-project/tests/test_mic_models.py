from __future__ import annotations

import json
import unittest

import numpy as np
import pandas as pd

from scripts import evaluate_mic_predictions as EVALUATE
from scripts import run_mic_models as MODULE
from tests_support_mic import small_cohort


class LearnerTests(unittest.TestCase):
    def test_xgboost_aft_learns_an_interaction_the_additive_model_cannot(self) -> None:
        # MIC is high only when BOTH mutations are present (epistasis).
        rng = np.random.default_rng(3)
        X = rng.binomial(1, 0.5, size=(2000, 2)).astype(float)
        latent = -3.0 + 5.0 * (X[:, 0] * X[:, 1]) + rng.normal(0, 0.4, 2000)
        d = np.ceil(latent)
        lo = np.where(d <= -2, -np.inf, np.where(d > 3, 4.0, d))
        hi = np.where(d <= -2, -2.0, np.where(d > 3, np.inf, d))
        xgb = MODULE.LEARNERS["xgboost_aft"]
        model = xgb.fit(X, lo, hi, {"max_depth": 2, "n_trees": 300, "aft_scale": 1.0}, seed=1)
        grid = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
        mu, sd = xgb.predict(model, grid)
        self.assertLess(mu[1] - mu[0], 1.0)          # one mutation alone: little effect
        self.assertLess(mu[2] - mu[0], 1.0)
        self.assertGreater(mu[3] - mu[0], 3.0)       # both together: large effect
        self.assertTrue(np.allclose(sd, 1.0 / np.log(2.0)))

    def test_grids_are_ordered_simplest_first(self) -> None:
        penalties = [p["penalty"] for p in MODULE.LEARNERS["ridge_censored"].grid]
        self.assertEqual(penalties, sorted(penalties, reverse=True))
        self.assertIn(0.001, penalties)
        depths = [p["max_depth"] for p in MODULE.LEARNERS["xgboost_aft"].grid]
        self.assertEqual(depths, sorted(depths))


class NestedTests(unittest.TestCase):
    def test_both_learners_run_nested_cv_without_leakage(self) -> None:
        cohort, splits, features = small_cohort()
        for name in ("ridge_censored", "xgboost_aft"):
            predictions, records = MODULE.run_nested(cohort, splits, features, "human_clinical",
                                                     MODULE.LEARNERS[name], seed=7)
            self.assertEqual(predictions["isolate_id"].nunique(), 600, name)
            fold0 = next(r for r in records if r["cv_fold"] == 0)
            self.assertEqual(fold0["n_features"], 3, name)   # fold-0-only feature excluded
            self.assertNotIn("fold0_only", [f["feature"] for f in fold0["top_features"]], name)
            result = EVALUATE.evaluate(predictions, iterations=200, seed=1)["per_antibiotic"]["ciprofloxacin"]
            self.assertGreater(result["essential_agreement"], 0.9, name)


if __name__ == "__main__":
    unittest.main()


def synthetic_embeddings(features: pd.DataFrame, isolates: list[str], dim: int = 12, seed: int = 5):
    """gyrA embedding encodes the gyrA mutation; acquired qnr has its own vector."""
    rng = np.random.default_rng(seed)
    carriers = {name: set(features.loc[features["element_symbol"] == name, "isolate_id"])
                for name in ("gyrA_S83L", "qnrS1")}
    hashes = {"gyrA_wt": 0, "gyrA_mut": 1, "parC_wt": 2, "qnr": 3}
    base = rng.normal(size=(4, dim)).astype(np.float32)
    layers = {}
    for k, layer in enumerate(("blocks.6.mlp.l3", "blocks.12.mlp.l3")):
        vectors = base.copy() if k == 0 else rng.normal(size=(4, dim)).astype(np.float32) * 0.01
        layers[layer] = ({h: i for h, i in hashes.items()}, vectors)
    mapping = {}
    for isolate in isolates:
        units = {"panel:gyrA": ["gyrA_mut" if isolate in carriers["gyrA_S83L"] else "gyrA_wt"],
                 "panel:parC": ["parC_wt"]}
        if isolate in carriers["qnrS1"]:
            units["acquired"] = ["qnr"]
        mapping[isolate] = units
    return {"mapping": mapping, "panel_loci": ["gyrA", "parC"], "layers": layers}


class EvoFeatureSetTests(unittest.TestCase):
    def test_evo2_feature_sets_learn_from_locus_embeddings_and_record_the_layer(self) -> None:
        cohort, splits, features = small_cohort()
        embeddings = synthetic_embeddings(features, cohort["isolate_id"].tolist())
        learner = MODULE.LEARNERS["ridge_censored"]
        for feature_set in ("evo2", "amr+evo2"):
            predictions, records = MODULE.run_nested(cohort, splits, features, "human_clinical", learner, 7,
                                                     feature_set, embeddings)
            result = EVALUATE.evaluate(predictions, iterations=200, seed=1)["per_antibiotic"]["ciprofloxacin"]
            self.assertGreater(result["essential_agreement"], 0.85, feature_set)
            layers = {r["selected"]["layer"] for r in records}
            self.assertTrue(layers <= {"blocks.6.mlp.l3", "blocks.12.mlp.l3"}, feature_set)
            self.assertEqual(len({(e["layer"], json.dumps(e["params"], sort_keys=True)) for e in records[0]["inner_ea"]}),
                             2 * len(learner.grid))

    def test_ablation_feature_sets_run_and_differ_in_width(self) -> None:
        cohort, splits, features = small_cohort()
        embeddings = synthetic_embeddings(features, cohort["isolate_id"].tolist())
        widths = {}
        for feature_set in ("amr+evo2_panel", "amr+evo2_acquired"):
            _, records = MODULE.run_nested(cohort, splits, features, "human_clinical",
                                           MODULE.LEARNERS["ridge_censored"], 7, feature_set, embeddings)
            widths[feature_set] = records[0]["n_features"] - records[0]["n_amr_features"]
        self.assertEqual(widths["amr+evo2_panel"], 2 * 17)
        self.assertEqual(widths["amr+evo2_acquired"], 33)

    def test_evo2_feature_set_requires_embeddings(self) -> None:
        cohort, splits, features = small_cohort()
        with self.assertRaises(Exception):
            MODULE.run_nested(cohort, splits, features, "human_clinical", MODULE.LEARNERS["ridge_censored"], 7, "evo2")
