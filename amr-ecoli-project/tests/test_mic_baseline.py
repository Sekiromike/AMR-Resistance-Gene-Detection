from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from scripts import evaluate_mic_predictions as EVALUATE
from scripts import run_mic_baseline as MODULE


def cohort_and_splits(n: int = 200, seed: int = 0):
    """Latent log2 MIC ~ N(-1, 1.5); panel reads 0.25..8 mg/L, censored outside."""
    rng = np.random.default_rng(seed)
    latent = rng.normal(-1.0, 1.5, size=n)
    reported = np.ceil(latent)  # a dilution d covers (d-1, d]
    rows = []
    for i, d in enumerate(reported):
        if d <= -2:
            sign, value = "<=", 0.25
        elif d > 3:
            sign, value = ">", 8.0
        else:
            sign, value = "=", float(2.0 ** d)
        rows.append({"isolate_id": f"I{i}", "antibiotic": "ciprofloxacin", "evaluation_split": "development",
                     "lineage_group": f"SLV:ST{i % 25}", "measurement_sign": sign, "ast_value": value,
                     "intended_use_population": "human_clinical" if i % 4 else "non_human_or_environmental"})
    rows.append({"isolate_id": "EXT", "antibiotic": "ciprofloxacin", "evaluation_split": "external",
                 "lineage_group": "SLV:ST1", "measurement_sign": "=", "ast_value": 1.0,
                 "intended_use_population": "human_clinical"})
    cohort = pd.DataFrame(rows)
    splits = pd.DataFrame({"isolate_id": cohort["isolate_id"], "evaluation_split": cohort["evaluation_split"],
                           "cv_fold": [i % 5 if s == "development" else pd.NA
                                       for i, s in enumerate(cohort["evaluation_split"])]})
    return cohort, splits


class BaselineTests(unittest.TestCase):
    def test_best_constant_maximises_training_agreement(self) -> None:
        # 6 rows "<=0.25" (d<=-2), 3 rows exact 1 mg/L (d=0), 1 row ">8" (d>=4).
        lo = np.array([-np.inf] * 6 + [0.0] * 3 + [4.0])
        hi = np.array([-2.0] * 6 + [0.0] * 3 + [np.inf])
        dilution, ea, ea_exact = MODULE.best_constant_dilution(lo, hi)
        self.assertEqual(dilution, -1)       # within one of both "<=0.25" and exact 1 mg/L
        self.assertAlmostEqual(ea, 0.9)
        self.assertAlmostEqual(ea_exact, 1.0)

    def test_ties_prefer_exact_agreement_then_the_lower_dilution(self) -> None:
        lo = np.array([-np.inf, 3.0])
        hi = np.array([-2.0, np.inf])
        dilution, ea, _ = MODULE.best_constant_dilution(lo, hi)
        self.assertEqual((dilution, ea), (-3, 0.5))
        with self.assertRaises(MODULE.BaselineError):
            MODULE.best_constant_dilution(np.array([-np.inf]), np.array([np.inf]))

    def test_out_of_fold_predictions_never_use_the_held_out_fold_or_external(self) -> None:
        cohort, splits = cohort_and_splits()
        predictions, fits = MODULE.run_baseline(cohort, splits, "all_sources")
        self.assertNotIn("EXT", set(predictions["isolate_id"]))
        self.assertEqual(len(predictions), 200)
        self.assertEqual(len(fits), 5)
        for fit in fits:
            self.assertEqual(fit["n_train"] + fit["n_test"], 200)
        # Each fold's prediction equals a fit that excluded that fold.
        by_fold = {fit["cv_fold"]: fit["dilution_log2"] for fit in fits}
        for fold, group in predictions.groupby("cv_fold"):
            self.assertTrue(np.allclose(group["predicted_log2_mic"], by_fold[fold]))

    def test_primary_population_filter(self) -> None:
        cohort, splits = cohort_and_splits()
        predictions, _ = MODULE.run_baseline(cohort, splits, "human_clinical")
        self.assertEqual(set(predictions["intended_use_population"]), {"human_clinical"})
        self.assertEqual(len(predictions), 150)

    def test_predictions_feed_the_evaluator(self) -> None:
        cohort, splits = cohort_and_splits()
        predictions, _ = MODULE.run_baseline(cohort, splits, "all_sources")
        result = EVALUATE.evaluate(predictions, iterations=200, seed=1)["per_antibiotic"]["ciprofloxacin"]
        self.assertNotIn("mean_interval_log_likelihood", result)  # point baseline, no sd
        self.assertGreater(result["essential_agreement"], 0.5)

    def test_missing_fold_is_refused(self) -> None:
        cohort, splits = cohort_and_splits()
        splits.loc[0, "cv_fold"] = pd.NA
        with self.assertRaises(MODULE.BaselineError):
            MODULE.run_baseline(cohort, splits, "all_sources")


if __name__ == "__main__":
    unittest.main()
