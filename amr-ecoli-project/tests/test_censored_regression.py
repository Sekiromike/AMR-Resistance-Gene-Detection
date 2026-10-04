from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from scripts import evaluate_mic_predictions as EVALUATE
from scripts import run_censored_regression as MODULE
from tests_support_mic import simulate, small_cohort


class LikelihoodTests(unittest.TestCase):
    def test_analytic_gradient_matches_finite_differences(self) -> None:
        rng = np.random.default_rng(1)
        X = rng.binomial(1, 0.4, size=(40, 3)).astype(float)
        lo = np.array([-np.inf, 0.0, 3.0, -1.0] * 10)
        hi = np.array([-2.0, 0.0, np.inf, 2.0] * 10)
        theta = np.array([0.3, np.log(1.3), 0.5, -0.7, 1.1])
        _, grad = MODULE.negative_log_likelihood(theta, X, lo, hi, 2.0)
        numeric = np.zeros_like(theta)
        for j in range(len(theta)):
            step = np.zeros_like(theta); step[j] = 1e-6
            numeric[j] = (MODULE.negative_log_likelihood(theta + step, X, lo, hi, 2.0)[0]
                          - MODULE.negative_log_likelihood(theta - step, X, lo, hi, 2.0)[0]) / 2e-6
        np.testing.assert_allclose(grad, numeric, rtol=1e-5, atol=1e-6)

    def test_extreme_tails_stay_finite(self) -> None:
        X = np.zeros((2, 1))
        lo = np.array([40.0, -np.inf]); hi = np.array([40.0, -40.0])
        value, grad = MODULE.negative_log_likelihood(np.array([0.0, 0.0, 0.0]), X, lo, hi, 0.0)
        self.assertTrue(np.isfinite(value) and np.all(np.isfinite(grad)))

    def test_badly_scaled_features_fit_without_overflow_and_sigma_stays_bounded(self) -> None:
        import warnings
        X, lo, hi = simulate(n=800, seed=2)
        X = np.hstack([X, np.random.default_rng(0).normal(0, 300, size=(800, 20))])
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            theta, converged, message = MODULE.fit_checked(X, lo, hi, penalty=0.001)
        self.assertTrue(np.all(np.isfinite(theta)))
        self.assertGreaterEqual(theta[1], MODULE.LOG_SIGMA_BOUNDS[0] - 1e-9)
        self.assertLessEqual(theta[1], MODULE.LOG_SIGMA_BOUNDS[1] + 1e-9)
        self.assertIsInstance(message, str)

    def test_fit_recovers_known_effects_from_censored_data(self) -> None:
        X, lo, hi = simulate()
        theta = MODULE.fit(X, lo, hi, penalty=0.01)
        np.testing.assert_allclose(theta[2:], [3.0, 2.0, 0.0], atol=0.15)
        self.assertAlmostEqual(theta[0], -3.0, delta=0.15)
        self.assertAlmostEqual(float(np.exp(theta[1])), 0.7, delta=0.1)


class NestedCrossValidationTests(unittest.TestCase):
    def test_nested_cv_predicts_every_row_once_and_beats_constant(self) -> None:
        cohort, splits, features = small_cohort()
        predictions, records = MODULE.run_model(cohort, splits, features, "human_clinical", seed=7)
        self.assertEqual(len(predictions), 600)
        self.assertEqual(predictions["isolate_id"].nunique(), 600)
        self.assertEqual(len(records), 5)
        result = EVALUATE.evaluate(predictions, iterations=200, seed=1)["per_antibiotic"]["ciprofloxacin"]
        self.assertGreater(result["essential_agreement"], 0.9)

    def test_test_fold_features_never_enter_the_model(self) -> None:
        cohort, splits, features = small_cohort()
        _, records = MODULE.run_model(cohort, splits, features, "human_clinical", seed=7)
        fold0 = next(r for r in records if r["cv_fold"] == 0)
        self.assertNotIn("fold0_only", [c["feature"] for c in fold0["top_coefficients"]])
        self.assertEqual(fold0["n_features"], 3)
        self.assertEqual(next(r for r in records if r["cv_fold"] == 1)["n_features"], 4)

    def test_inner_folds_keep_groups_whole_and_are_deterministic(self) -> None:
        groups = pd.Series([f"g{i % 13}" for i in range(200)])
        a = MODULE.inner_folds(groups, 4, seed=3)
        b = MODULE.inner_folds(groups, 4, seed=3)
        np.testing.assert_array_equal(a, b)
        self.assertTrue((pd.DataFrame({"g": groups, "f": a}).groupby("g")["f"].nunique() == 1).all())

    def test_external_rows_are_never_modelled(self) -> None:
        cohort, splits, features = small_cohort()
        cohort.loc[0, "evaluation_split"] = "external"
        predictions, _ = MODULE.run_model(cohort, splits, features, "human_clinical", seed=7)
        self.assertNotIn("I0", set(predictions["isolate_id"]))


if __name__ == "__main__":
    unittest.main()
