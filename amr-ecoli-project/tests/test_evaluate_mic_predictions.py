from __future__ import annotations

import math
import unittest

import numpy as np
import pandas as pd

from scripts import evaluate_mic_predictions as MODULE


def frame(rows: list[tuple[str, float, float]], *, split: str = "development", lineages=None) -> pd.DataFrame:
    """rows: (sign, reference mg/L, predicted log2)."""
    return pd.DataFrame({
        "isolate_id": [f"I{i}" for i in range(len(rows))],
        "antibiotic": "ciprofloxacin",
        "evaluation_split": split,
        "lineage_group": lineages or [f"ST{i % 4}" for i in range(len(rows))],
        "measurement_sign": [r[0] for r in rows],
        "ast_value": [r[1] for r in rows],
        "predicted_log2_mic": [r[2] for r in rows],
    })


class DilutionTests(unittest.TestCase):
    def test_conventional_labels_snap_and_off_scale_values_are_refused(self) -> None:
        self.assertEqual(MODULE.dilution_index(0.12), -3)
        self.assertEqual(MODULE.dilution_index(0.06), -4)
        self.assertEqual(MODULE.dilution_index(0.015), -6)
        self.assertEqual(MODULE.dilution_index(0.008), -7)
        self.assertEqual(MODULE.dilution_index(4), 2)
        for value in (0.2, 2.5, 200, 3):  # the three off-scale development values, and 3
            with self.assertRaises(MODULE.MicEvaluationError):
                MODULE.dilution_index(value)
        for value in (0, -1, float("nan")):
            with self.assertRaises(MODULE.MicEvaluationError):
                MODULE.dilution_index(value)

    def test_reference_sets_for_each_comparator(self) -> None:
        self.assertEqual(MODULE.reference_bounds("=", 1), (0, 0))
        self.assertEqual(MODULE.reference_bounds("<=", 0.25), (-math.inf, -2))
        self.assertEqual(MODULE.reference_bounds("<", 0.25), (-math.inf, -3))
        self.assertEqual(MODULE.reference_bounds(">=", 4), (2, math.inf))
        self.assertEqual(MODULE.reference_bounds(">", 4), (3, math.inf))
        with self.assertRaises(MODULE.MicEvaluationError):
            MODULE.reference_bounds("~", 1)

    def test_rounding_is_nearest_with_halves_up(self) -> None:
        np.testing.assert_array_equal(MODULE.predicted_dilution(np.array([-2.5, -0.49, 0.5, 1.49])),
                                      np.array([-2, 0, 1, 1]))


class ScoringTests(unittest.TestCase):
    def test_exact_reference_agreement(self) -> None:
        scored = MODULE.score_rows(frame([("=", 1, 0.0), ("=", 1, 1.0), ("=", 1, -2.0)]))
        self.assertEqual(scored["error_dilutions"].tolist(), [0, 1, 2])
        self.assertEqual(scored["essential_agreement"].tolist(), [True, True, False])
        self.assertEqual(scored["exact_agreement"].tolist(), [True, False, False])
        self.assertEqual(scored["under_call_2plus"].tolist(), [False, False, True])

    def test_censored_reference_agreement(self) -> None:
        # "<=0.25" (d <= -2): predictions -5, -2 consistent; -1 within one; 0 two above.
        scored = MODULE.score_rows(frame([("<=", 0.25, -5.0), ("<=", 0.25, -2.0),
                                          ("<=", 0.25, -1.0), ("<=", 0.25, 0.0)]))
        self.assertEqual(scored["error_dilutions"].tolist(), [0, 0, 1, 2])
        self.assertEqual(scored["over_call_2plus"].tolist(), [False, False, False, True])
        # ">4" means d >= 3 (8 mg/L or more): prediction 2 (4 mg/L) is one step off.
        scored = MODULE.score_rows(frame([(">", 4, 2.0), (">", 4, 1.0), (">", 4, 6.0)]))
        self.assertEqual(scored["error_dilutions"].tolist(), [1, 2, 0])
        self.assertEqual(scored["under_call_2plus"].tolist(), [False, True, False])

    def test_summary_reports_all_rows_and_exact_only(self) -> None:
        scored = MODULE.score_rows(frame([("=", 1, 0.0), ("=", 1, 3.0), ("<=", 0.25, -4.0), ("<=", 0.25, -3.0)]))
        summary = MODULE.summarize(scored)
        self.assertEqual(summary["n_exact_reference"], 2)
        self.assertAlmostEqual(summary["essential_agreement"], 0.75)
        self.assertAlmostEqual(summary["essential_agreement_exact_only"], 0.5)
        self.assertAlmostEqual(summary["bias_exact_only"], 1.5)

    def test_contract_violations(self) -> None:
        f = frame([("=", 1, 0.0), ("=", 1, 0.0)])
        f.loc[1, "isolate_id"] = "I0"
        with self.assertRaisesRegex(MODULE.MicEvaluationError, "duplicate"):
            MODULE.score_rows(f)
        f = frame([("=", 1, float("nan"))])
        with self.assertRaisesRegex(MODULE.MicEvaluationError, "prediction"):
            MODULE.score_rows(f)
        with self.assertRaisesRegex(MODULE.MicEvaluationError, "columns"):
            MODULE.score_rows(frame([("=", 1, 0.0)]).drop(columns="lineage_group"))


class LikelihoodTests(unittest.TestCase):
    def test_interval_likelihood_matches_normal_probabilities(self) -> None:
        # Exact d=0 -> latent (-1, 0]; censored <=-2 -> (-inf, -2]; >=2 -> (1, inf).
        ll = MODULE.interval_log_likelihood(np.array([0, -np.inf, 2]), np.array([0, -2, np.inf]),
                                            np.zeros(3), np.ones(3))
        expected = np.log([0.5 - 0.15865525393145707, 0.022750131948179195, 0.15865525393145707])
        np.testing.assert_allclose(ll, expected, rtol=1e-9)

    def test_far_upper_tail_stays_finite_and_ordered(self) -> None:
        ll = MODULE.interval_log_likelihood(np.array([30.0, 40.0]), np.array([30.0, 40.0]),
                                            np.zeros(2), np.ones(2))
        self.assertTrue(np.all(np.isfinite(ll)))

    def test_empty_sd_column_means_point_prediction_but_partial_is_refused(self) -> None:
        # Regression: external scoring job 65020797 stopped on baseline files
        # whose sd column is entirely empty.
        f = frame([("=", 1, 0.0), ("<=", 0.25, -3.0)])
        f["predicted_log2_sd"] = [float("nan"), float("nan")]
        self.assertNotIn("mean_interval_log_likelihood", MODULE.summarize(MODULE.score_rows(f)))
        f["predicted_log2_sd"] = [1.0, float("nan")]
        with self.assertRaises(MODULE.MicEvaluationError):
            MODULE.score_rows(f)

    def test_sd_column_adds_log_likelihood_and_must_be_positive(self) -> None:
        f = frame([("=", 1, 0.0), ("<=", 0.25, -3.0)])
        f["predicted_log2_sd"] = [1.0, 1.0]
        self.assertIn("mean_interval_log_likelihood", MODULE.summarize(MODULE.score_rows(f)))
        f["predicted_log2_sd"] = [1.0, 0.0]
        with self.assertRaises(MODULE.MicEvaluationError):
            MODULE.score_rows(f)


class EvaluateTests(unittest.TestCase):
    def test_external_rows_are_refused_unless_locked(self) -> None:
        f = frame([("=", 1, 0.0)] * 4, split="external")
        with self.assertRaisesRegex(MODULE.MicEvaluationError, "locked-external"):
            MODULE.evaluate(f, iterations=200, seed=1)
        mixed = pd.concat([frame([("=", 1, 0.0)] * 4), f.assign(isolate_id=lambda d: "E" + d["isolate_id"])])
        with self.assertRaises(MODULE.MicEvaluationError):
            MODULE.evaluate(mixed, iterations=200, seed=1, locked_external=True)

    def test_lineage_bootstrap_is_deterministic_and_brackets_the_estimate(self) -> None:
        rows = [("=", 1, 0.0)] * 30 + [("=", 1, 3.0)] * 10
        lineages = [f"ST{i % 8}" for i in range(40)]
        a = MODULE.evaluate(frame(rows, lineages=lineages), iterations=500, seed=7)["per_antibiotic"]
        b = MODULE.evaluate(frame(rows, lineages=lineages), iterations=500, seed=7)["per_antibiotic"]
        self.assertEqual(a, b)
        result = a["ciprofloxacin"]
        low, high = result["bootstrap_95ci"]["essential_agreement"]
        self.assertLessEqual(low, result["essential_agreement"])
        self.assertGreaterEqual(high, result["essential_agreement"])
        self.assertEqual(result["bootstrap_95ci"]["groups"], 8)

    def test_whole_lineage_resampling_widens_intervals_for_clonal_errors(self) -> None:
        # All errors sit in one lineage: resampling lineages must reflect that.
        clonal = ["ST0"] * 4 + [f"ST{i}" for i in range(1, 37)]
        rows = [("=", 1, 3.0)] * 4 + [("=", 1, 0.0)] * 36
        low, high = MODULE.evaluate(frame(rows, lineages=clonal), iterations=2000, seed=3)[
            "per_antibiotic"]["ciprofloxacin"]["bootstrap_95ci"]["essential_agreement"]
        self.assertLess(low, 0.9)
        self.assertEqual(high, 1.0)


class ReferenceClassTests(unittest.TestCase):
    def test_classes_and_the_resistant_under_call(self) -> None:
        # low: "<=0.25"; exact: "=1"; high: ">4" (d >= 3).
        rows = [("<=", 0.25, -4.0)] * 4 + [("=", 1, 0.0)] * 2 + [(">", 4, 5.0)] * 3 + [(">", 4, 0.0)]
        lineages = [f"ST{i % 5}" for i in range(10)]
        result = MODULE.evaluate(frame(rows, lineages=lineages), iterations=300, seed=1)[
            "per_antibiotic"]["ciprofloxacin"]["by_reference_class"]
        self.assertEqual({k: v["n"] for k, v in result.items()},
                         {"left_censored": 4, "exact": 2, "right_censored": 4})
        high = result["right_censored"]
        self.assertAlmostEqual(high["essential_agreement"], 0.75)
        self.assertAlmostEqual(high["under_call_2plus"], 0.25)
        low, up = high["under_call_2plus_95ci"]
        self.assertLessEqual(low, 0.25)
        self.assertGreaterEqual(up, 0.25)
        self.assertNotIn("under_call_2plus_95ci", result["exact"])

    def test_constant_low_prediction_gets_all_left_and_no_right_censored_credit(self) -> None:
        rows = [("<=", 0.25, -5.0)] * 6 + [(">", 4, -5.0)] * 4
        result = MODULE.evaluate(frame(rows), iterations=200, seed=1)["per_antibiotic"]["ciprofloxacin"]
        self.assertEqual(result["by_reference_class"]["left_censored"]["essential_agreement"], 1.0)
        self.assertEqual(result["by_reference_class"]["right_censored"]["under_call_2plus"], 1.0)


class StratifiedReportingTests(unittest.TestCase):
    def test_dominant_lineage_and_the_rest_are_reported_separately(self) -> None:
        rows = [("=", 1, 0.0)] * 10 + [("=", 1, 3.0)] * 10
        lineages = ["SLV:ST131"] * 10 + [f"SLV:ST{i}" for i in range(10)]
        f = frame(rows, lineages=lineages)
        result = MODULE.evaluate(f, iterations=300, seed=2, dominant_lineage="SLV:ST131")[
            "per_antibiotic"]["ciprofloxacin"]
        self.assertAlmostEqual(result["essential_agreement"], 0.5)
        self.assertAlmostEqual(result["strata"]["SLV:ST131"]["essential_agreement"], 1.0)
        self.assertEqual(result["strata"]["SLV:ST131"]["bootstrap_unit"], "isolate (single lineage)")
        self.assertAlmostEqual(result["strata"]["other_lineages"]["essential_agreement"], 0.0)
        self.assertEqual(result["strata"]["other_lineages"]["bootstrap_unit"], "lineage_group")

    def test_absent_dominant_lineage_is_reported_empty(self) -> None:
        result = MODULE.evaluate(frame([("=", 1, 0.0)] * 8), iterations=200, seed=2,
                                 dominant_lineage="SLV:ST131")["per_antibiotic"]["ciprofloxacin"]
        self.assertEqual(result["strata"]["SLV:ST131"], {"n": 0})


if __name__ == "__main__":
    unittest.main()
