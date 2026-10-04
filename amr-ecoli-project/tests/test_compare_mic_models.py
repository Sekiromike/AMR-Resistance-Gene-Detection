from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from scripts import compare_mic_models as MODULE


def predictions(correct: list[bool], lineages: list[str]) -> pd.DataFrame:
    return pd.DataFrame({
        "isolate_id": [f"I{i}" for i in range(len(correct))], "antibiotic": "ciprofloxacin",
        "evaluation_split": "development", "lineage_group": lineages, "measurement_sign": "=", "ast_value": 1.0,
        "predicted_log2_mic": [0.0 if ok else 3.0 for ok in correct],
    })


class PairedDifferenceTests(unittest.TestCase):
    def test_difference_and_paired_interval(self) -> None:
        # Contiguous blocks: the gain is concentrated in lineages L6-L8, so it varies by lineage.
        lineages = [f"L{i // 10}" for i in range(100)]
        model = predictions([True] * 90 + [False] * 10, lineages)
        baseline = predictions([True] * 60 + [False] * 40, lineages)
        result = MODULE.paired_difference(model, baseline, iterations=500, seed=1)
        self.assertAlmostEqual(result["difference"], 0.30)
        low, high = result["difference_95ci"]
        # A gain carried by 3 of 10 lineages is uncertain: resampling lineages must show it.
        self.assertLess(low, 0.30)
        self.assertGreater(high, 0.30)
        self.assertGreaterEqual(low, 0.0)

    def test_identical_models_have_a_degenerate_zero_interval(self) -> None:
        lineages = [f"L{i % 7}" for i in range(50)]
        frame = predictions([i % 3 == 0 for i in range(50)], lineages)
        result = MODULE.paired_difference(frame, frame, iterations=200, seed=1)
        self.assertEqual(result["difference"], 0.0)
        self.assertEqual(result["difference_95ci"], [0.0, 0.0])

    def test_rows_must_match(self) -> None:
        lineages = [f"L{i % 5}" for i in range(20)]
        with self.assertRaises(MODULE.ComparisonError):
            MODULE.paired_difference(predictions([True] * 20, lineages), predictions([True] * 19, lineages[:19]),
                                     iterations=100, seed=1)


class SourceRobustnessTests(unittest.TestCase):
    def test_gain_confined_to_one_source_vanishes_when_it_is_excluded(self) -> None:
        lineages = [f"L{i % 10}" for i in range(100)]
        # The model beats the baseline only on isolates I0-I19, all from project P1.
        model = predictions([True] * 100, lineages)
        baseline = predictions([False] * 20 + [True] * 80, lineages)
        sources = pd.Series({f"I{i}": ("P1" if i < 20 else ("P2" if i < 60 else f"small{i}")) for i in range(100)})
        result = MODULE.source_robustness(model, baseline, sources, iterations=300, seed=1)
        self.assertEqual(set(result), {"P1", "P2"})              # sources with >= 10% of rows
        self.assertAlmostEqual(result["P1"]["share_excluded"], 0.20)
        self.assertEqual(result["P1"]["difference"], 0.0)          # gain disappears without P1
        self.assertAlmostEqual(result["P2"]["difference"], 20 / 60)

    def test_single_source_set_is_skipped_not_emptied(self) -> None:
        lineages = [f"L{i % 10}" for i in range(50)]
        model, baseline = predictions([True] * 50, lineages), predictions([False] * 50, lineages)
        sources = pd.Series({f"I{i}": "PRJDB10842" for i in range(50)})
        result = MODULE.source_robustness(model, baseline, sources, iterations=100, seed=1)
        self.assertEqual(result["PRJDB10842"]["skipped"], "too few rows remain")


class PrecisionPlanningTests(unittest.TestCase):
    def test_icc_is_high_when_errors_are_lineage_bound_and_zero_when_not(self) -> None:
        groups = np.repeat(np.arange(20), 10)
        clustered = np.repeat([1.0, 0.0] * 10, 10)
        self.assertGreater(MODULE.anova_icc(clustered, groups), 0.95)
        rng = np.random.default_rng(0)
        self.assertLess(MODULE.anova_icc(rng.binomial(1, 0.1, 200).astype(float), groups), 0.05)

    def test_plan_uses_external_group_sizes_and_labels_power(self) -> None:
        lineages = [f"L{i % 10}" for i in range(200)]
        model = predictions([i % 10 != 0 for i in range(200)], lineages)  # all errors in lineage L0
        many_small = pd.Series([f"E{i}" for i in range(1500)])
        one_big = pd.Series(["E0"] * 1500)
        good = MODULE.plan_precision(model, many_small, 0.9, 0.05)
        bad = MODULE.plan_precision(model, one_big, 0.9, 0.05)
        self.assertAlmostEqual(good["design_effect"], 1.0)
        self.assertEqual(good["status"], "adequate")
        self.assertEqual(bad["status"], "underpowered")
        self.assertGreater(bad["development_error_icc"], 0.9)


if __name__ == "__main__":
    unittest.main()
