from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).parents[1] / "scripts" / "evaluate_predictions.py"
SPEC = importlib.util.spec_from_file_location("evaluate_predictions", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def predictions() -> pd.DataFrame:
    rows = []
    values = {
        "validation": [(0, 0.1), (0, 0.2), (1, 0.8), (1, 0.9)],
        "external": [(0, 0.1), (0, 0.3), (1, 0.7), (1, 0.9)],
    }
    for partition, pairs in values.items():
        for index, (label, score) in enumerate(pairs):
            rows.append(
                {
                    "isolate_id": f"{partition}-{index}",
                    "antibiotic": "ciprofloxacin",
                    "model": "test-model",
                    "partition": partition,
                    "y_true": label,
                    "y_score": score,
                    "genomic_cluster": f"{partition}-cluster-{index}",
                    "lineage_group": f"{partition}-lineage-{index}",
                }
            )
    return pd.DataFrame(rows)


class LockedEvaluationTests(unittest.TestCase):
    def test_threshold_comes_from_validation_and_external_is_reported(self) -> None:
        thresholds, metrics = MODULE.evaluate(
            predictions(), min_sensitivity=1.0, bootstrap_iterations=200, seed=9
        )
        self.assertEqual(thresholds.loc[0, "selection_partition"], "validation")
        self.assertEqual(set(metrics["partition"]), {"external"})
        sensitivity = metrics.loc[metrics["metric"].eq("sensitivity"), "estimate"].iloc[0]
        # The threshold is locked on validation; the lower-scoring external
        # resistant isolate is therefore a real external false-susceptible call.
        self.assertEqual(sensitivity, 0.5)

    def test_duplicate_prediction_fails(self) -> None:
        data = predictions()
        data = pd.concat([data, data.iloc[[0]]], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            MODULE.evaluate(data, bootstrap_iterations=100)

    def test_fractional_binary_labels_fail_before_integer_coercion(self) -> None:
        data = predictions()
        data["y_true"] = data["y_true"].astype(float)
        data.loc[0, "y_true"] = 0.5
        with self.assertRaisesRegex(ValueError, "binary labels"):
            MODULE.evaluate(data, bootstrap_iterations=100)

    def test_model_comparison_requires_identical_isolates(self) -> None:
        data = predictions()
        second = data.copy()
        second["model"] = "other-model"
        second.loc[0, "isolate_id"] = "different-isolate"
        with self.assertRaisesRegex(ValueError, "identical isolates"):
            MODULE.evaluate(pd.concat([data, second], ignore_index=True), bootstrap_iterations=100)


if __name__ == "__main__":
    unittest.main()


class BootstrapUnitTests(unittest.TestCase):
    """Amendment 005: resampling is by lineage, never by chained ANI components."""

    def test_a_single_chained_cluster_does_not_break_the_bootstrap(self) -> None:
        frame = predictions()
        external = frame.loc[frame["partition"].eq("external")].copy()
        external["genomic_cluster"] = "one-giant-chained-component"
        # Resampling ANI components would see a single unit and refuse.
        with self.assertRaises(ValueError):
            MODULE.cluster_bootstrap(external, 0.5, 200, 1, unit="genomic_cluster")
        intervals = MODULE.cluster_bootstrap(external, 0.5, 200, 1)
        self.assertIn("sensitivity", intervals)

    def test_default_unit_and_reported_unit_are_lineage(self) -> None:
        self.assertEqual(MODULE.BOOTSTRAP_UNIT, "lineage_group")
        _, metrics = MODULE.evaluate(predictions(), min_sensitivity=1.0, bootstrap_iterations=200, seed=9)
        self.assertEqual(set(metrics["bootstrap_unit"]), {"lineage_group"})
