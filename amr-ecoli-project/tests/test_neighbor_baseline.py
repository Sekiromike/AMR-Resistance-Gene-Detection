from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts import run_neighbor_baseline as MODULE

HEADER = "development_sequence_id\texternal_sequence_id\tani_percent\taligned_fraction_development\taligned_fraction_external\n"


def comparisons(root: Path, pairs: list[tuple[str, str, float, float]]) -> Path:
    path = root / "comparisons.tsv"
    lines = [HEADER]
    for a, b, ani, af in pairs:
        for q, r in ((a, b), (b, a)):
            lines.append(f"development::{q}\tdevelopment::{r}\t{ani}\t{af}\t{af}\n")
    for isolate in {x for pair in pairs for x in pair[:2]}:
        lines.append(f"development::{isolate}\tdevelopment::{isolate}\t100\t1\t1\n")
    path.write_text("".join(lines), encoding="utf-8")
    return path


def cohort() -> pd.DataFrame:
    rows = []
    spec = {"A": ("human_clinical", "=", 1.0), "B": ("human_clinical", "<=", 0.25),
            "C": ("non_human_or_environmental", ">", 4.0), "D": ("human_clinical", ">=", 8.0)}
    for isolate, (pop, sign, value) in spec.items():
        rows.append({"isolate_id": isolate, "antibiotic": "ciprofloxacin", "evaluation_split": "development",
                     "lineage_group": f"SLV:{isolate}", "intended_use_population": pop,
                     "measurement_sign": sign, "ast_value": value})
    return pd.DataFrame(rows)


FOLD = {"A": 0, "B": 1, "C": 2, "D": 0}


class NeighborBaselineTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def test_prediction_is_the_reference_set_boundary(self) -> None:
        self.assertEqual(MODULE.boundary_prediction("=", 1.0), 0.0)
        self.assertEqual(MODULE.boundary_prediction("<=", 0.25), -2.0)
        self.assertEqual(MODULE.boundary_prediction("<", 0.25), -3.0)
        self.assertEqual(MODULE.boundary_prediction(">=", 8.0), 3.0)
        self.assertEqual(MODULE.boundary_prediction(">", 4.0), 3.0)

    def test_same_fold_relatives_are_never_used_and_population_is_respected(self) -> None:
        # A's closest genome is D, but D shares A's fold. Next is C (not human clinical), then B.
        path = comparisons(self.root, [("A", "D", 99.99, 0.99), ("A", "C", 99.5, 0.95), ("A", "B", 99.0, 0.95),
                                       ("B", "C", 98.0, 0.9), ("C", "D", 97.0, 0.9), ("B", "D", 96.0, 0.9)])
        clinical = {"A", "B", "D"}
        neighbors = MODULE.collect_neighbors(path, FOLD, clinical)
        self.assertEqual([r for _, _, r in neighbors["all_sources"]["A"]], ["C", "B"])
        self.assertEqual([r for _, _, r in neighbors["human_clinical"]["A"]], ["B"])
        clin = MODULE.predict(cohort(), FOLD, neighbors, "human_clinical").set_index("isolate_id")
        self.assertEqual(clin.loc["A", "neighbor_isolate_id"], "B")
        self.assertEqual(clin.loc["A", "predicted_log2_mic"], -2.0)   # B is "<= 0.25"
        allsrc = MODULE.predict(cohort(), FOLD, neighbors, "all_sources").set_index("isolate_id")
        self.assertEqual(allsrc.loc["A", "neighbor_isolate_id"], "C")
        self.assertEqual(allsrc.loc["A", "predicted_log2_mic"], 3.0)  # C is "> 4"

    def test_ties_break_on_aligned_fraction_then_isolate_id(self) -> None:
        path = comparisons(self.root, [("A", "C", 99.0, 0.95), ("A", "B", 99.0, 0.95)])
        neighbors = MODULE.collect_neighbors(path, FOLD, {"A", "B", "C"})
        self.assertEqual([r for _, _, r in neighbors["all_sources"]["A"]], ["B", "C"])

    def test_isolate_without_any_relative_fails_closed(self) -> None:
        path = comparisons(self.root, [("A", "B", 99.0, 0.95)])
        neighbors = MODULE.collect_neighbors(path, FOLD, {"A", "B", "D"})
        with self.assertRaisesRegex(MODULE.NeighborError, "no training relative"):
            MODULE.predict(cohort(), FOLD, neighbors, "all_sources")
