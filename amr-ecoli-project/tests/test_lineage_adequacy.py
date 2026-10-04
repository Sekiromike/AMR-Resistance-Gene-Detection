from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts import assess_lineage_adequacy as MODULE

ROOT = Path(__file__).parents[1]
CRITERIA = {
    "min_typing_coverage": 0.95,
    "max_largest_lineage_fraction": 0.40,
    "min_inverse_simpson": 10,
    "folds": 5,
}


def diverse(n_lineages: int, per: int) -> list[str]:
    return [f"ST{i}" for i in range(n_lineages) for _ in range(per)]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class MeasureTests(unittest.TestCase):
    def test_diverse_fully_typed_population_keeps_st(self) -> None:
        result = MODULE.measure(diverse(20, 5), 100, CRITERIA)
        self.assertEqual(result["outcome"], "keep_st")
        self.assertAlmostEqual(result["metrics"]["inverse_simpson"], 20.0)
        self.assertEqual(result["metrics"]["lineages_covering_half"], 10)

    def test_coverage_boundary_is_inclusive(self) -> None:
        lineages = diverse(19, 5) + [""] * 5  # exactly 95 of 100 resolved
        self.assertEqual(MODULE.measure(lineages, 100, CRITERIA)["outcome"], "keep_st")
        lineages = diverse(19, 5)[:-1] + [""] * 6  # 94 of 100
        self.assertEqual(MODULE.measure(lineages, 100, CRITERIA)["outcome"], "escalate_cgmlst")

    def test_escalation_takes_precedence_over_concentration(self) -> None:
        lineages = ["ST131"] * 50 + [""] * 50
        self.assertEqual(MODULE.measure(lineages, 100, CRITERIA)["outcome"], "escalate_cgmlst")

    def test_dominant_lineage_reopens(self) -> None:
        lineages = ["ST131"] * 41 + diverse(59, 1)
        result = MODULE.measure(lineages, 100, CRITERIA)
        self.assertFalse(result["checks"]["concentration"])
        self.assertEqual(result["outcome"], "reopen_amendment_004")

    def test_low_effective_diversity_reopens(self) -> None:
        # Nine equal lineages: largest 11% passes, inverse Simpson 9 fails.
        result = MODULE.measure(diverse(9, 10), 90, CRITERIA)
        self.assertTrue(result["checks"]["concentration"])
        self.assertFalse(result["checks"]["effective_diversity"])
        self.assertEqual(result["outcome"], "reopen_amendment_004")

    def test_lineage_above_fold_share_is_reported_not_failed(self) -> None:
        lineages = ["ST131"] * 25 + diverse(75, 1)
        result = MODULE.measure(lineages, 100, CRITERIA)
        self.assertEqual(result["outcome"], "keep_st")
        self.assertTrue(result["metrics"]["largest_lineage_exceeds_fold_share"])

    def test_diversity_boundary_is_exact(self) -> None:
        # Ten equal lineages give inverse Simpson exactly 10, which passes.
        result = MODULE.measure(diverse(10, 7), 70, CRITERIA)
        self.assertEqual(result["metrics"]["inverse_simpson"], 10.0)
        self.assertEqual(result["outcome"], "keep_st")

    def test_empty_population_escalates(self) -> None:
        self.assertEqual(MODULE.measure([], 0, CRITERIA)["outcome"], "escalate_cgmlst")


class AssessTests(unittest.TestCase):
    def _write(self, root: Path, genomes: list[dict[str, str]], attributes: list[dict[str, str]]) -> None:
        write_tsv(root / "genomes.tsv", genomes)
        write_tsv(root / "attributes.tsv", attributes)

    def test_uses_only_eligible_development_and_worst_population_wins(self) -> None:
        genomes, attributes = [], []
        # Human clinical: diverse. Non-human: all one ST, which makes the
        # all-sources population too concentrated.
        for i in range(40):
            genomes.append({"isolate_id": f"H{i}", "cohort": "development",
                            "lineage_group": f"ST{i % 20}", "eligible_for_modeling": "true"})
            attributes.append({"isolate_id": f"H{i}", "cohort": "development",
                               "intended_use_population": "human_clinical"})
        for i in range(60):
            genomes.append({"isolate_id": f"N{i}", "cohort": "development",
                            "lineage_group": "ST10", "eligible_for_modeling": "true"})
            attributes.append({"isolate_id": f"N{i}", "cohort": "development",
                               "intended_use_population": "non_human_or_environmental"})
        # Ineligible and external genomes never count, even when unresolved.
        genomes.append({"isolate_id": "X", "cohort": "development", "lineage_group": "",
                        "eligible_for_modeling": "false"})
        genomes.extend({"isolate_id": f"E{i}", "cohort": "external", "lineage_group": "",
                        "eligible_for_modeling": "true"} for i in range(50))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._write(root, genomes, attributes)
            result = MODULE.assess(root / "genomes.tsv", root / "attributes.tsv", CRITERIA)
        self.assertEqual(result["eligible_development_genomes"], 100)
        self.assertEqual(result["populations"]["human_clinical"]["outcome"], "keep_st")
        self.assertEqual(result["populations"]["all_sources"]["outcome"], "reopen_amendment_004")
        self.assertEqual(result["outcome"], "reopen_amendment_004")

    def test_missing_source_attributes_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._write(
                root,
                [{"isolate_id": "A", "cohort": "development", "lineage_group": "ST1",
                  "eligible_for_modeling": "true"}],
                [{"isolate_id": "B", "cohort": "development",
                  "intended_use_population": "human_clinical"}],
            )
            with self.assertRaises(MODULE.LineageAdequacyError):
                MODULE.assess(root / "genomes.tsv", root / "attributes.tsv", CRITERIA)

    def test_study_config_carries_the_adopted_thresholds(self) -> None:
        criteria = MODULE.load_criteria(ROOT / "config" / "study.json")
        self.assertEqual(criteria["min_typing_coverage"], 0.95)
        self.assertEqual(criteria["max_largest_lineage_fraction"], 0.40)
        self.assertEqual(criteria["min_inverse_simpson"], 10)
        self.assertEqual(criteria["folds"], 5)
        config = json.loads((ROOT / "config" / "study.json").read_text(encoding="utf-8"))
        self.assertEqual(config["splitting"]["group"], "lineage_group")


if __name__ == "__main__":
    unittest.main()
