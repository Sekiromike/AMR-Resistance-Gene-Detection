from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

MODULE_PATH = Path(__file__).parents[1] / "scripts" / "build_source_attributes.py"
SPEC = importlib.util.spec_from_file_location("build_source_attributes", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ClassificationTests(unittest.TestCase):
    """One test per row of the amendment 006 classification table."""

    def test_clinical_human_is_human_clinical(self) -> None:
        self.assertEqual(MODULE.classify_development("clinical", "Homo sapiens")[0], "human_clinical")

    def test_clinical_without_host_is_undetermined_not_assumed_human(self) -> None:
        self.assertEqual(MODULE.classify_development("clinical", "")[0], "undetermined")

    def test_veterinary_clinical_is_not_human(self) -> None:
        population, basis = MODULE.classify_development("clinical", "Canis lupus familiaris")
        self.assertEqual(population, "non_human_or_environmental")
        self.assertIn("Canis lupus familiaris", basis)

    def test_environmental_with_human_host_is_a_conflict(self) -> None:
        self.assertEqual(MODULE.classify_development("environmental/other", "Homo sapiens")[0], "undetermined")

    def test_environmental_is_not_human(self) -> None:
        self.assertEqual(MODULE.classify_development("environmental/other", "")[0], "non_human_or_environmental")

    def test_absent_epi_type_is_undetermined(self) -> None:
        self.assertEqual(MODULE.classify_development("", "Homo sapiens")[0], "undetermined")

    def test_host_match_is_exact_and_case_insensitive(self) -> None:
        self.assertEqual(MODULE.classify_development("clinical", "homo SAPIENS")[0], "human_clinical")
        self.assertNotEqual(MODULE.classify_development("clinical", "human-associated")[0], "human_clinical")


class BuildTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def _write(self, name, fields, rows):
        with (self.root / name).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def _fixture(self, development_rows, jarbs_rows):
        self._write("dev.csv", ["isolate_id", "antibiotic", "epi_type", "host", "isolation_source"], development_rows)
        self._write("jarbs.csv", ["isolate_id", "isolation_source_raw"], jarbs_rows)

    def _build(self):
        return MODULE.build_source_attributes(
            self.root / "dev.csv", self.root / "jarbs.csv", self.root / "out.tsv", self.root / "m.json"
        )

    def test_counts_and_external_basis(self) -> None:
        self._fixture(
            [
                {"isolate_id": "D1", "antibiotic": "cip", "epi_type": "clinical", "host": "Homo sapiens", "isolation_source": "urine"},
                {"isolate_id": "D1", "antibiotic": "gen", "epi_type": "clinical", "host": "Homo sapiens", "isolation_source": "urine"},
                {"isolate_id": "D2", "antibiotic": "cip", "epi_type": "environmental/other", "host": "", "isolation_source": "ground turkey"},
                {"isolate_id": "D3", "antibiotic": "cip", "epi_type": "clinical", "host": "", "isolation_source": "not collected"},
            ],
            [{"isolate_id": "JB-1", "isolation_source_raw": "venous blood"}],
        )
        payload = self._build()
        self.assertEqual(payload["population_counts"]["development"],
                         {"human_clinical": 1, "non_human_or_environmental": 1, "undetermined": 1})
        self.assertEqual(payload["population_counts"]["external"], {"human_clinical": 1})
        self.assertEqual(payload["specimen_source_absent"]["development"], 1)  # "not collected"
        rows = {r["isolate_id"]: r for r in csv.DictReader((self.root / "out.tsv").open(encoding="utf-8"), delimiter="\t")}
        self.assertEqual(len(rows), 4)  # one row per isolate, not per AST row
        self.assertIn("study design", rows["JB-1"]["population_basis"])
        self.assertFalse(payload["scientific_boundary"]["phenotypes_read"])

    def test_conflicting_values_across_an_isolates_rows_fail_closed(self) -> None:
        self._fixture(
            [
                {"isolate_id": "D1", "antibiotic": "cip", "epi_type": "clinical", "host": "Homo sapiens", "isolation_source": "urine"},
                {"isolate_id": "D1", "antibiotic": "gen", "epi_type": "clinical", "host": "Bos taurus", "isolation_source": "urine"},
            ],
            [{"isolate_id": "JB-1", "isolation_source_raw": "urine"}],
        )
        with self.assertRaisesRegex(MODULE.AttributeError_, "conflicting host"):
            self._build()

    def test_identifier_in_both_cohorts_fails_closed(self) -> None:
        self._fixture(
            [{"isolate_id": "X", "antibiotic": "cip", "epi_type": "clinical", "host": "Homo sapiens", "isolation_source": "urine"}],
            [{"isolate_id": "X", "isolation_source_raw": "urine"}],
        )
        with self.assertRaisesRegex(MODULE.AttributeError_, "both cohorts"):
            self._build()
