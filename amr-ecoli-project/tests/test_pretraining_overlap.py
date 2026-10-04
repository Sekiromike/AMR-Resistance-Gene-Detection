from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts import audit_pretraining_overlap as MODULE


class OverlapTests(unittest.TestCase):
    def test_classification_thresholds(self) -> None:
        self.assertEqual(MODULE.classify(99.0, 95.0), "near_identical")
        self.assertEqual(MODULE.classify(98.9, 100.0), "close")
        self.assertEqual(MODULE.classify(99.5, 90.0), "close")
        self.assertEqual(MODULE.classify(89.9, 100.0), "not_found")
        self.assertEqual(MODULE.classify(95.0, 79.0), "not_found")

    def test_best_class_per_unit_and_row_weighting(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "blast.tsv").write_text(
                "h1\tp1\t95.0\t900\t1000\t90\t800\n"
                "h1\tp2\t99.8\t1000\t1000\t100\t1800\n"
                "h2\tp3\t92.0\t1000\t1000\t85\t900\n", encoding="utf-8")
            (root / "map.tsv").write_text(
                "isolate_id\tcohort\tunit\tstatus\tlength\tsequence_sha256\n"
                "A\tdevelopment\tacquired:blaCTX-M-15\tfound\t1000\th1\n"
                "B\tdevelopment\tacquired:blaCTX-M-15\tfound\t1000\th1\n"
                "B\tdevelopment\tacquired:aac(3)-IId\tfound\t1000\th2\n"
                "C\texternal\tacquired:qnrS1#2\tfound\t1000\th3\n"
                "C\texternal\tpanel:gyrA\tfound\t3000\th4\n", encoding="utf-8")
            summary = MODULE.summarize(MODULE.best_hits(root / "blast.tsv"), root / "map.tsv")
        self.assertEqual(summary["unique_acquired_units"], {"near_identical": 1, "close": 1, "not_found": 1})
        self.assertEqual(summary["acquired_unit_rows"], {"near_identical": 2, "close": 1, "not_found": 1})
        self.assertEqual(summary["most_common_genes"]["qnrS1"]["not_found"], 1)
