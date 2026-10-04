from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location("normalize_jarbs_ast", ROOT / "scripts" / "normalize_jarbs_ast.py")
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
VOCAB = ROOT / "config" / "insdc_geo_loc_name_vocabulary.tsv"
FIELDS = ["isolate_id", "genome_species", "biosample_accession", "primary_short_read_run", "collection_date_raw",
          "geo_loc_name_raw", "isolation_source_raw", "ceftriaxone_mic_raw", "ciprofloxacin_mic_raw", "gentamicin_mic_raw"]


class JarbsNormalizationTests(unittest.TestCase):
    def _run(self, rows):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "jarbs.csv"
            with path.open("w", encoding="utf-8", newline="") as h:
                w = csv.DictWriter(h, fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
            return MODULE.normalize_jarbs(path, VOCAB, "test-release")

    @staticmethod
    def row(isolate="JB-1", cro="<=1", cip=">2", gen="4", date="2019-05-01", species="Escherichia coli"):
        return {"isolate_id": isolate, "genome_species": species, "biosample_accession": "SAMD00000001",
                "primary_short_read_run": "DRR000001", "collection_date_raw": date,
                "geo_loc_name_raw": "Japan", "isolation_source_raw": "urine",
                "ceftriaxone_mic_raw": cro, "ciprofloxacin_mic_raw": cip, "gentamicin_mic_raw": gen}

    def test_signs_values_and_study_level_fields(self) -> None:
        frame, counts = self._run([self.row()])
        by_drug = frame.set_index("antibiotic")
        self.assertEqual((by_drug.loc["ceftriaxone", "measurement_sign"], by_drug.loc["ceftriaxone", "ast_value"]), ("<=", "1"))
        self.assertEqual((by_drug.loc["ciprofloxacin", "measurement_sign"], by_drug.loc["ciprofloxacin", "ast_value"]), (">", "2"))
        self.assertEqual((by_drug.loc["gentamicin", "measurement_sign"], by_drug.loc["gentamicin", "ast_value"]), ("=", "4"))
        self.assertTrue((frame["country"] == "Japan").all())
        self.assertTrue((frame["country_source"] == "insdc_geo_loc_name_prefix").all())
        self.assertTrue((frame["assembly_accession"] == "").all())
        self.assertTrue((frame["ast_testing_date"] == "").all())  # not published; never invented
        self.assertEqual(counts["mic_eligible_rows"], 3)

    def test_untested_drug_emits_no_row(self) -> None:
        frame, counts = self._run([self.row(cro="", cip="", gen="")])
        self.assertEqual(len(frame), 0)
        self.assertEqual(counts["not_tested_by_drug"], {"ceftriaxone": 1, "ciprofloxacin": 1, "gentamicin": 1})

    def test_unparseable_value_is_ineligible_not_dropped(self) -> None:
        frame, counts = self._run([self.row(cro="see note")])
        bad = frame[frame["antibiotic"] == "ceftriaxone"].iloc[0]
        self.assertEqual(bad["mic_eligible"], "false")
        self.assertEqual(bad["mic_ineligibility_reasons"], "unparseable_mic")

    def test_missing_date_token_becomes_blank(self) -> None:
        frame, _ = self._run([self.row(date="missing")])
        self.assertTrue((frame["collection_date"] == "").all())

    def test_record_ids_are_unique_and_traceable(self) -> None:
        frame, _ = self._run([self.row("JB-1"), self.row("JB-2")])
        self.assertFalse(frame["source_ast_record_id"].duplicated().any())
        self.assertIn("JARBS:JB-2:gentamicin", set(frame["source_ast_record_id"]))

    def test_non_ecoli_species_fails_closed(self) -> None:
        with self.assertRaises(MODULE.JarbsNormalizationError):
            self._run([self.row(species="Klebsiella pneumoniae")])

    def test_duplicate_isolate_fails_closed(self) -> None:
        with self.assertRaises(MODULE.JarbsNormalizationError):
            self._run([self.row("JB-1"), self.row("JB-1")])
