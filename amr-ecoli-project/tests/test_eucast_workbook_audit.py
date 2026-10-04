from __future__ import annotations

import importlib.util
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).parents[1]
MODULE_PATH = ROOT / "scripts" / "audit_eucast_workbook.py"
SPEC = importlib.util.spec_from_file_location("audit_eucast_workbook", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
WORKBOOK = ROOT / "data" / "reference" / "EUCAST_v16.1_Breakpoint_Tables.xlsx"
MANIFEST = ROOT / "data" / "reference" / "EUCAST_v16.1_Breakpoint_Tables.manifest.json"


class EucastWorkbookAuditTests(unittest.TestCase):
    def test_rich_text_preserves_superscript_runs(self) -> None:
        node = ET.fromstring(
            f"""
            <si xmlns="{MODULE.MAIN_NS}">
              <r><t>(2)</t></r>
              <r><rPr><vertAlign val="superscript"/><b/></rPr><t>1</t></r>
            </si>
            """
        )
        parsed = MODULE._rich_text(node)
        self.assertEqual(parsed["text"], "(2)1")
        self.assertEqual(parsed["runs"][1]["vertical_alignment"], "superscript")
        self.assertTrue(parsed["runs"][1]["bold"])

    def test_range_intersection_uses_exact_cell_coordinates(self) -> None:
        self.assertTrue(MODULE._range_contains("I97", "I95:I101"))
        self.assertFalse(MODULE._range_contains("H97", "I95:I101"))
        with self.assertRaises(MODULE.EucastAuditError):
            MODULE._cell_coordinates("A0")

    @unittest.skipUnless(WORKBOOK.is_file() and MANIFEST.is_file(), "frozen EUCAST source absent")
    def test_frozen_workbook_contract_preserves_atu_and_note_markers(self) -> None:
        payload = MODULE.audit_workbook(WORKBOOK, MANIFEST)
        self.assertEqual(payload["source"]["sha256"], "ee5709eb5b7c9beb1c9ba8182a4e7032e71963aaca30e618b9f76dc53f6bab73")
        self.assertEqual(payload["findings"]["ciprofloxacin_non_meningitis_has_mic_atu"], "0.5")
        self.assertEqual(payload["findings"]["ciprofloxacin_non_meningitis_has_zone_atu"], "22-24")
        cells = {row["reference"]: row for row in payload["cells"]}
        self.assertEqual(cells["A83"]["rich_text_runs"][-1]["vertical_alignment"], "superscript")
        self.assertEqual(cells["F83"]["rich_text_runs"][-1]["text"], "B")
        self.assertEqual(cells["B97"]["rich_text_runs"][-1]["vertical_alignment"], "superscript")
        self.assertIn("1/A.", cells["I95"]["value"])
        self.assertEqual(payload["approval_status"], "INDEPENDENT_REVIEW_REQUIRED")


if __name__ == "__main__":
    unittest.main()
