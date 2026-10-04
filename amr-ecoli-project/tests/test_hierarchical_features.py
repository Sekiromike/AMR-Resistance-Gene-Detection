from __future__ import annotations

import unittest

import pandas as pd

from scripts import build_hierarchical_features as MODULE


class FamilyTests(unittest.TestCase):
    def test_family_rules(self) -> None:
        cases = {
            "blaCTX-M-15": "blaCTX-M", "blaCTX-M-2": "blaCTX-M", "blaCTX-M-27": "blaCTX-M",
            "blaTEM-1": "blaTEM", "blaIMP-1": "blaIMP", "blaNDM-5": "blaNDM", "blaOXA-1": "blaOXA",
            "blaCMY-2": "blaCMY", "blaEC-5": "blaEC", "blaEC": "blaEC", "blaCTX-M": "blaCTX-M",
            "qnrS1": "qnrS", "qnrB19": "qnrB", "aac(3)-IId": "aac(3)-II", "aac(3)-IIe": "aac(3)-II",
            "aac(3)-VIa": "aac(3)-VI", "ant(2'')-Ia": "ant(2'')-I", "sul1": "sul", "tet(A)": "tet(A)",
            "armA": "armA", "rmtB1": "rmtB",
        }
        for symbol, family in cases.items():
            self.assertEqual(MODULE.gene_family(symbol), family, symbol)
        self.assertIsNone(MODULE.gene_family("gyrA_S83L"))

    def test_hierarchy_rows_are_added_without_changing_allele_rows(self) -> None:
        features = pd.DataFrame([
            {"isolate_id": "A", "cohort": "external", "element_symbol": "blaCTX-M-2", "subtype": "AMR",
             "class": "BETA-LACTAM", "subclass": "CEPHALOSPORIN", "method": "EXACTX"},
            {"isolate_id": "A", "cohort": "external", "element_symbol": "blaTEM-1", "subtype": "AMR",
             "class": "BETA-LACTAM", "subclass": "BETA-LACTAM", "method": "EXACTX"},
            {"isolate_id": "A", "cohort": "external", "element_symbol": "gyrA_S83L", "subtype": "POINT",
             "class": "QUINOLONE", "subclass": "QUINOLONE", "method": "POINTX"},
            {"isolate_id": "B", "cohort": "development", "element_symbol": "blaCTX-M-15", "subtype": "AMR",
             "class": "BETA-LACTAM", "subclass": "CEPHALOSPORIN", "method": "EXACTX"},
        ])
        out = MODULE.add_hierarchy(features)
        a = set(out.loc[out.isolate_id == "A", "element_symbol"])
        self.assertTrue({"blaCTX-M-2", "blaTEM-1", "gyrA_S83L"} <= a)
        self.assertTrue({"family:blaCTX-M", "family:blaTEM", "subclass:CEPHALOSPORIN",
                         "subclass:BETA-LACTAM", "subclass:QUINOLONE"} <= a)
        self.assertNotIn("family:gyrA", {s.split("_")[0] for s in a})
        # An unseen allele shares its family and class with a common one.
        b = set(out.loc[out.isolate_id == "B", "element_symbol"])
        self.assertIn("family:blaCTX-M", b & a)
        self.assertEqual(len(out), len(out.drop_duplicates(["isolate_id", "element_symbol"])))
