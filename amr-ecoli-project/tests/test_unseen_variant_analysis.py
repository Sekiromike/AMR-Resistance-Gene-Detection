from __future__ import annotations

import unittest

import pandas as pd

from scripts import unseen_variant_analysis as MODULE


class UnseenVariantTests(unittest.TestCase):
    def test_share_on_resistant_side_by_carrier_group(self) -> None:
        scored = pd.DataFrame({
            "isolate_id": ["U1", "U2", "C1", "C2", "S1"], "antibiotic": "ceftriaxone",
            "evaluation_split": "external", "lineage_group": ["L1", "L2", "L3", "L4", "L5"],
            "measurement_sign": [">", ">", ">", ">", "<="], "ast_value": [64.0, 64.0, 64.0, 64.0, 0.5],
            "predicted_log2_mic": [3.0, -2.0, 5.0, 4.0, -3.0]})
        features = pd.DataFrame({"isolate_id": ["U1", "U2", "C1", "C2", "C2"],
                                 "element_symbol": ["blaCTX-M-2", "blaIMP-1", "blaCTX-M-15", "blaCTX-M-27", "blaTEM-1"]})
        out = MODULE.read_out(scored, features)
        self.assertEqual((out["unseen_n"], out["common_n"]), (2, 2))
        self.assertEqual(out["unseen_resistant_side"], 0.5)
        self.assertEqual(out["common_resistant_side"], 1.0)
