from __future__ import annotations

import csv
import hashlib
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from scripts import attach_external_references as ATTACH
from scripts import predict_external as MODULE
from tests_support_mic import small_cohort

SENTINEL = "SEALED-NOT-A-NUMBER"


def sealed_cohort(root: Path) -> tuple[Path, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Development from the shared fixture plus external rows whose MIC columns are unparseable."""
    cohort, splits, features = small_cohort()
    development = cohort.assign(evaluation_split="development")
    external = development.head(40).copy()
    external["isolate_id"] = "X" + external["isolate_id"]
    external["evaluation_split"] = "external"
    external["measurement_sign"] = SENTINEL
    external["ast_value"] = SENTINEL
    extra = features[features["isolate_id"].isin(development.head(40)["isolate_id"])].copy()
    extra["isolate_id"] = "X" + extra["isolate_id"]
    path = root / "cohort.csv"
    pd.concat([development, external]).to_csv(path, index=False)
    return path, splits, pd.concat([features, extra], ignore_index=True), external


class SealTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def test_external_rows_are_reduced_to_genotype_columns_on_read(self) -> None:
        path, *_ = sealed_cohort(self.root)
        development, external = MODULE.load_cohort_sealed(path)
        self.assertEqual(list(external.columns), list(MODULE.GENOTYPE_COLUMNS))
        self.assertNotIn(SENTINEL, development.to_numpy().ravel().tolist())
        self.assertEqual(len(external), 40)

    def test_all_models_predict_without_touching_external_references(self) -> None:
        # The sentinel would raise if any external MIC were parsed.
        path, splits, features, _ = sealed_cohort(self.root)
        development, external = MODULE.load_cohort_sealed(path)
        for model in ("constant", "ridge_censored:amr", "xgboost_aft:amr"):
            predictions, records = MODULE.predict(model, development, external, splits, features,
                                                  "human_clinical", seed=7)
            self.assertEqual(len(predictions), 40, model)
            self.assertTrue(np.isfinite(predictions["predicted_log2_mic"]).all(), model)
            self.assertNotIn("ast_value", predictions.columns)
            self.assertNotIn("measurement_sign", predictions.columns)

    def test_neighbor_uses_the_best_development_relative(self) -> None:
        path, splits, features, _ = sealed_cohort(self.root)
        development, external = MODULE.load_cohort_sealed(path)
        comparisons = self.root / "comparisons.tsv"
        lines = ["development_sequence_id\texternal_sequence_id\tani_percent\taligned_fraction_development\taligned_fraction_external"]
        for x in external["isolate_id"].unique():
            twin = x[1:]
            lines.append(f"development::{twin}\texternal::{x}\t99.95\t0.99\t0.99")
            lines.append(f"development::I599\texternal::{x}\t99.10\t0.99\t0.99")
        comparisons.write_text("\n".join(lines) + "\n", encoding="utf-8")
        predictions, _ = MODULE.predict("neighbor", development, external, splits, features, "human_clinical",
                                        seed=7, comparisons=comparisons)
        dev = development.set_index("isolate_id")
        for row in predictions.itertuples():
            twin = dev.loc[row.isolate_id[1:]]
            expected = MODULE.boundary_prediction(twin["measurement_sign"], float(twin["ast_value"]))
            self.assertEqual(row.predicted_log2_mic, expected)


class AttachTests(unittest.TestCase):
    def test_references_join_only_to_locked_predictions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cohort = root / "cohort.csv"
            pd.DataFrame({"isolate_id": ["E1", "E2", "D1"], "antibiotic": "gentamicin",
                          "evaluation_split": ["external", "external", "development"],
                          "measurement_sign": ["<=", ">", "="], "ast_value": ["1", "8", "2"]}).to_csv(cohort, index=False)
            predictions = root / "predictions_constant_human_clinical.csv"
            pd.DataFrame({"isolate_id": ["E1", "E2"], "antibiotic": "gentamicin", "evaluation_split": "external",
                          "lineage_group": ["SLV:ST1", "SLV:ST2"], "predicted_log2_mic": [0.0, 0.0]}).to_csv(predictions, index=False)
            lock = root / "lock.tsv"
            digest = hashlib.sha256(predictions.read_bytes()).hexdigest()
            lock.write_text(f"file\tsha256\n{predictions.name}\t{digest}\n", encoding="utf-8")
            merged = ATTACH.attach(predictions, cohort, lock)
            self.assertEqual(merged["ast_value"].tolist(), [1.0, 8.0])
            # Any change after locking is refused.
            predictions.write_text(predictions.read_text(encoding="utf-8").replace("0.0", "1.0"), encoding="utf-8")
            with self.assertRaisesRegex(ATTACH.LockError, "differs"):
                ATTACH.attach(predictions, cohort, lock)
            lock.write_text("file\tsha256\n", encoding="utf-8")
            with self.assertRaisesRegex(ATTACH.LockError, "not in the lock"):
                ATTACH.attach(predictions, cohort, lock)
