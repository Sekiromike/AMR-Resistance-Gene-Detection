from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).parents[1] / "scripts" / "build_split_manifest.py"
SPEC = importlib.util.spec_from_file_location("build_split_manifest", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def cohort() -> pd.DataFrame:
    rows = []
    for index in range(12):
        rows.append(
            {
                "isolate_id": f"iso-{index}",
                "antibiotic": "ciprofloxacin",
                "ast_category": "R" if index % 2 else "S",
                "evaluation_split": "external" if index >= 10 else "development",
                "genomic_cluster": f"cluster-{index}",
                "deduplication_group": f"dedup-{index}",
                "lineage_group": f"ST{index}",
            }
        )
    return pd.DataFrame(rows)


class SplitManifestTests(unittest.TestCase):
    def test_is_deterministic_and_external_has_no_fold(self) -> None:
        first = MODULE.build_manifest(cohort(), n_folds=5, seed=7)
        second = MODULE.build_manifest(cohort(), n_folds=5, seed=7)
        pd.testing.assert_frame_equal(first, second)
        self.assertTrue(first.loc[first["evaluation_split"].eq("external"), "cv_fold"].isna().all())
        self.assertEqual(
            first.loc[first["evaluation_split"].eq("development"), "cv_fold"].nunique(),
            5,
        )

    def test_cluster_crossing_external_boundary_fails(self) -> None:
        data = cohort()
        data.loc[10, "genomic_cluster"] = data.loc[0, "genomic_cluster"]
        with self.assertRaisesRegex(ValueError, "cross development/external"):
            MODULE.build_manifest(data, n_folds=5)

    def test_lineage_never_crosses_development_folds(self) -> None:
        data = cohort()
        data.loc[data["evaluation_split"].eq("development"), "lineage_group"] = [
            f"ST{index % 5}" for index in range(10)
        ]
        result = MODULE.build_manifest(data, n_folds=5, seed=3)
        development = result.loc[result["evaluation_split"].eq("development")]
        self.assertTrue((development.groupby("lineage_group")["cv_fold"].nunique() == 1).all())


if __name__ == "__main__":
    unittest.main()


class GroupingVariableTests(unittest.TestCase):
    """Amendment 004: fold grouping is by lineage, not by ANI components."""

    def test_chained_genomic_cluster_stops_the_split_instead_of_collapsing_folds(self) -> None:
        """Amendment 007: clusters join the CV unit, so chaining must fail loudly."""
        data = cohort()
        development = data["evaluation_split"].eq("development")
        # One chained component spanning every development isolate, as observed
        # in the external set where the largest component held 42% of the cohort.
        data.loc[development, "genomic_cluster"] = "chained-component"
        with self.assertRaisesRegex(ValueError, "chained"):
            MODULE.build_manifest(data, n_folds=5, seed=3)

    def test_a_large_lineage_alone_is_not_mistaken_for_chaining(self) -> None:
        data = cohort()
        development = data["evaluation_split"].eq("development")
        data.loc[development & data["isolate_id"].isin([f"iso-{i}" for i in range(4)]), "lineage_group"] = "ST131"
        result = MODULE.build_manifest(data, n_folds=5, seed=3).set_index("isolate_id")
        self.assertEqual(result.loc[[f"iso-{i}" for i in range(4)], "cv_fold"].nunique(), 1)

    def test_near_neighbours_in_different_lineages_share_a_fold(self) -> None:
        data = cohort()
        data.loc[data["isolate_id"].isin(["iso-1", "iso-2"]), "genomic_cluster"] = "shared-cluster"
        result = MODULE.build_manifest(data, n_folds=5, seed=3).set_index("isolate_id")
        self.assertEqual(result.loc["iso-1", "cv_fold"], result.loc["iso-2", "cv_fold"])
        self.assertEqual(result.loc["iso-1", "cv_group"], result.loc["iso-2", "cv_group"])

    def test_deduplication_group_still_binds_isolates_together(self) -> None:
        data = cohort()
        development = data.index[data["evaluation_split"].eq("development")]
        # Two isolates with different lineages but the same duplicate group.
        data.loc[development[0], "deduplication_group"] = "shared-dedup"
        data.loc[development[1], "deduplication_group"] = "shared-dedup"
        groups = MODULE.development_cv_groups(data)
        first = data.loc[development[0], "isolate_id"]
        second = data.loc[development[1], "isolate_id"]
        self.assertEqual(groups[first], groups[second])

    def test_missing_grouping_column_is_reported_clearly(self) -> None:
        data = cohort().drop(columns=["deduplication_group"])
        with self.assertRaisesRegex(ValueError, "deduplication_group"):
            MODULE.development_cv_groups(data)


def mic_cohort() -> pd.DataFrame:
    """MIC-endpoint cohort: no ast_category, comparators and populations instead."""
    signs = ["<=", "=", ">", "<=", "=", ">", "<=", "=", ">", "<="]
    rows = []
    for index in range(40):
        rows.append({
            "isolate_id": f"iso-{index}",
            "antibiotic": "ciprofloxacin" if index % 4 else "gentamicin",
            "measurement_sign": signs[index % len(signs)],
            "intended_use_population": "human_clinical" if index % 3 else "non_human_or_environmental",
            "evaluation_split": "external" if index >= 34 else "development",
            "genomic_cluster": f"cluster-{index}",
            "deduplication_group": f"dedup-{index}",
            "lineage_group": f"ST{index % 17}",
        })
    return pd.DataFrame(rows)


class MicEndpointSplitTests(unittest.TestCase):
    def test_mic_cohort_without_ast_category_is_split_by_lineage(self) -> None:
        frame = mic_cohort()
        table = MODULE.build_manifest(frame, n_folds=5, seed=20260815)
        development = table[table["evaluation_split"].eq("development")]
        self.assertTrue(development["cv_fold"].notna().all())
        self.assertTrue(table.loc[table["evaluation_split"].eq("external"), "cv_fold"].isna().all())
        self.assertTrue((development.groupby("lineage_group")["cv_fold"].nunique() == 1).all())

    def test_mic_stratum_is_censoring_class_within_population(self) -> None:
        stratum = MODULE.balance_stratum(mic_cohort().head(3))
        self.assertEqual(stratum.tolist(), ["non_human_or_environmental|left", "human_clinical|exact",
                                            "human_clinical|right"])
        bad = mic_cohort().head(1).assign(measurement_sign="~")
        with self.assertRaises(ValueError):
            MODULE.balance_stratum(bad)
        with self.assertRaises(ValueError):
            MODULE.balance_stratum(mic_cohort().drop(columns="measurement_sign"))

    def test_manifest_never_reports_external_strata(self) -> None:
        import json
        import subprocess
        import sys
        import tempfile
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mic_cohort().to_csv(root / "cohort.csv", index=False)
            subprocess.run([sys.executable, str(MODULE_PATH), "--cohort", str(root / "cohort.csv"),
                            "--output", str(root / "splits.csv"), "--manifest", str(root / "m.json")],
                           check=True, capture_output=True)
            manifest = json.loads((root / "m.json").read_text(encoding="utf-8"))
        self.assertTrue(all(record["cv_fold"] is not None for record in manifest["distribution"]))
        self.assertEqual(sum(manifest["external_rows_by_antibiotic"].values()), 6)
        self.assertNotIn("external", json.dumps(manifest["distribution"]))
