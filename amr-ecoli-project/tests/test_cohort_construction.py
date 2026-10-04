from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "scripts" / "construct_research_cohort.py"
SPEC = importlib.util.spec_from_file_location("construct_research_cohort", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def ast_row(record: str = "ast-1", isolate: str = "PDT1") -> dict[str, str]:
    row = {field: "reported" for field in MODULE.AST_REQUIRED}
    row.update(
        {
            "source_dataset": "NCBI Pathogen Detection AST",
            "source_release": "snapshot",
            "source_fingerprint": "a" * 64,
            "source_row_number": "2",
            "source_ast_record_id": record,
            "isolate_id": isolate,
            "biosample_accession": "SAMN00000001",
            "assembly_accession": "GCA_000000001.1",
            "scientific_name": "Escherichia coli",
            "antibiotic": "ciprofloxacin",
            "submitted_ast_category": "R",
            "ast_measurement_type": "MIC",
            "measurement_sign": "=",
            "raw_ast_value": "0.25",
            "ast_value": "0.25",
            "ast_unit": "mg/L",
            "ast_method": "broth microdilution",
            "ast_testing_date": "2025-01-02",
            "collection_date": "2025-01-01",
            "country": "USA",
            "bioproject_accession": "PRJNA1",
            "clinical_indication": "non_meningitis",
            "ast_category": "S",
            "breakpoint_rule_id": "rule-1",
            "breakpoint_rule_row_sha256": "b" * 64,
            "breakpoint_artifact_sha256": "c" * 64,
            "breakpoint_standard": "EUCAST",
            "breakpoint_version": "16.1",
        }
    )
    return row


def metadata_row(isolate: str = "PDT1") -> dict[str, str]:
    return {
        "isolate_id": isolate,
        "biosample_accession": "SAMN00000001",
        "assembly_accession": "GCA_000000001.1",
        "species_method": "ANI",
        "species_result": "Escherichia coli",
        "genome_qc_status": "PASS",
        "site": "site-a",
        "surveillance_network": "network-a",
        "specimen_source": "blood",
        "deduplication_group": "patient-1",
        "lineage_group": "ST131",
        "genomic_cluster": "cluster-1",
        "evaluation_split": "development",
        "intended_use_population": "human_clinical",
        "genome_source": "registered_assembly",
        "genome_source_accession": "GCA_000000001.1",
        "genome_sha256": "d" * 64,
        "clinical_indication": "non_meningitis",
    }


def genome_row(eligible: str = "true") -> dict[str, str]:
    return {
        "isolate_id": "PDT1",
        "biosample_accession": "SAMN00000001",
        "assembly_accession": "GCA_000000001.1",
        "eligible_for_modeling": eligible,
        "overall_status": "PASS" if eligible == "true" else "EXCLUDE",
        "exclusion_reasons": "" if eligible == "true" else "ASSEMBLY_LENGTH_BELOW_MIN",
    }


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str] | None = None) -> None:
    fieldnames = fields or list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


class CohortConstructionTests(unittest.TestCase):
    def test_output_phenotype_and_rule_hash_come_from_interpreted_ast(self) -> None:
        accepted, excluded = MODULE.construct_rows([ast_row()], [metadata_row()])
        self.assertEqual(excluded, [])
        self.assertEqual(accepted[0]["ast_category"], "S")
        self.assertEqual(accepted[0]["submitted_ast_category"], "R")
        self.assertEqual(accepted[0]["breakpoint_rule_sha256"], "b" * 64)
        self.assertEqual(accepted[0]["species_method"], "ANI")

    def test_identity_mismatch_is_excluded(self) -> None:
        metadata = metadata_row()
        metadata["assembly_accession"] = "GCA_000000002.1"
        accepted, excluded = MODULE.construct_rows([ast_row()], [metadata])
        self.assertEqual(accepted, [])
        self.assertEqual(excluded[0]["cohort_exclusion_reason"], "identity_mismatch")

    def test_every_duplicate_isolate_drug_measurement_is_excluded(self) -> None:
        accepted, excluded = MODULE.construct_rows(
            [ast_row("ast-1"), ast_row("ast-2")], [metadata_row()]
        )
        self.assertEqual(accepted, [])
        self.assertEqual(len(excluded), 2)
        self.assertEqual(
            {row["cohort_exclusion_reason"] for row in excluded},
            {"duplicate_isolate_antibiotic_unreconciled"},
        )

    def test_metadata_must_be_one_row_per_isolate(self) -> None:
        with self.assertRaisesRegex(MODULE.CohortConstructionError, "not unique"):
            MODULE.construct_rows([ast_row()], [metadata_row(), metadata_row()])

    def test_executed_genome_qc_controls_modeling_eligibility(self) -> None:
        accepted, excluded = MODULE.construct_rows(
            [ast_row()], [metadata_row()], [genome_row("false")]
        )
        self.assertEqual(accepted, [])
        self.assertEqual(excluded[0]["cohort_exclusion_reason"], "genome_not_eligible")

        accepted, excluded = MODULE.construct_rows(
            [ast_row()], [metadata_row()], [genome_row("true")]
        )
        self.assertEqual(excluded, [])
        self.assertEqual(accepted[0]["genome_qc_status"], "PASS")

    def test_end_to_end_manifest_links_both_source_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ast = root / "ast.csv"
            metadata = root / "metadata.csv"
            output = root / "cohort.csv"
            exclusions = root / "exclusions.csv"
            manifest = root / "manifest.json"
            write_csv(ast, [ast_row()], list(MODULE.AST_REQUIRED))
            write_csv(metadata, [metadata_row()], list(MODULE.METADATA_REQUIRED))
            payload = MODULE.run_construction(ast, metadata, output, exclusions, manifest)
            self.assertEqual(payload["outputs"]["cohort"]["rows"], 1)
            self.assertEqual(payload["sources"]["interpreted_ast"]["sha256"], MODULE.sha256_file(ast))
            persisted = json.loads(manifest.read_text(encoding="utf-8"))
            self.assertEqual(persisted["scientific_status"], "SOURCE_LINKED_COHORT_CONSTRUCTION")


if __name__ == "__main__":
    unittest.main()


def mic_row(record: str, isolate: str, drug: str, sign: str, value: str, row_number: int,
            measurement: str = "MIC", unit: str = "mg/L") -> dict[str, str]:
    """A normalized-AST row carrying only what the MIC endpoint needs."""
    return {
        "source_dataset": "NCBI Pathogen Detection AST",
        "source_release": "snapshot",
        "source_fingerprint": "a" * 64,
        "source_row_number": str(row_number),
        "source_ast_record_id": record,
        "isolate_id": isolate,
        "biosample_accession": f"SAMN{isolate[-1]}0000001",
        "assembly_accession": f"GCA_00000000{isolate[-1]}.1",
        "scientific_name": "Escherichia coli",
        "antibiotic": drug,
        "ast_measurement_type": measurement,
        "measurement_sign": sign,
        "raw_ast_value": f"{sign}{value}",
        "ast_value": value,
        "ast_unit": unit,
        "collection_date": "2020-01-01",
        "country": "USA",
        "bioproject_accession": "PRJNA1",
    }


class RepeatReconciliationTests(unittest.TestCase):
    """Amendment 005 rule R, one test per case in the policy."""

    def _reconcile(self, specs):
        rows = [mic_row(f"r{i}", "PDT1", "ciprofloxacin", s, v, i + 2) for i, (s, v) in enumerate(specs)]
        return MODULE.reconcile_mic_repeats(rows)

    def test_open_and_closed_lower_bounds_intersect(self) -> None:
        row, reason = self._reconcile([(">", "32.0"), (">=", "64.0")])
        self.assertEqual(reason, "")
        self.assertEqual((row["measurement_sign"], row["ast_value"]), (">=", "64.0"))

    def test_upper_bounds_take_the_tightest(self) -> None:
        row, _ = self._reconcile([("<=", "4.0"), ("<=", "1.0"), ("<=", "1.0")])
        self.assertEqual((row["measurement_sign"], row["ast_value"]), ("<=", "1.0"))

    def test_equal_open_bound_beats_closed(self) -> None:
        row, _ = self._reconcile([(">", "2.0"), (">=", "2.0")])
        self.assertEqual((row["measurement_sign"], row["ast_value"]), (">", "2.0"))

    def test_matching_exact_values_reconcile(self) -> None:
        row, _ = self._reconcile([("=", "2.0"), ("=", "2.0"), ("<=", "4.0")])
        self.assertEqual((row["measurement_sign"], row["ast_value"]), ("=", "2.0"))

    def test_conflicting_exact_values_are_excluded(self) -> None:
        row, reason = self._reconcile([("=", "2.0"), ("=", "8.0")])
        self.assertIsNone(row)
        self.assertEqual(reason, "discordant_repeat_measurements")

    def test_disjoint_censored_intervals_are_excluded(self) -> None:
        _, reason = self._reconcile([("<=", "1.0"), (">", "4.0")])
        self.assertEqual(reason, "discordant_repeat_measurements")

    def test_touching_open_bounds_are_excluded(self) -> None:
        _, reason = self._reconcile([("<", "2.0"), (">=", "2.0")])
        self.assertEqual(reason, "discordant_repeat_measurements")

    def test_two_sided_interval_is_not_representable(self) -> None:
        _, reason = self._reconcile([("<=", "4.0"), (">=", "2.0")])
        self.assertEqual(reason, "repeat_not_representable")

    def test_unit_conflict_is_excluded(self) -> None:
        rows = [
            mic_row("r1", "PDT1", "ciprofloxacin", "=", "2.0", 2),
            mic_row("r2", "PDT1", "ciprofloxacin", "=", "2.0", 3, unit="ug/mL"),
        ]
        self.assertEqual(MODULE.reconcile_mic_repeats(rows), (None, "repeat_unit_conflict"))

    def test_provenance_keeps_every_contributing_record(self) -> None:
        row, _ = self._reconcile([(">=", "64.0"), (">", "32.0")])
        self.assertEqual(row["reconciled_from"], "r0;r1")
        self.assertEqual(row["reconciliation"], "interval_intersection")
        self.assertEqual(row["source_row_number"], "2")  # lowest source row carries it
        self.assertEqual(row["raw_ast_value"], ">=64.0;>32.0")

    def test_zone_rows_leave_the_mic_endpoint_before_reconciliation(self) -> None:
        rows = [
            mic_row("m", "PDT1", "gentamicin", "<=", "1.0", 2),
            mic_row("z", "PDT1", "gentamicin", "=", "20.0", 3, measurement="ZONE", unit="mm"),
        ]
        kept, excluded = MODULE.prepare_mic_rows(rows)
        self.assertEqual([r["source_ast_record_id"] for r in kept], ["m"])
        self.assertEqual(kept[0]["reconciliation"], "single")
        self.assertEqual(excluded[0]["cohort_exclusion_reason"], "not_mic_endpoint_measurement")

    def test_the_frozen_snapshot_repeat_groups_all_reconcile(self) -> None:
        """The six real MIC repeat groups from the development snapshot."""
        groups = {
            "ceftriaxone": [(">=", "64.0"), (">", "32.0")],
            "ciprofloxacin": [(">=", "4.0"), (">", "2.0"), (">", "2.0")],
            "gentamicin": [("<=", "4.0"), ("<=", "1.0"), ("<=", "1.0")],
        }
        expected = {"ceftriaxone": (">=", "64.0"), "ciprofloxacin": (">=", "4.0"), "gentamicin": ("<=", "1.0")}
        for isolate in ("PDT000041778.1", "PDT000041780.1"):
            for drug, specs in groups.items():
                rows = [mic_row(f"{isolate}-{drug}-{i}", isolate, drug, s, v, i + 2) for i, (s, v) in enumerate(specs)]
                row, reason = MODULE.reconcile_mic_repeats(rows)
                self.assertEqual(reason, "", (isolate, drug))
                self.assertEqual((row["measurement_sign"], row["ast_value"]), expected[drug])


class MicEndpointConstructionTests(unittest.TestCase):
    """End-to-end MIC build: would have caught the site/breakpoint KeyErrors."""

    def test_mic_profile_builds_without_categorical_fields(self) -> None:
        ast = [
            mic_row("a1", "PDT1", "ciprofloxacin", ">", "32.0", 2),
            mic_row("a2", "PDT1", "ciprofloxacin", ">=", "64.0", 3),
            mic_row("a3", "PDT2", "ciprofloxacin", "<=", "0.25", 4),
        ]
        metadata = []
        for isolate in ("PDT1", "PDT2"):
            metadata.append({
                "isolate_id": isolate,
                "biosample_accession": f"SAMN{isolate[-1]}0000001",
                "assembly_accession": f"GCA_00000000{isolate[-1]}.1",
                "species_method": "ANI", "species_result": "Escherichia coli",
                "genome_qc_status": "PASS", "specimen_source": "urine",
                "deduplication_group": f"dedup-{isolate}", "lineage_group": "ST131",
                "genomic_cluster": f"cluster-{isolate}", "evaluation_split": "development",
                "intended_use_population": "human_clinical",
                "genome_source": "registered_assembly",
                "genome_source_accession": f"GCA_00000000{isolate[-1]}.1",
                "genome_sha256": "d" * 64,
            })
        accepted, excluded = MODULE.construct_rows(ast, metadata, profile="mic_regression")
        self.assertEqual(excluded, [])
        self.assertEqual(len(accepted), 2)
        by_isolate = {row["isolate_id"]: row for row in accepted}
        self.assertEqual(by_isolate["PDT1"]["measurement_sign"], ">=")
        self.assertEqual(by_isolate["PDT1"]["reconciled_from"], "a1;a2")
        self.assertNotIn("site", by_isolate["PDT1"])
        self.assertNotIn("breakpoint_rule_sha256", by_isolate["PDT1"])

    def test_categorical_profile_keeps_the_blanket_duplicate_exclusion(self) -> None:
        """Rule R is scoped to the MIC endpoint; the blocked secondary is unchanged."""
        rows = [ast_row("ast-1"), ast_row("ast-2")]
        rows[1]["source_row_number"] = "3"
        accepted, excluded = MODULE.construct_rows(rows, [metadata_row()])
        self.assertEqual(accepted, [])
        self.assertTrue(all(r["cohort_exclusion_reason"] == "duplicate_isolate_antibiotic_unreconciled" for r in excluded))


class MultiSourceConstructionTests(unittest.TestCase):
    def test_development_and_external_sources_combine_with_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dev = [mic_row("d1", "PDT1", "ciprofloxacin", "<=", "0.25", 2)]
            ext = [mic_row("JARBS:JB2:ciprofloxacin", "PDT2", "ciprofloxacin", ">", "2", 5)]
            ext[0]["source_dataset"] = "JARBS-GNR Supplementary Data 6"
            write_csv(root / "dev.csv", dev)
            write_csv(root / "ext.csv", ext)
            metadata = []
            for isolate, split in (("PDT1", "development"), ("PDT2", "external")):
                metadata.append({
                    "isolate_id": isolate, "biosample_accession": f"SAMN{isolate[-1]}0000001",
                    "assembly_accession": f"GCA_00000000{isolate[-1]}.1", "species_method": "ANI",
                    "species_result": "Escherichia coli", "genome_qc_status": "PASS",
                    "specimen_source": "urine", "deduplication_group": f"dedup-{isolate}",
                    "lineage_group": "ST131", "genomic_cluster": f"c-{isolate}", "evaluation_split": split,
                    "intended_use_population": "human_clinical", "genome_source": "registered_assembly",
                    "genome_source_accession": f"GCA_00000000{isolate[-1]}.1", "genome_sha256": "d" * 64,
                })
            write_csv(root / "meta.csv", metadata)
            payload = MODULE.run_construction(
                [root / "dev.csv", root / "ext.csv"], root / "meta.csv", root / "cohort.csv",
                root / "exclusions.csv", root / "manifest.json", profile="mic_regression")
            self.assertEqual(payload["outputs"]["cohort"]["rows"], 2)
            self.assertEqual(len(payload["sources"]["ast_sources"]), 2)
            self.assertEqual(payload["sources"]["interpreted_ast"]["path"], str(root / "dev.csv"))

    def test_same_file_twice_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_csv(root / "dev.csv", [mic_row("d1", "PDT1", "ciprofloxacin", "=", "1", 2)])
            write_csv(root / "meta.csv", [{"isolate_id": "PDT1"}])
            with self.assertRaises(MODULE.CohortConstructionError):
                MODULE.run_construction([root / "dev.csv", root / "dev.csv"], root / "meta.csv",
                                        root / "c.csv", root / "e.csv", root / "m.json", profile="mic_regression")



class DoublingScaleTests(unittest.TestCase):
    def test_off_scale_values_are_excluded_and_labels_snap(self) -> None:
        rows = [mic_row("a", "PDT1", "ciprofloxacin", "=", "0.12", 2),
                mic_row("b", "PDT2", "ciprofloxacin", "=", "2.5", 3),
                mic_row("c", "PDT3", "ciprofloxacin", "<=", "0.2", 4),
                mic_row("d", "PDT4", "ciprofloxacin", ">", "200", 5)]
        kept, excluded = MODULE.prepare_mic_rows(rows)
        self.assertEqual([r["isolate_id"] for r in kept], ["PDT1"])
        self.assertEqual({r["cohort_exclusion_reason"] for r in excluded}, {"mic_value_off_doubling_scale"})

    def test_off_scale_row_does_not_join_a_repeat_group(self) -> None:
        rows = [mic_row("a", "PDT1", "ciprofloxacin", "=", "1", 2),
                mic_row("b", "PDT1", "ciprofloxacin", "=", "2.5", 3)]
        kept, excluded = MODULE.prepare_mic_rows(rows)
        self.assertEqual(len(kept), 1)
        self.assertEqual(kept[0]["reconciliation"], "single")

    def test_tolerance_matches_the_evaluator(self) -> None:
        from scripts import evaluate_mic_predictions
        self.assertEqual(MODULE.DOUBLING_SCALE_TOLERANCE, evaluate_mic_predictions.SNAP_TOLERANCE)
