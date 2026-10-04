from __future__ import annotations

import importlib.util
import csv
import hashlib
import tempfile
import unittest
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).parents[1] / "scripts" / "validate_research_cohort.py"
SPEC = importlib.util.spec_from_file_location("validate_research_cohort", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def valid_cohort() -> pd.DataFrame:
    rows = []
    for i, (category, split) in enumerate(
        [("S", "development"), ("R", "development"), ("S", "external"), ("R", "external")],
        start=1,
    ):
        rows.append(
            {
                "isolate_id": f"isolate-{i}",
                "biosample_accession": f"SAMN{i:08d}",
                "assembly_accession": f"GCA_{i:09d}.1",
                "scientific_name": "Escherichia coli",
                "species_method": "ANI",
                "species_result": "Escherichia coli",
                "genome_qc_status": "PASS",
                "antibiotic": "ciprofloxacin",
                "source_dataset": "curated-test",
                "source_release": "2026-08-15",
                "source_fingerprint": "a" * 64,
                "source_row_number": i,
                "source_ast_record_id": f"ast-{i}",
                "ast_measurement_type": "MIC",
                "measurement_sign": "=",
                "raw_ast_value": "0.25" if category == "S" else "4",
                "ast_value": 0.25 if category == "S" else 4.0,
                "ast_unit": "mg/L",
                "ast_method": "broth microdilution",
                "ast_testing_date": f"202{min(i, 4)}-01-02",
                "submitted_ast_category": category,
                "ast_category": category,
                "breakpoint_standard": "CLSI",
                "breakpoint_version": "M100-Ed36",
                "breakpoint_rule_id": "test-ciprofloxacin-mic",
                "breakpoint_rule_sha256": "b" * 64,
                "breakpoint_artifact_sha256": "c" * 64,
                "collection_date": f"202{min(i, 4)}-01-01",
                "site": "site-a" if split == "development" else "site-b",
                "country": "USA",
                "bioproject_accession": "PRJNA000001" if split == "development" else "PRJNA000002",
                "clinical_indication": "non_meningitis",
                "surveillance_network": "network-a" if split == "development" else "network-b",
                "specimen_source": "blood",
                "deduplication_group": f"dedup-{i}",
                "lineage_group": f"lineage-{i}",
                "genomic_cluster": f"cluster-{i}",
                "evaluation_split": split,
                "intended_use_population": "human_clinical",
                "genome_source": "registered_assembly",
                "genome_source_accession": f"GCA_{i:09d}.1",
                "genome_sha256": "d" * 64,
            }
        )
    return pd.DataFrame(rows)


class ResearchCohortTests(unittest.TestCase):
    def test_empty_cohort_fails(self) -> None:
        report = MODULE.validate_cohort(valid_cohort().iloc[0:0], min_per_category=1)
        self.assertEqual(report["status"], "FAIL")
        self.assertIn("Cohort table has no rows.", report["failures"])

    def test_valid_curated_cohort_passes(self) -> None:
        report = MODULE.validate_cohort(valid_cohort(), min_per_category=1)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["n_isolates"], 4)

    def test_mixed_breakpoint_versions_fail(self) -> None:
        cohort = valid_cohort()
        cohort.loc[0, "breakpoint_version"] = "M100-Ed35"
        report = MODULE.validate_cohort(cohort, min_per_category=1)
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("Mixed breakpoint" in failure for failure in report["failures"]))

    def test_missing_ast_provenance_fails(self) -> None:
        cohort = valid_cohort().drop(columns=["ast_method"])
        report = MODULE.validate_cohort(cohort, min_per_category=1)
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("Missing required columns" in failure for failure in report["failures"]))

    def test_near_neighbor_cluster_crossing_external_boundary_fails(self) -> None:
        cohort = valid_cohort()
        cohort.loc[2, "genomic_cluster"] = cohort.loc[0, "genomic_cluster"]
        report = MODULE.validate_cohort(cohort, min_per_category=1)
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("Near-neighbor" in failure for failure in report["failures"]))

    def test_measurement_units_must_match_type(self) -> None:
        cohort = valid_cohort()
        cohort.loc[0, "ast_unit"] = "mm"
        report = MODULE.validate_cohort(cohort, min_per_category=1)
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("type/unit" in failure for failure in report["failures"]))

    def test_external_requires_both_classes(self) -> None:
        cohort = valid_cohort()
        cohort.loc[cohort["evaluation_split"].eq("external"), "ast_category"] = "S"
        report = MODULE.validate_cohort(cohort, min_per_category=1)
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("both S and R" in failure for failure in report["failures"]))

    def test_frozen_rules_reject_an_inverted_cohort_label(self) -> None:
        cohort = valid_cohort().iloc[[0]].copy()
        cohort["breakpoint_standard"] = "EUCAST"
        cohort["breakpoint_version"] = "16.1"
        cohort["breakpoint_rule_id"] = "EUCAST-16.1-CIP-MIC"
        cohort["ast_category"] = "R"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifact = root / "breakpoints.xlsx"
            artifact.write_bytes(b"test-only frozen artifact\n")
            artifact_sha = hashlib.sha256(artifact.read_bytes()).hexdigest()
            rules = root / "rules.csv"
            fields = [
                "breakpoint_rule_id",
                "breakpoint_artifact_sha256",
                "breakpoint_standard",
                "breakpoint_version",
                "antibiotic",
                "ast_measurement_type",
                "ast_unit",
                "susceptible_breakpoint",
                "resistant_breakpoint",
                "ast_method",
                "disk_potency",
                "site",
            ]
            row = {
                "breakpoint_rule_id": "EUCAST-16.1-CIP-MIC",
                "breakpoint_artifact_sha256": artifact_sha,
                "breakpoint_standard": "EUCAST",
                "breakpoint_version": "16.1",
                "antibiotic": "ciprofloxacin",
                "ast_measurement_type": "MIC",
                "ast_unit": "mg/L",
                "susceptible_breakpoint": "0.25",
                "resistant_breakpoint": "0.5",
                "ast_method": "broth microdilution",
                "disk_potency": "",
                "site": "",
            }
            with rules.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerow(row)
            loaded = MODULE.load_breakpoint_rules(rules, artifact)[0]
            cohort["breakpoint_rule_sha256"] = loaded.row_sha256
            cohort["breakpoint_artifact_sha256"] = artifact_sha
            failures = MODULE.validate_breakpoint_derivation(cohort, rules, artifact)
            self.assertTrue(any("disagree" in failure for failure in failures))


if __name__ == "__main__":
    unittest.main()


def mic_cohort() -> pd.DataFrame:
    """A cohort carrying only what the interval-censored MIC endpoint needs.

    It deliberately omits assay method, AST testing date, clinical indication,
    laboratory site, surveillance network, and every breakpoint field, because
    docs/ENDPOINT_AMENDMENT.md records that no public source supplies them.
    """
    rows = []
    for i, (value, split) in enumerate(
        [(0.25, "development"), (4.0, "development"), (0.25, "external"), (4.0, "external")],
        start=1,
    ):
        rows.append(
            {
                "isolate_id": f"isolate-{i}",
                "biosample_accession": f"SAMN{i:08d}",
                "assembly_accession": f"GCA_{i:09d}.1",
                "scientific_name": "Escherichia coli",
                "species_method": "ANI",
                "species_result": "Escherichia coli",
                "genome_qc_status": "PASS",
                "antibiotic": "ciprofloxacin",
                "source_dataset": "curated-test",
                "source_release": "2026-09-20",
                "source_fingerprint": "a" * 64,
                "source_row_number": i,
                "source_ast_record_id": f"ast-{i}",
                "ast_measurement_type": "MIC",
                "measurement_sign": "=",
                "raw_ast_value": str(value),
                "ast_value": value,
                "ast_unit": "mg/L",
                "collection_date": f"202{min(i, 4)}-01-01",
                "country": "USA",
                "bioproject_accession": "PRJNA000001" if split == "development" else "PRJNA000002",
                "specimen_source": "blood",
                "deduplication_group": f"dedup-{i}",
                "lineage_group": f"lineage-{i}",
                "genomic_cluster": f"cluster-{i}",
                "evaluation_split": split,
                "intended_use_population": "human_clinical",
                "genome_source": "registered_assembly",
                "genome_source_accession": f"GCA_{i:09d}.1",
                "genome_sha256": "d" * 64,
            }
        )
    return pd.DataFrame(rows)


class EndpointProfileTests(unittest.TestCase):
    """Guards for the two-profile contract in docs/ENDPOINT_AMENDMENT.md."""

    def test_mic_profile_accepts_cohort_without_interpretation_metadata(self) -> None:
        report = MODULE.validate_cohort(
            mic_cohort(), min_per_category=1, profile="mic_regression"
        )
        self.assertEqual(report["status"], "PASS", report["failures"])
        self.assertEqual(report["endpoint_profile"], "mic_regression")
        self.assertEqual(report["n_isolates"], 4)

    def test_categorical_profile_still_fails_closed_on_the_same_cohort(self) -> None:
        """The strict contract must not be weakened by the new profile."""
        report = MODULE.validate_cohort(
            mic_cohort(), min_per_category=1, profile="categorical_sir"
        )
        self.assertEqual(report["status"], "FAIL")
        missing = [f for f in report["failures"] if "Missing required columns" in f]
        self.assertTrue(missing, report["failures"])
        for field in (
            "ast_method",
            "ast_testing_date",
            "clinical_indication",
            "site",
            "breakpoint_version",
            "ast_category",
        ):
            self.assertIn(field, missing[0])

    def test_default_profile_is_the_strict_contract(self) -> None:
        report = MODULE.validate_cohort(mic_cohort(), min_per_category=1)
        self.assertEqual(report["status"], "FAIL")
        self.assertEqual(report["endpoint_profile"], "categorical_sir")

    def test_unknown_profile_is_rejected(self) -> None:
        report = MODULE.validate_cohort(
            mic_cohort(), min_per_category=1, profile="anything_goes"
        )
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("Unknown endpoint profile" in f for f in report["failures"]))

    def test_absent_site_separation_is_warned_not_silently_passed(self) -> None:
        report = MODULE.validate_cohort(
            mic_cohort(), min_per_category=1, profile="mic_regression"
        )
        self.assertTrue(
            any("separation could not be checked" in w for w in report["warnings"]),
            report["warnings"],
        )

    def test_mic_profile_still_enforces_near_neighbor_separation(self) -> None:
        cohort = mic_cohort()
        cohort.loc[2, "genomic_cluster"] = cohort.loc[0, "genomic_cluster"]
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("Near-neighbor" in f for f in report["failures"]))

    def test_mic_profile_still_enforces_deduplication_separation(self) -> None:
        cohort = mic_cohort()
        cohort.loc[2, "deduplication_group"] = cohort.loc[0, "deduplication_group"]
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("deduplication" in f for f in report["failures"]))

    def test_mic_profile_rejects_invalid_comparator(self) -> None:
        cohort = mic_cohort()
        cohort.loc[0, "measurement_sign"] = "~"
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("comparator" in f for f in report["failures"]))

    def test_fully_censored_partition_fails(self) -> None:
        cohort = mic_cohort()
        cohort["measurement_sign"] = "<="
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("no exact MIC measurement" in f for f in report["failures"]))

    def test_censoring_profile_is_reported(self) -> None:
        cohort = mic_cohort()
        cohort.loc[1, "measurement_sign"] = ">"
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertIn("censoring", report)
        development = [
            row for row in report["censoring"] if row["evaluation_split"] == "development"
        ]
        self.assertEqual(len(development), 1)
        self.assertEqual(development[0]["n_censored"], 1)
        self.assertEqual(development[0]["fraction_censored"], 0.5)

    def test_external_censoring_is_checked_but_never_reported(self) -> None:
        cohort = mic_cohort()
        cohort.loc[3, "measurement_sign"] = ">"
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertEqual(report["status"], "PASS")
        external = [row for row in report["censoring"] if row["evaluation_split"] == "external"]
        self.assertEqual(len(external), 1)
        self.assertEqual(external[0]["n"], 2)
        self.assertEqual(external[0]["n_censored"], "sealed")
        self.assertEqual(external[0]["fraction_censored"], "sealed")
        # The fully-censored structural check still runs on the sealed split.
        cohort.loc[2, "measurement_sign"] = "<="
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertEqual(report["status"], "FAIL")

    def test_partial_insdc_dates_are_valid_without_filling_month_or_day(self) -> None:
        today = MODULE.date(2026, 9, 27)
        for value in ("2015", "2015-06", "2015-06-30", "2026-09-27"):
            self.assertTrue(MODULE.valid_partial_date(value, today), value)
        for value in ("2015-13", "2015-02-30", "1899", "2027", "2026-10", "15-06-2015", "2015/2016", "missing"):
            self.assertFalse(MODULE.valid_partial_date(value, today), value)
        cohort = mic_cohort()
        cohort["collection_date"] = ["2015", "2015-06", "2015-06-30", "2019"]
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertEqual(report["status"], "PASS", report["failures"])
        cohort.loc[0, "collection_date"] = "2015-02-30"
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertTrue(any("Invalid collection dates: 1" in f for f in report["failures"]))

    def test_absent_country_and_date_are_warnings_under_mic_but_failures_under_categorical(self) -> None:
        cohort = mic_cohort()
        cohort.loc[0, "country"] = ""
        cohort.loc[1, "collection_date"] = ""
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertEqual(report["status"], "PASS", report["failures"])
        self.assertTrue(any("'collection_date': 1, 'country': 1" in w for w in report["warnings"]))
        cohort.loc[3, "specimen_source"] = ""
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertEqual(report["status"], "PASS", report["failures"])
        self.assertTrue(any("'specimen_source': 1" in w for w in report["warnings"]))
        cohort.loc[2, "bioproject_accession"] = ""
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertTrue(any("bioproject_accession" in f for f in report["failures"]))

    def test_profile_resolution_prefers_override_then_config_then_strict(self) -> None:
        self.assertEqual(
            MODULE.resolve_endpoint_profile({"primary_endpoint": "mic_regression"}),
            "mic_regression",
        )
        self.assertEqual(
            MODULE.resolve_endpoint_profile(
                {"primary_endpoint": "mic_regression"}, "categorical_sir"
            ),
            "categorical_sir",
        )
        self.assertEqual(MODULE.resolve_endpoint_profile(None), "categorical_sir")
        self.assertEqual(MODULE.resolve_endpoint_profile({}), "categorical_sir")
        self.assertEqual(
            MODULE.resolve_endpoint_profile({"primary_endpoint": "nonsense"}),
            "categorical_sir",
        )

    def test_categorical_endpoints_read_from_secondary_block(self) -> None:
        config = {
            "secondary_endpoint": {
                "antibiotic_endpoints": {"ciprofloxacin": {"clinical_indication": "non_meningitis"}}
            }
        }
        self.assertEqual(
            MODULE.configured_antibiotic_endpoints(config),
            {"ciprofloxacin": {"clinical_indication": "non_meningitis"}},
        )
        legacy = {"antibiotic_endpoints": {"gentamicin": {"clinical_indication": "urinary_tract"}}}
        self.assertEqual(
            MODULE.configured_antibiotic_endpoints(legacy),
            {"gentamicin": {"clinical_indication": "urinary_tract"}},
        )

    def test_mic_profile_skips_breakpoint_artifact_validation(self) -> None:
        failures = MODULE.validate_breakpoint_artifact(
            mic_cohort(), {"ast": {}}, Path("."), profile="mic_regression"
        )
        self.assertEqual(failures, [])

    def test_categorical_profile_still_requires_frozen_breakpoint_artifact(self) -> None:
        failures = MODULE.validate_breakpoint_artifact(
            valid_cohort(), {"ast": {}}, Path("."), profile="categorical_sir"
        )
        self.assertTrue(failures)


class GenomeProvenanceTests(unittest.TestCase):
    """Raw-read genomes have no NCBI assembly accession; the contract must say so honestly."""

    def _raw(self, cohort):
        cohort = cohort.copy()
        cohort.loc[2, "genome_source"] = "raw_reads"
        cohort.loc[2, "assembly_accession"] = ""
        cohort.loc[2, "genome_source_accession"] = "DRR388395"
        return cohort

    def test_raw_read_genome_with_run_accession_passes(self) -> None:
        report = MODULE.validate_cohort(self._raw(mic_cohort()), min_per_category=1, profile="mic_regression")
        self.assertEqual(report["status"], "PASS", report["failures"])

    def test_raw_read_genome_needs_a_valid_run_accession(self) -> None:
        cohort = self._raw(mic_cohort())
        cohort.loc[2, "genome_source_accession"] = "not-a-run"
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertTrue(any("valid run accession" in f for f in report["failures"]))

    def test_registered_genome_still_needs_its_assembly_accession(self) -> None:
        cohort = mic_cohort()
        cohort.loc[0, "assembly_accession"] = ""
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertEqual(report["status"], "FAIL")

    def test_registered_source_accession_must_match(self) -> None:
        cohort = mic_cohort()
        cohort.loc[0, "genome_source_accession"] = "GCA_999999999.1"
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertTrue(any("source accession differs" in f for f in report["failures"]))

    def test_genome_checksum_must_be_sha256(self) -> None:
        cohort = mic_cohort()
        cohort.loc[1, "genome_sha256"] = "short"
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertTrue(any("genome SHA-256" in f for f in report["failures"]))

    def test_unknown_genome_source_is_rejected(self) -> None:
        cohort = mic_cohort()
        cohort.loc[1, "genome_source"] = "downloaded_somehow"
        report = MODULE.validate_cohort(cohort, min_per_category=1, profile="mic_regression")
        self.assertTrue(any("Invalid genome sources" in f for f in report["failures"]))
