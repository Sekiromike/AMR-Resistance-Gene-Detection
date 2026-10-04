from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "scripts" / "apply_breakpoints.py"
SPEC = importlib.util.spec_from_file_location("apply_breakpoints", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


AST_FIELDS = [
    "source_ast_record_id",
    "antibiotic",
    "submitted_ast_category",
    "ast_measurement_type",
    "measurement_sign",
    "ast_value",
    "ast_unit",
    "ast_method",
    "disk_content",
    "ast_testing_date",
    "breakpoint_eligible",
    "breakpoint_ineligibility_reasons",
    "clinical_indication",
]

RULE_FIELDS = [
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
    "disk_content",
    "clinical_indication",
]


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_normalization_manifest(
    path: Path, ast_path: Path, rows: list[dict[str, str]]
) -> None:
    path.write_text(
        json.dumps(
            {
                "scientific_status": "SOURCE_NORMALIZATION_ONLY",
                "mapping_version": "2.0.0",
                "normalizer": {"script_sha256": "1" * 64},
                "acquisition_provenance": {"sha256": "2" * 64},
                "n_breakpoint_eligible_rows": sum(
                    str(row["breakpoint_eligible"]).casefold() == "true" for row in rows
                ),
                "outputs": {
                    "normalized": {
                        "sha256": hashlib.sha256(ast_path.read_bytes()).hexdigest(),
                        "rows": len(rows),
                        "bytes": ast_path.stat().st_size,
                    }
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )


def ast_row(
    record_id: str,
    measurement_type: str,
    comparator: str,
    value: str,
    submitted: str,
    unit: str,
    **overrides: str,
) -> dict[str, str]:
    row = {
        "source_ast_record_id": record_id,
        "antibiotic": "ciprofloxacin",
        "submitted_ast_category": submitted,
        "ast_measurement_type": measurement_type,
        "measurement_sign": comparator,
        "ast_value": value,
        "ast_unit": unit,
        "ast_method": "disk diffusion" if measurement_type == "ZONE" else "broth microdilution",
        "disk_content": "5 ug" if measurement_type == "ZONE" else "",
        "ast_testing_date": "2024-02-04",
        "breakpoint_eligible": "true",
        "breakpoint_ineligibility_reasons": "",
        "clinical_indication": "",
    }
    row.update(overrides)
    return row


def rule_row(
    artifact_sha: str,
    measurement_type: str,
    susceptible: str,
    resistant: str,
    unit: str,
    **overrides: str,
) -> dict[str, str]:
    row = {
        "breakpoint_rule_id": f"EUCAST-16.1-CIP-{measurement_type}",
        "breakpoint_artifact_sha256": artifact_sha,
        "breakpoint_standard": "EUCAST",
        "breakpoint_version": "16.1",
        "antibiotic": "ciprofloxacin",
        "ast_measurement_type": measurement_type,
        "ast_unit": unit,
        "susceptible_breakpoint": susceptible,
        "resistant_breakpoint": resistant,
        "ast_method": "",
        "disk_content": "",
        "clinical_indication": "",
    }
    row.update(overrides)
    return row


class BreakpointInterpretationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.artifact = self.root / "EUCAST_v16.1.xlsx"
        self.artifact.write_bytes(b"frozen breakpoint artifact fixture\n")
        self.artifact_sha = hashlib.sha256(self.artifact.read_bytes()).hexdigest()

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def load_rules(self, rows: list[dict[str, str]]) -> list[object]:
        rules_path = self.root / "rules.csv"
        write_csv(rules_path, RULE_FIELDS, rows)
        return MODULE.load_breakpoint_rules(rules_path, self.artifact)

    def test_mic_interval_semantics_and_inverted_submitted_label(self) -> None:
        rules = self.load_rules([rule_row(self.artifact_sha, "MIC", "1", "2", "mg/L")])
        rows = [
            ast_row("exact-s", "MIC", "=", "1", "R", "mg/L"),
            ast_row("exact-i", "MIC", "=", "1.5", "I", "mg/L"),
            ast_row("r-boundary-is-i", "MIC", "=", "2", "I", "mg/L"),
            ast_row("censored-r", "MIC", ">", "2", "S", "mg/L"),
            ast_row("censored-s", "MIC", "<=", "1", "S", "mg/L"),
            ast_row("cross-s-i", "MIC", "<", "2", "I", "mg/L"),
            ast_row("cross-i-r", "MIC", ">=", "2", "R", "mg/L"),
            ast_row("cross-wide", "MIC", ">", "1", "R", "mg/L"),
        ]

        interpreted, excluded = MODULE.interpret_rows(rows, rules)
        categories = {row["source_ast_record_id"]: row["ast_category"] for row in interpreted}
        self.assertEqual(
            categories,
            {
                "exact-s": "S",
                "exact-i": "I",
                "r-boundary-is-i": "I",
                "censored-r": "R",
                "censored-s": "S",
            },
        )
        exact_s = next(row for row in interpreted if row["source_ast_record_id"] == "exact-s")
        self.assertEqual(exact_s["submitted_ast_category"], "R")
        self.assertEqual(exact_s["submitted_category_concordance"], "DISCORDANT")
        self.assertEqual(exact_s["measurement_sign"], "=")
        self.assertEqual(exact_s["ast_value"], "1")
        self.assertEqual(exact_s["breakpoint_rule_id"], "EUCAST-16.1-CIP-MIC")
        self.assertRegex(exact_s["breakpoint_rule_row_sha256"], r"^[0-9a-f]{64}$")
        self.assertEqual(
            {row["exclusion_reason"] for row in excluded},
            {"censored_measurement_crosses_breakpoints"},
        )
        self.assertEqual(len(excluded), 3)

    def test_zone_interval_semantics_at_censored_boundaries(self) -> None:
        rules = self.load_rules([rule_row(self.artifact_sha, "ZONE", "21", "17", "mm")])
        rows = [
            ast_row("exact-s", "ZONE", "=", "21", "S", "mm"),
            ast_row("censored-s", "ZONE", ">=", "21", "S", "mm"),
            ast_row("censored-r", "ZONE", "<", "17", "R", "mm"),
            ast_row("exact-i", "ZONE", "=", "17", "I", "mm"),
            ast_row("bounded-r", "ZONE", "<=", "16", "R", "mm"),
            ast_row("cross-r-i", "ZONE", "<=", "17", "R", "mm"),
            ast_row("cross-i-s", "ZONE", ">", "17", "S", "mm"),
        ]

        interpreted, excluded = MODULE.interpret_rows(rows, rules)
        categories = {row["source_ast_record_id"]: row["ast_category"] for row in interpreted}
        self.assertEqual(
            categories,
            {
                "exact-s": "S",
                "censored-s": "S",
                "censored-r": "R",
                "exact-i": "I",
                "bounded-r": "R",
            },
        )
        self.assertEqual(
            {row["source_ast_record_id"] for row in excluded},
            {"cross-r-i", "cross-i-s"},
        )

    def test_method_disk_content_indication_and_comparator_are_hard_constraints(self) -> None:
        rules = self.load_rules(
            [
                rule_row(
                    self.artifact_sha,
                    "ZONE",
                    "21",
                    "17",
                    "mm",
                    ast_method="disk diffusion",
                    disk_content="5 ug",
                    clinical_indication="uncomplicated urinary tract",
                )
            ]
        )
        common = {
            "ast_method": "Disk  Diffusion",
            "disk_content": "5 UG",
            "clinical_indication": "Uncomplicated Urinary Tract",
        }
        rows = [
            ast_row("match", "ZONE", "=", "25", "S", "mm", **common),
            ast_row(
                "missing-potency",
                "ZONE",
                "=",
                "25",
                "S",
                "mm",
                ast_method="disk diffusion",
                disk_content="",
                breakpoint_eligible="false",
                breakpoint_ineligibility_reasons="missing_disk_content",
                clinical_indication="uncomplicated urinary tract",
            ),
            ast_row(
                "wrong-indication",
                "ZONE",
                "=",
                "25",
                "S",
                "mm",
                ast_method="disk diffusion",
                disk_content="5 ug",
                clinical_indication="meningitis",
            ),
            ast_row("missing-comparator", "ZONE", "", "25", "S", "mm", **common),
        ]

        interpreted, excluded = MODULE.interpret_rows(rows, rules)
        self.assertEqual([row["source_ast_record_id"] for row in interpreted], ["match"])
        reasons = {row["source_ast_record_id"]: row["exclusion_reason"] for row in excluded}
        self.assertEqual(reasons["missing-potency"], "source_stage_breakpoint_ineligible")
        self.assertEqual(reasons["wrong-indication"], "no_matching_breakpoint_rule")
        self.assertEqual(reasons["missing-comparator"], "missing_measurement_sign")

    def test_every_rule_must_match_frozen_artifact_sha(self) -> None:
        rules_path = self.root / "rules.csv"
        valid = rule_row(self.artifact_sha, "MIC", "1", "2", "mg/L")
        invalid = rule_row("0" * 64, "ZONE", "21", "17", "mm")
        write_csv(rules_path, RULE_FIELDS, [valid, invalid])

        with self.assertRaisesRegex(MODULE.BreakpointError, "artifact SHA-256 mismatch"):
            MODULE.load_breakpoint_rules(rules_path, self.artifact)

    def test_end_to_end_writes_provenance_manifest_and_exclusion_ledger(self) -> None:
        ast_path = self.root / "normalized.csv"
        normalization_manifest_path = self.root / "normalization_manifest.json"
        rules_path = self.root / "rules.csv"
        output_path = self.root / "interpreted.csv"
        exclusions_path = self.root / "excluded.csv"
        manifest_path = self.root / "manifest.json"
        ast_rows = [
            ast_row("resolved", "MIC", "<=", "1", "S", "mg/L"),
            ast_row("unresolved", "MIC", ">=", "2", "R", "mg/L"),
        ]
        write_csv(ast_path, AST_FIELDS, ast_rows)
        write_normalization_manifest(normalization_manifest_path, ast_path, ast_rows)
        write_csv(rules_path, RULE_FIELDS, [rule_row(self.artifact_sha, "MIC", "1", "2", "mg/L")])

        manifest = MODULE.run_interpretation(
            ast_path,
            normalization_manifest_path,
            rules_path,
            self.artifact,
            output_path,
            exclusions_path,
            manifest_path,
        )

        self.assertEqual(manifest["outputs"]["interpreted"]["rows"], 1)
        self.assertEqual(manifest["outputs"]["exclusions"]["rows"], 1)
        self.assertEqual(
            manifest["sources"]["frozen_breakpoint_artifact"]["sha256"], self.artifact_sha
        )
        self.assertEqual(
            manifest["sources"]["normalization_manifest"]["normalized_ast_sha256"],
            hashlib.sha256(ast_path.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            manifest["exclusion_counts"], {"censored_measurement_crosses_breakpoints": 1}
        )
        persisted = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(persisted["category_counts"], {"I": 0, "R": 0, "S": 1})
        with output_path.open("r", encoding="utf-8", newline="") as handle:
            output_rows = list(csv.DictReader(handle))
        self.assertEqual(output_rows[0]["measurement_sign"], "<=")
        self.assertEqual(output_rows[0]["submitted_ast_category"], "S")
        with exclusions_path.open("r", encoding="utf-8", newline="") as handle:
            excluded_rows = list(csv.DictReader(handle))
        self.assertEqual(excluded_rows[0]["measurement_sign"], ">=")

    def test_overlapping_rules_are_rejected_per_row_as_ambiguous(self) -> None:
        first = rule_row(self.artifact_sha, "MIC", "1", "2", "mg/L")
        second = dict(first)
        second["breakpoint_rule_id"] = "EUCAST-16.1-CIP-MIC-DUPLICATE"
        rules = self.load_rules([first, second])

        interpreted, excluded = MODULE.interpret_rows(
            [ast_row("ambiguous", "MIC", "=", "1", "S", "mg/L")], rules
        )
        self.assertEqual(interpreted, [])
        self.assertEqual(excluded[0]["exclusion_reason"], "ambiguous_breakpoint_rules")

    def test_source_stage_ineligible_row_cannot_match_unconstrained_rule(self) -> None:
        rules = self.load_rules([rule_row(self.artifact_sha, "MIC", "1", "2", "mg/L")])
        row = ast_row(
            "ineligible",
            "MIC",
            "=",
            "4",
            "R",
            "mg/L",
            ast_method="",
            ast_testing_date="",
            breakpoint_eligible="false",
            breakpoint_ineligibility_reasons="missing_ast_method;missing_ast_testing_date",
        )

        interpreted, excluded = MODULE.interpret_rows([row], rules)

        self.assertEqual(interpreted, [])
        self.assertEqual(excluded[0]["exclusion_reason"], "source_stage_breakpoint_ineligible")

    def test_forged_eligible_flag_cannot_bypass_assay_metadata_validation(self) -> None:
        rules = self.load_rules([rule_row(self.artifact_sha, "MIC", "1", "2", "mg/L")])
        row = ast_row(
            "forged",
            "MIC",
            "=",
            "4",
            "R",
            "mg/L",
            ast_method="unknown-proxy",
            ast_testing_date="not-a-date",
        )

        interpreted, excluded = MODULE.interpret_rows([row], rules)

        self.assertEqual(interpreted, [])
        self.assertEqual(excluded[0]["exclusion_reason"], "invalid_breakpoint_eligibility_state")

    def test_clinical_indication_context_requires_exact_identity(self) -> None:
        source = ast_row("context", "MIC", "=", "0.25", "S", "mg/L")
        source.update(
            {
                "isolate_id": "PDT1",
                "biosample_accession": "SAMN00000001",
                "assembly_accession": "GCA_000000001.1",
            }
        )
        context = {
            "isolate_id": "PDT1",
            "biosample_accession": "SAMN00000001",
            "assembly_accession": "GCA_000000001.1",
            "clinical_indication": "non_meningitis",
        }
        enriched, excluded = MODULE.enrich_with_context([source], [context])
        self.assertEqual(excluded, [])
        self.assertEqual(enriched[0]["clinical_indication"], "non_meningitis")

        context["assembly_accession"] = "GCA_000000002.1"
        enriched, excluded = MODULE.enrich_with_context([source], [context])
        self.assertEqual(enriched, [])
        self.assertEqual(
            excluded[0]["exclusion_reason"], "breakpoint_context_identity_mismatch"
        )


if __name__ == "__main__":
    unittest.main()
