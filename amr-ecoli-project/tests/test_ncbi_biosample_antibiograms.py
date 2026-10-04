from __future__ import annotations

import importlib.util
import http.client
import json
import sys
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
MODULE_PATH = SCRIPTS / "acquire_ncbi_biosample_antibiograms.py"
SPEC = importlib.util.spec_from_file_location(
    "acquire_ncbi_biosample_antibiograms", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def biosample_xml(*, duplicate: bool = False, omit_second: bool = False) -> bytes:
    second = "" if omit_second else """
  <BioSample submission_date="2020-02-01" last_update="2021-02-01">
    <Ids><Id db="BioSample">SAMN00000002</Id></Ids>
    <Description><Comment><Paragraph>No antibiogram supplied.</Paragraph></Comment></Description>
  </BioSample>"""
    first_accession = "SAMN00000002" if duplicate else "SAMN00000001"
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<BioSampleSet>
  <BioSample submission_date="2020-01-01" last_update="2021-01-01">
    <Ids><Id db="BioSample">{first_accession}</Id><Id db="Sample name">isolate-1</Id></Ids>
    <Description><Comment>
      <Table>
        <Caption>Antibiogram</Caption>
        <Header>
          <Cell>Antibiotic</Cell><Cell>Resistance phenotype</Cell>
          <Cell>Measurement sign</Cell><Cell>Measurement</Cell>
          <Cell>Measurement units</Cell><Cell>Laboratory typing method</Cell>
          <Cell>Laboratory typing platform</Cell><Cell>Vendor</Cell>
          <Cell>Laboratory typing method version or reagent</Cell>
          <Cell>Testing standard</Cell><Cell>Unexpected provenance</Cell>
        </Header>
        <Body>
          <Row><Cell>ciprofloxacin</Cell><Cell>resistant</Cell><Cell>&gt;=</Cell>
            <Cell>4</Cell><Cell>mg/L</Cell><Cell>MIC</Cell><Cell>Vitek</Cell>
            <Cell>Biomérieux</Cell><Cell>AST-N391</Cell><Cell>CLSI</Cell><Cell>raw-x</Cell></Row>
          <Row><Cell>gentamicin</Cell><Cell>susceptible</Cell><Cell>&lt;=</Cell>
            <Cell>1</Cell><Cell>mg/L</Cell><Cell>MIC</Cell><Cell>Phoenix</Cell>
            <Cell>Becton Dickinson</Cell><Cell></Cell><Cell>CLSI</Cell><Cell></Cell></Row>
        </Body>
      </Table>
    </Comment></Description>
  </BioSample>{second}
</BioSampleSet>
""".encode("utf-8")


def overlay_row(**updates: str) -> dict[str, str]:
    row = {field: "" for field in MODULE.OVERLAY_COLUMNS}
    row.update(
        {
            "biosample_accession": "SAMN00000001",
            "antibiogram_row_id": "row-1",
            "antibiotic_raw": "ciprofloxacin",
            "resistance_phenotype_raw": "resistant",
            "measurement_sign_raw": ">=",
            "measurement_raw": "4",
            "measurement_unit_raw": "mg/L",
            "laboratory_typing_method_raw": "MIC",
            "laboratory_typing_platform_raw": "Vitek",
            "vendor_raw": "Biomérieux",
            "laboratory_typing_method_version_or_reagent_raw": "AST-N391",
            "testing_standard_raw": "CLSI",
        }
    )
    row.update(updates)
    return row


def parent_row(**updates: str) -> dict[str, str]:
    row = {
        "id": "ast-1",
        "checksum": "checksum-1",
        "target_acc": "PDT000000001.1",
        "biosample_acc": "SAMN00000001",
        "antibiotic": "ciprofloxacin",
        "phenotype": "resistant",
        "measurement_sign": ">=",
        "mic": "4.0",
        "disk_diffusion": "",
        "standard": "CLSI",
        "platform": "Vitek",
        "vendor": "Biomérieux",
        "reagent": "AST-N391",
    }
    row.update(updates)
    return row


class NcbiBioSampleAntibiogramTests(unittest.TestCase):
    def test_request_plan_is_sorted_batched_and_hash_frozen(self) -> None:
        plan = MODULE.build_request_plan(
            ["SAMN00000003", "SAMN00000001", "SAMN00000002"], 2
        )
        self.assertEqual(plan["batch_count"], 2)
        self.assertEqual(
            plan["batches"][0]["accessions"], ["SAMN00000001", "SAMN00000002"]
        )
        self.assertEqual(plan["batches"][1]["accessions"], ["SAMN00000003"])
        self.assertEqual(
            plan["accession_list_sha256"],
            MODULE.accession_list_sha256(
                ["SAMN00000001", "SAMN00000002", "SAMN00000003"]
            ),
        )
        self.assertTrue(
            plan["parameters_without_private_values"]["email_supplied_at_runtime"]
        )
        self.assertNotIn("researcher@example.org", json.dumps(plan))

    def test_request_plan_rejects_duplicates_and_oversized_batches(self) -> None:
        with self.assertRaises(MODULE.BioSampleOverlayError):
            MODULE.build_request_plan(["SAMN00000001", "SAMN00000001"], 2)
        with self.assertRaises(MODULE.BioSampleOverlayError):
            MODULE.build_request_plan(["SAMN00000001"], MODULE.MAX_BATCH_SIZE + 1)

    def test_incomplete_chunked_response_is_retryable_acquisition_error(self) -> None:
        response = mock.MagicMock()
        response.__enter__.return_value.read.side_effect = http.client.IncompleteRead(
            b"partial"
        )
        with mock.patch.object(MODULE.urllib.request, "urlopen", return_value=response):
            with self.assertRaisesRegex(
                MODULE.BioSampleOverlayError, "NCBI EFetch request failed"
            ):
                MODULE.fetch_batch(
                    ["SAMN00000001"],
                    email="contact@example.org",
                    api_key="",
                    timeout=1.0,
                )

    def test_xml_parser_preserves_raw_rows_and_missing_nonproxy_fields(self) -> None:
        rows, summary = MODULE.parse_biosample_xml(
            biosample_xml(),
            expected_accessions=["SAMN00000001", "SAMN00000002"],
            raw_batch_file="raw_xml/batch_00001.xml",
        )
        self.assertEqual(summary["biosample_records"], 2)
        self.assertEqual(summary["biosamples_with_antibiograms"], 1)
        self.assertEqual(summary["antibiogram_rows"], 2)
        self.assertEqual(summary["unknown_headers"], {"Unexpected provenance": 1})
        self.assertEqual(rows[0]["measurement_raw"], "4")
        self.assertEqual(rows[0]["laboratory_typing_method_raw"], "MIC")
        self.assertEqual(
            rows[0]["laboratory_typing_method_version_or_reagent_raw"], "AST-N391"
        )
        self.assertEqual(rows[0]["testing_standard_version_raw"], "")
        self.assertEqual(rows[0]["ast_testing_date_raw"], "")
        self.assertEqual(
            json.loads(rows[0]["unmapped_cells_json"]),
            {"Unexpected provenance": "raw-x"},
        )

    def test_xml_parser_fails_closed_on_accession_mismatch(self) -> None:
        with self.assertRaisesRegex(
            MODULE.BioSampleOverlayError, "batch identity mismatch"
        ):
            MODULE.parse_biosample_xml(
                biosample_xml(omit_second=True),
                expected_accessions=["SAMN00000001", "SAMN00000002"],
                raw_batch_file="raw_xml/batch_00001.xml",
            )
        with self.assertRaises(MODULE.BioSampleOverlayError):
            MODULE.parse_biosample_xml(
                biosample_xml(duplicate=True),
                expected_accessions=["SAMN00000001", "SAMN00000002"],
                raw_batch_file="raw_xml/batch_00001.xml",
            )

    def test_reconciliation_matches_numeric_values_without_selecting_ambiguity(self) -> None:
        rows, summary = MODULE.reconcile_ast_rows(
            [parent_row()], [overlay_row()]
        )
        self.assertEqual(rows[0]["reconciliation_status"], "unique_core_match")
        self.assertEqual(rows[0]["unique_core_match_row_id"], "row-1")
        self.assertEqual(rows[0]["standard_exact_when_present"], "True")
        self.assertEqual(summary["status_counts"], {"unique_core_match": 1})

        ambiguous, summary = MODULE.reconcile_ast_rows(
            [parent_row()],
            [overlay_row(), overlay_row(antibiogram_row_id="row-2")],
        )
        self.assertEqual(ambiguous[0]["reconciliation_status"], "ambiguous_core_match")
        self.assertEqual(ambiguous[0]["unique_core_match_row_id"], "")
        self.assertEqual(
            json.loads(ambiguous[0]["core_match_row_ids_json"]), ["row-1", "row-2"]
        )
        self.assertEqual(summary["status_counts"], {"ambiguous_core_match": 1})

    def test_reconciliation_does_not_use_collection_or_record_dates_as_ast_date(self) -> None:
        source = overlay_row(
            biosample_record_submission_date_raw="2020-01-01",
            biosample_record_last_update_raw="2021-01-01",
        )
        coverage = MODULE.field_coverage([source])
        self.assertEqual(coverage["ast_testing_date_raw"]["present"], 0)
        self.assertEqual(coverage["testing_standard_version_raw"]["present"], 0)


if __name__ == "__main__":
    unittest.main()
