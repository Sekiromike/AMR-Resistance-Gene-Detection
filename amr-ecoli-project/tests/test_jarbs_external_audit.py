from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
MODULE_PATH = SCRIPTS / "audit_jarbs_external.py"
SPEC = importlib.util.spec_from_file_location("audit_jarbs_external", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


HEADER = [
    "isolate ID",
    "genome_species",
    "CTRX",
    "CPFX",
    "GM",
    "Mean_total_coverage",
    "Num_contigs",
    "Total_bases",
    "MLST",
]


def supplement_row(
    isolate_id: str,
    species: str = "Escherichia coli",
    ctrx: str = "<=1",
    cpfx: str = ">4",
    gm: str = "2",
) -> list[str]:
    return [isolate_id, species, ctrx, cpfx, gm, "50", "100", "5000000", "131"]


def run_row(
    isolate_id: str,
    run: str,
    platform: str = "ILLUMINA",
    biosample: str = "SAMD00000001",
) -> dict[str, str]:
    return {
        "isolate_id": isolate_id,
        "scientific_name": "Escherichia coli",
        "taxid": "562",
        "biosample_accession": biosample,
        "sra_sample_accession": "DRS000001",
        "experiment_accession": "DRX000001",
        "run_accession": run,
        "bioproject_accession": "PRJDB10842",
        "submission_accession": "DRA000001",
        "platform": platform,
        "instrument_model": "instrument",
        "library_layout": "PAIRED" if platform == "ILLUMINA" else "SINGLE",
        "collection_date_raw": "2019",
        "geo_loc_name_raw": "Japan:Tokyo",
        "isolation_source_raw": "urine",
        "run_published_raw": "2023-01-01",
        "run_total_bases": "1000",
        "run_file_md5": "a" * 32,
    }


class JarbsExternalAuditTests(unittest.TestCase):
    def test_supplement_counts_target_mics_without_filling_missing_values(self) -> None:
        sheets = {
            "Dataset 6": [
                HEADER,
                supplement_row("A"),
                supplement_row("B", ctrx="", cpfx="", gm=""),
                supplement_row("C", species="Klebsiella pneumoniae"),
                ["", "", "", "", "", "", "", "", "n.d.: not determined"],
            ]
        }
        rows, summary = MODULE.audit_supplement(sheets)
        self.assertEqual(len(rows), 3)
        self.assertEqual(summary["ecoli_isolates"], 2)
        self.assertEqual(summary["ecoli_all_three_target_mics_present"], 1)
        self.assertEqual(summary["ecoli_target_drugs"]["ceftriaxone"]["missing"], 1)
        self.assertEqual(summary["ecoli_target_drugs"]["ceftriaxone"]["censored"], 1)
        self.assertFalse(summary["explicit_field_presence"]["ast_testing_date"])

    def test_supplement_rejects_unexpected_target_mic_text(self) -> None:
        sheets = {"Dataset 6": [HEADER, supplement_row("A", ctrx="resistant")]}
        with self.assertRaisesRegex(MODULE.JarbsAuditError, "Unexpected target-drug MIC"):
            MODULE.audit_supplement(sheets)

    def test_linkage_separates_short_and_long_reads(self) -> None:
        supplement = [
            dict(zip(HEADER, supplement_row("A"))),
            dict(zip(HEADER, supplement_row("B", species="Klebsiella pneumoniae"))),
        ]
        runs = [
            run_row("A", "DRR000001"),
            run_row("B", "DRR000002", biosample="SAMD00000002"),
            run_row("B", "DRR000003", platform="OXFORD_NANOPORE", biosample="SAMD00000002"),
        ]
        manifest, summary = MODULE.build_isolate_sequence_manifest(supplement, runs)
        self.assertEqual(manifest[0]["primary_short_read_run"], "DRR000001")
        self.assertEqual(json.loads(manifest[1]["long_read_runs_json"]), ["DRR000003"])
        self.assertEqual(summary["ecoli_link_status_counts"], {"linked_one_short_read_run": 1})

    def test_exact_overlap_audit_covers_biosample_run_assembly_and_isolate(self) -> None:
        sequence = [
            {
                "isolate_id": "JARBS-1",
                "biosample_accession": "SAMD00000001",
                "all_runs_json": '["DRR000001"]',
            }
        ]
        assemblies = [{"assembly_accession": "GCA_000000001.1", "strain_raw": "JARBS-1"}]
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            isolates_path = base / "isolates.csv"
            enrichment_path = base / "enrichment.csv"
            with isolates_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["target_acc", "biosample_acc", "asm_acc", "strain"],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "target_acc": "PDT000000001.1",
                        "biosample_acc": "SAMD00000001",
                        "asm_acc": "GCA_000000001.1",
                        "strain": "JARBS-1",
                    }
                )
            with enrichment_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "target_acc",
                        "sra_run_accessions_raw",
                        "alternative_isolate_identifiers_json",
                    ],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "target_acc": "PDT000000001.1",
                        "sra_run_accessions_raw": "DRR000001",
                        "alternative_isolate_identifiers_json": '["JARBS-1"]',
                    }
                )
            overlaps, summary = MODULE.audit_development_overlap(
                sequence, assemblies, isolates_path, enrichment_path
            )
        self.assertEqual(
            {row["identifier_type"] for row in overlaps},
            {
                "biosample_accession",
                "sra_run_accession",
                "assembly_accession",
                "isolate_identifier",
            },
        )
        self.assertEqual(summary["external_lock_status"], "blocked_exact_identifier_overlap")

    def test_registered_assembly_report_preserves_coverage_boundary(self) -> None:
        payload = {
            "reports": [
                {
                    "accession": "GCA_000000001.1",
                    "paired_accession": "GCF_000000001.1",
                    "source_database": "SOURCE_DATABASE_GENBANK",
                    "organism": {
                        "tax_id": 562,
                        "organism_name": "Escherichia coli",
                        "infraspecific_names": {"strain": "JARBS-1"},
                    },
                    "assembly_info": {
                        "assembly_level": "Complete Genome",
                        "assembly_status": "current",
                        "release_date": "2023-01-01",
                        "bioproject_accession": "PRJDB10842",
                    },
                }
            ]
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "report.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            rows, summary = MODULE.parse_assembly_report(path)
        self.assertEqual(rows[0]["assembly_accession"], "GCA_000000001.1")
        self.assertEqual(summary["ecoli_assembly_accessions"], 1)
        self.assertIn("not_all_supplement", summary["coverage_boundary"])


if __name__ == "__main__":
    unittest.main()
