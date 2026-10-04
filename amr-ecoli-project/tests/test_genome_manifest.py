from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "scripts" / "build_genome_manifest.py"
SPEC = importlib.util.spec_from_file_location("build_genome_manifest", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class GenomeRequestTests(unittest.TestCase):
    def test_request_deduplicates_antibiotic_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cohort = Path(tmp) / "cohort.csv"
            write_csv(
                cohort,
                [
                    {
                        "isolate_id": "iso-1",
                        "biosample_accession": "SAMN00000001",
                        "assembly_accession": "GCF_000000001.1",
                        "antibiotic": "ciprofloxacin",
                    },
                    {
                        "isolate_id": "iso-1",
                        "biosample_accession": "SAMN00000001",
                        "assembly_accession": "GCF_000000001.1",
                        "antibiotic": "ceftriaxone",
                    },
                ],
            )
            requests, exclusions = MODULE.build_requests(cohort)
            self.assertEqual(len(requests), 1)
            self.assertEqual(requests[0]["n_ast_records"], 2)
            self.assertEqual(requests[0]["assembly_accession"], "GCF_000000001.1")
            self.assertEqual(exclusions, [])

    def test_ambiguous_isolate_and_assembly_mappings_are_excluded(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cohort = Path(tmp) / "cohort.csv"
            write_csv(
                cohort,
                [
                    {
                        "isolate_id": "iso-1",
                        "assembly_accession": "GCF_000000001.1",
                    },
                    {
                        "isolate_id": "iso-1",
                        "assembly_accession": "GCF_000000002.1",
                    },
                    {
                        "isolate_id": "iso-2",
                        "assembly_accession": "GCF_000000002.1",
                    },
                ],
            )
            requests, exclusions = MODULE.build_requests(cohort)
            self.assertEqual(requests, [])
            reasons = ";".join(row["reason"] for row in exclusions)
            self.assertIn("ISOLATE_MULTIPLE_ASSEMBLIES", reasons)
            self.assertIn("ASSEMBLY_MULTIPLE_ISOLATES", reasons)


class GenomeStagingTests(unittest.TestCase):
    def _request_table(self, root: Path, accessions: list[str]) -> Path:
        table = root / "requests.tsv"
        write_tsv(
            table,
            [
                {
                    "assembly_accession": accession,
                    "isolate_id": f"iso-{index}",
                    "biosample_accession": f"SAMN{index:08d}",
                    "n_ast_records": 1,
                }
                for index, accession in enumerate(accessions, start=1)
            ],
            MODULE.REQUEST_FIELDS,
        )
        return table

    def test_stage_validates_official_layout_and_records_missing_accession(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            present = "GCF_000000001.1"
            missing = "GCF_000000002.1"
            requests = self._request_table(root, [present, missing])
            data = root / "package" / "ncbi_dataset" / "data"
            assembly_dir = data / present
            assembly_dir.mkdir(parents=True)
            report = data / "assembly_data_report.jsonl"
            report.write_text(json.dumps({"accession": present}) + "\n", encoding="utf-8")
            (data / "dataset_catalog.json").write_text("{}\n", encoding="utf-8")
            fasta = assembly_dir / f"{present}_genomic.fna"
            fasta.write_text(">contig-1\nACGTNN\n>contig-2\nACGT\n", encoding="ascii")

            staged, exclusions, metadata = MODULE.stage_package(
                requests, root / "package", root / "staged"
            )
            by_accession = {row["assembly_accession"]: row for row in staged}
            self.assertEqual(by_accession[present]["stage_status"], "STAGED")
            self.assertEqual(by_accession[present]["fasta_contigs"], 2)
            self.assertEqual(by_accession[present]["fasta_total_bases"], 10)
            self.assertEqual(by_accession[present]["fasta_ambiguous_bases"], 2)
            self.assertEqual(by_accession[present]["fasta_sha256"], MODULE.sha256_file(fasta))
            self.assertEqual(exclusions[0]["reason"], "ACCESSION_ABSENT_FROM_ASSEMBLY_REPORT")
            self.assertEqual(metadata["reported_assembly_accessions"], 1)

    def test_stage_rejects_non_dataset_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            requests = self._request_table(root, ["GCF_000000001.1"])
            with self.assertRaisesRegex(ValueError, "NCBI Datasets genome package"):
                MODULE.stage_package(requests, root / "empty", root / "staged")


class GenomeFinalizationTests(unittest.TestCase):
    def _stage_table(self, root: Path, accession: str) -> Path:
        fasta = root / "staged" / accession / f"{accession}.fna"
        fasta.parent.mkdir(parents=True)
        fasta.write_text(">contig\nACGT\n", encoding="ascii")
        table = root / "staging.tsv"
        write_tsv(
            table,
            [
                {
                    "assembly_accession": accession,
                    "isolate_id": "iso-1",
                    "biosample_accession": "SAMN00000001",
                    "n_ast_records": 1,
                    "source_fasta": str(fasta),
                    "staged_fasta": str(fasta),
                    "fasta_sha256": MODULE.sha256_file(fasta),
                    "fasta_size_bytes": fasta.stat().st_size,
                    "fasta_contigs": 1,
                    "fasta_total_bases": 4,
                    "fasta_ambiguous_bases": 0,
                    "stage_status": "STAGED",
                    "exclusion_reason": "",
                }
            ],
            MODULE.STAGE_FIELDS,
        )
        return table

    def _tool_outputs(self, root: Path, accession: str, mlst_status: str = "PERFECT") -> None:
        quast = root / "quast" / accession
        quast.mkdir(parents=True)
        (quast / "exit_code.txt").write_text("0\n", encoding="ascii")
        (quast / "report.tsv").write_text(
            "Assembly\tgenome\n"
            "# contigs\t80\n"
            "Total length\t5000000\n"
            "N50\t75000\n"
            "GC (%)\t50.6\n"
            "# N's per 100 kbp\t2.0\n",
            encoding="utf-8",
        )
        mlst = root / "mlst" / accession
        mlst.mkdir(parents=True)
        (mlst / "exit_code.txt").write_text("0\n", encoding="ascii")
        (mlst / "mlst.tsv").write_text(
            "FILE\tSCHEME\tST\tSTATUS\tSCORE\tALLELES\n"
            f"genome.fna\tecoli_achtman_4\t131\t{mlst_status}\t100\tadk(53);fumC(40)\n",
            encoding="utf-8",
        )
        amr = root / "amr" / accession
        amr.mkdir(parents=True)
        (amr / "exit_code.txt").write_text("0\n", encoding="ascii")
        (amr / "amrfinder.tsv").write_text(
            "Name\tProtein identifier\tGene symbol\n"
            f"{accession}\tprot-1\tblaCTX-M\n",
            encoding="utf-8",
        )

    def test_mlst_nonperfect_is_a_warning_not_a_genomic_exclusion(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            accession = "GCF_000000001.1"
            stage = self._stage_table(root, accession)
            self._tool_outputs(root, accession, mlst_status="NOVEL")
            rows = MODULE.finalize_analyses(
                stage,
                root / "quast",
                root / "mlst",
                root / "amr",
                min_total_length=4_000_000,
                max_total_length=6_500_000,
                max_contigs=500,
                min_n50=20_000,
                max_ns_per_100kb=1000.0,
            )
            self.assertEqual(rows[0]["overall_status"], "WARN")
            self.assertEqual(rows[0]["eligible_for_modeling"], "true")
            self.assertEqual(rows[0]["warning_reasons"], "MLST_NOT_PERFECT")
            self.assertEqual(rows[0]["amrfinder_hit_count"], 1)

    def test_quast_threshold_and_amrfinder_failure_are_explicit_exclusions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            accession = "GCF_000000001.1"
            stage = self._stage_table(root, accession)
            self._tool_outputs(root, accession)
            quast_report = root / "quast" / accession / "report.tsv"
            quast_report.write_text(
                "# contigs\t700\nTotal length\t3000000\nN50\t5000\n# N's per 100 kbp\t2\n",
                encoding="utf-8",
            )
            (root / "amr" / accession / "exit_code.txt").write_text("7\n", encoding="ascii")
            rows = MODULE.finalize_analyses(
                stage,
                root / "quast",
                root / "mlst",
                root / "amr",
                min_total_length=4_000_000,
                max_total_length=6_500_000,
                max_contigs=500,
                min_n50=20_000,
                max_ns_per_100kb=1000.0,
            )
            reasons = rows[0]["exclusion_reasons"]
            self.assertEqual(rows[0]["overall_status"], "EXCLUDE")
            self.assertEqual(rows[0]["eligible_for_modeling"], "false")
            self.assertIn("ASSEMBLY_LENGTH_BELOW_MIN", reasons)
            self.assertIn("ASSEMBLY_CONTIG_COUNT_ABOVE_MAX", reasons)
            self.assertIn("ASSEMBLY_N50_BELOW_MIN", reasons)
            self.assertIn("AMRFINDER_EXECUTION_FAILED", reasons)


if __name__ == "__main__":
    unittest.main()
