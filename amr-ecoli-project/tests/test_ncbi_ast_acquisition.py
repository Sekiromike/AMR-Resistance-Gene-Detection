from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).parents[1]
MODULE_PATH = ROOT / "scripts" / "acquire_ncbi_ast.py"
SPEC = importlib.util.spec_from_file_location("acquire_ncbi_ast", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class NcbiAstAcquisitionTests(unittest.TestCase):
    def test_queries_use_documented_tables_and_target_join(self) -> None:
        queries = MODULE.build_queries()
        self.assertIn("`ncbi-pathogen-detect.pdbrowser.ast`", queries["ast"])
        self.assertIn("`ncbi-pathogen-detect.pdbrowser.isolates`", queries["isolates"])
        self.assertIn("USING (target_acc)", queries["isolates"])
        self.assertIn("i.creation_date AS target_creation_date", queries["isolates"])
        self.assertIn("@species", queries["ast"])
        self.assertIn("@antibiotics_csv", queries["ast"])
        self.assertNotIn("breakpoint_version", queries["ast"])

    def test_offline_dry_run_needs_no_project_bq_or_credentials(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary) / "plan"
            environment = os.environ.copy()
            environment.pop("GOOGLE_CLOUD_PROJECT", None)
            environment["PATH"] = ""
            result = subprocess.run(
                [
                    sys.executable,
                    str(MODULE_PATH),
                    "--dry-run",
                    "--output-dir",
                    str(output_dir),
                    "--antibiotic",
                    "Ciprofloxacin",
                ],
                cwd=ROOT,
                env=environment,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = json.loads((output_dir / MODULE.MANIFEST_NAME).read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "planned")
            self.assertIsNone(manifest["retrieval"])
            self.assertEqual(
                manifest["selection"]["antibiotics_exact_case_insensitive"], ["ciprofloxacin"]
            )
            query = (output_dir / MODULE.AST_QUERY_NAME).read_text(encoding="utf-8")
            self.assertEqual(manifest["queries"]["ast"]["sha256"], MODULE.sha256_text(query))
            self.assertFalse((output_dir / MODULE.AST_OUTPUT_NAME).exists())

    def test_execution_without_project_fails_before_looking_for_bq(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, mock.patch.dict(
            os.environ, {"GOOGLE_CLOUD_PROJECT": ""}, clear=False
        ), mock.patch.object(MODULE.shutil, "which") as which:
            code = MODULE.main(["--output-dir", temporary])
            self.assertEqual(code, 2)
            which.assert_not_called()

    def test_missing_bq_error_is_actionable(self) -> None:
        queries = MODULE.build_queries()
        manifest = MODULE.base_manifest(
            species="Escherichia coli",
            antibiotics=[],
            queries=queries,
            generated_at="2026-08-15T00:00:00Z",
        )
        with tempfile.TemporaryDirectory() as temporary, mock.patch.object(
            MODULE.shutil, "which", return_value=None
        ):
            with self.assertRaisesRegex(MODULE.AcquisitionError, "bq.*not installed"):
                MODULE.acquire(
                    output_dir=Path(temporary),
                    project="example-project",
                    species="Escherichia coli",
                    antibiotics=[],
                    queries=queries,
                    manifest=manifest,
                )

    def test_cloud_cli_environment_forces_utf8_output(self) -> None:
        environment = MODULE.command_environment()
        self.assertEqual(environment["CLOUDSDK_ENCODING"], "utf-8")
        self.assertEqual(environment["PYTHONIOENCODING"], "utf-8")

    def test_completed_acquisition_records_source_and_output_hashes(self) -> None:
        queries = MODULE.build_queries()
        manifest = MODULE.base_manifest(
            species="Escherichia coli",
            antibiotics=["ciprofloxacin"],
            queries=queries,
            generated_at="2026-08-15T00:00:00Z",
        )
        table_metadata = json.dumps({"id": "source-table", "lastModifiedTime": "123456789"})
        captures = [
            subprocess.CompletedProcess([], 0, stdout=table_metadata, stderr=""),
            subprocess.CompletedProcess([], 0, stdout=table_metadata, stderr=""),
            subprocess.CompletedProcess([], 0, stdout="bq 2.1.0\n", stderr=""),
        ]

        def fake_export(command: list[str], destination: Path, sql: str) -> None:
            self.assertIn("--parameter=species:STRING:Escherichia coli", command)
            self.assertNotIn(sql, command)
            if destination.name == MODULE.AST_OUTPUT_NAME:
                self.assertEqual(sql, queries["ast"])
                destination.write_text(
                    "target_acc,scientific_name,antibiotic,mic\n"
                    "PDT1,Escherichia coli,ciprofloxacin,1\n",
                    encoding="utf-8",
                )
            else:
                self.assertEqual(sql, queries["isolates"])
                destination.write_text(
                    "target_acc,biosample_acc,asm_acc\nPDT1,SAMN1,GCA_000000001.1\n",
                    encoding="utf-8",
                )

        private_bq_path = r"C:\Users\private-user\GoogleCloudSDK\bin\bq.CMD"
        with tempfile.TemporaryDirectory() as temporary, mock.patch.object(
            MODULE.shutil, "which", return_value=private_bq_path
        ), mock.patch.object(MODULE, "run_capture", side_effect=captures), mock.patch.object(
            MODULE, "run_query_to_file", side_effect=fake_export
        ):
            output_dir = Path(temporary) / "raw"
            completed = MODULE.acquire(
                output_dir=output_dir,
                project="example-project",
                species="Escherichia coli",
                antibiotics=["ciprofloxacin"],
                queries=queries,
                manifest=manifest,
            )
            self.assertEqual(completed["status"], "complete")
            self.assertEqual(completed["retrieval"]["outputs"]["ast"]["rows"], 1)
            ast_path = output_dir / MODULE.AST_OUTPUT_NAME
            self.assertEqual(
                completed["retrieval"]["outputs"]["ast"]["sha256"],
                MODULE.sha256_file(ast_path),
            )
            source_path = output_dir / MODULE.SOURCE_METADATA_NAME
            self.assertEqual(
                completed["retrieval"]["source_table_metadata"]["sha256"],
                MODULE.sha256_file(source_path),
            )
            on_disk = json.loads((output_dir / MODULE.MANIFEST_NAME).read_text(encoding="utf-8"))
            self.assertEqual(on_disk["retrieval"]["query_project"], "example-project")
            self.assertEqual(on_disk["retrieval"]["client"]["executable"], "bq.CMD")
            self.assertNotIn("private-user", json.dumps(on_disk))

    def test_csv_row_counter_handles_embedded_newlines_and_schema_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "ast.csv"
            path.write_text(
                'target_acc,scientific_name,antibiotic,note\n'
                'PDT1,Escherichia coli,ciprofloxacin,"line one\nline two"\n',
                encoding="utf-8",
            )
            self.assertEqual(
                MODULE.count_csv_rows(path, {"target_acc", "scientific_name", "antibiotic"}), 1
            )
            with self.assertRaisesRegex(MODULE.AcquisitionError, "schema may have changed"):
                MODULE.count_csv_rows(path, {"mic"})

    def test_antibiotic_filter_rejects_delimiter_injection(self) -> None:
        with self.assertRaises(ValueError):
            MODULE.normalize_antibiotics(["ciprofloxacin|meropenem"])
        with self.assertRaises(ValueError):
            MODULE.normalize_antibiotics(["ciprofloxacin,meropenem"])

    def test_manifest_antibiotics_parameter_matches_bq_command(self) -> None:
        antibiotics = MODULE.normalize_antibiotics(
            ["Gentamicin", "ciprofloxacin", "Ceftriaxone"]
        )
        queries = MODULE.build_queries()
        manifest = MODULE.base_manifest(
            species="Escherichia coli",
            antibiotics=antibiotics,
            queries=queries,
            generated_at="2026-08-15T00:00:00Z",
        )
        command = MODULE.bq_query_command(
            bq="bq.cmd",
            project="example-project",
            species="Escherichia coli",
            antibiotics=antibiotics,
        )
        command_parameter = next(
            argument.removeprefix("--parameter=antibiotics_csv:STRING:")
            for argument in command
            if argument.startswith("--parameter=antibiotics_csv:STRING:")
        )
        manifest_parameter = manifest["queries"]["parameters"]["antibiotics_csv"]
        self.assertEqual(manifest_parameter, "ceftriaxone,ciprofloxacin,gentamicin")
        self.assertEqual(command_parameter, manifest_parameter)

    def test_executable_identity_removes_windows_and_posix_directories(self) -> None:
        self.assertEqual(
            MODULE.executable_identity(r"C:\Users\private-user\bin\bq.CMD"),
            "bq.CMD",
        )
        self.assertEqual(MODULE.executable_identity("/opt/google-cloud-sdk/bin/bq"), "bq")

    def test_bq_command_keeps_sql_out_of_windows_batch_arguments(self) -> None:
        command = MODULE.bq_query_command(
            bq="bq.cmd",
            project="example-project",
            species="Escherichia coli",
            antibiotics=["ciprofloxacin", "ceftriaxone", "gentamicin"],
        )
        self.assertNotIn("|", " ".join(command))
        self.assertIn(
            "--parameter=antibiotics_csv:STRING:ciprofloxacin,ceftriaxone,gentamicin",
            command,
        )
        self.assertFalse(any("SELECT" in argument for argument in command))

    def test_refuses_to_overwrite_acquisition_without_force(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary)
            (output_dir / MODULE.MANIFEST_NAME).write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(MODULE.AcquisitionError, "Refusing to overwrite"):
                MODULE.check_destinations(output_dir, force=False)
            MODULE.check_destinations(output_dir, force=True)


if __name__ == "__main__":
    unittest.main()
