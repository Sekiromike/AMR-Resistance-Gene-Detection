from __future__ import annotations

import csv
import importlib.util
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).parents[1]
MODULE_PATH = ROOT / "scripts" / "acquire_ncbi_isolate_enrichment.py"
SPEC = importlib.util.spec_from_file_location("acquire_ncbi_isolate_enrichment", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class NcbiIsolateEnrichmentTests(unittest.TestCase):
    def make_parent(
        self,
        root: Path,
        *,
        duplicate_target: bool = False,
        invalid_target: bool = False,
    ) -> tuple[Path, Path, Path]:
        parent_dir = root / "parent"
        parent_dir.mkdir(parents=True)
        ast_path = parent_dir / "ncbi_ast.csv"
        isolates_path = parent_dir / "ncbi_isolates.csv"
        manifest_path = parent_dir / "ncbi_acquisition_manifest.json"
        ast_query_path = parent_dir / "ncbi_ast_query.sql"
        isolates_query_path = parent_dir / "ncbi_isolates_query.sql"
        source_metadata_path = parent_dir / "ncbi_source_tables.json"
        ast_path.write_text(
            "target_acc,scientific_name,antibiotic,mic\n"
            "PDT000000002.1,Escherichia coli,ciprofloxacin,1\n"
            "PDT000000001.1,Escherichia coli,gentamicin,2\n",
            encoding="utf-8",
        )
        target_two = "invalid-target" if invalid_target else "PDT000000002.1"
        rows = [
            {
                "target_acc": target_two,
                "biosample_acc": "SAMN00000002",
                "asm_acc": "",
                "scientific_name": "Escherichia coli",
                "species_taxid": "562",
                "bioproject_acc": "PRJNA2",
                "collection_date": "2022",
                "geo_loc_name": "Canada",
            },
            {
                "target_acc": target_two if duplicate_target else "PDT000000001.1",
                "biosample_acc": "SAMN00000001",
                "asm_acc": "GCA_000000001.1",
                "scientific_name": "Escherichia coli",
                "species_taxid": "562",
                "bioproject_acc": "PRJNA1",
                "collection_date": "2021-01-02",
                "geo_loc_name": "USA: MA",
            },
        ]
        with isolates_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        ast_query_path.write_text("SELECT * FROM `example.ast`\n", encoding="utf-8")
        isolates_query_path.write_text(
            "SELECT * FROM `example.isolates`\n", encoding="utf-8"
        )
        source_metadata_path.write_text(
            json.dumps({"tables": ["example.ast", "example.isolates"]}),
            encoding="utf-8",
        )
        manifest = {
            "schema_version": "1.0.0",
            "status": "complete",
            "queries": {
                "ast": {
                    "file": ast_query_path.name,
                    "sha256": MODULE.sha256_file(ast_query_path),
                },
                "isolates": {
                    "file": isolates_query_path.name,
                    "sha256": MODULE.sha256_file(isolates_query_path),
                },
                "parameters": {},
            },
            "retrieval": {
                "client": {"executable": "bq"},
                "source_table_metadata": {
                    "file": source_metadata_path.name,
                    "sha256": MODULE.sha256_file(source_metadata_path),
                },
                "outputs": {
                    "ast": {
                        "file": ast_path.name,
                        "rows": 2,
                        "bytes": ast_path.stat().st_size,
                        "sha256": MODULE.sha256_file(ast_path),
                    },
                    "isolates": {
                        "file": isolates_path.name,
                        "rows": 2,
                        "bytes": isolates_path.stat().st_size,
                        "sha256": MODULE.sha256_file(isolates_path),
                    },
                }
            },
        }
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        return isolates_path, ast_path, manifest_path

    def load_context(self, root: Path, **kwargs: object) -> dict[str, object]:
        isolates, ast, manifest = self.make_parent(root, **kwargs)
        return MODULE.load_parent_context(
            isolates_path=isolates,
            ast_path=ast,
            acquisition_manifest_path=manifest,
        )

    def enrichment_row(
        self,
        target: str,
        biosample: str,
        assembly: str,
        *,
        match_count: str = "1",
    ) -> dict[str, str]:
        row = {field: "" for field in MODULE.OUTPUT_COLUMNS}
        row.update(
            {
                "target_acc": target,
                "source_target_acc": target,
                "source_match_count": match_count,
                "source_biosample_accession": biosample,
                "source_assembly_accession": assembly,
                "source_scientific_name": "Escherichia coli",
                "source_species_taxonomy_id": "562",
                "sra_run_accessions_raw": "SRR123456",
                "submitted_assembly_length_bp": "5000000",
                "submitted_assembly_contig_count": "100",
                "submitted_assembly_contig_n50": "100000",
                "alternative_isolate_identifiers_json": '["public-id"]',
                "source_isolate_checksum": "source-checksum",
            }
        )
        if target == "PDT000000001.1":
            row.update(
                {
                    "source_bioproject_accession": "PRJNA1",
                    "source_collection_date_raw": "2021-01-02",
                    "source_geographic_location_raw": "USA: MA",
                }
            )
        else:
            row.update(
                {
                    "source_bioproject_accession": "PRJNA2",
                    "source_collection_date_raw": "2022",
                    "source_geographic_location_raw": "Canada",
                }
            )
        return row

    def write_enrichment(self, path: Path, rows: list[dict[str, str]]) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MODULE.OUTPUT_COLUMNS)
            writer.writeheader()
            writer.writerows(rows)

    def source_metadata(self) -> dict[str, object]:
        return {
            "type": "TABLE",
            "id": "ncbi-pathogen-detect:pdbrowser.isolates",
            "location": "US",
            "schema": {
                "fields": [
                    {"name": name, "type": field_type, "mode": mode}
                    for name, (field_type, mode) in MODULE.SOURCE_SCHEMA.items()
                ]
            },
        }

    def job_metadata(
        self,
        *,
        query: str,
        job_id: str,
        project: str,
        source_as_of_utc: str,
        maximum_bytes_billed: int,
    ) -> dict[str, object]:
        return {
            "user_email": "private-account@example.org",
            "selfLink": "https://example.invalid/private-job-link",
            "jobReference": {
                "projectId": project,
                "jobId": job_id,
                "location": "US",
            },
            "status": {"state": "DONE"},
            "configuration": {
                "query": {
                    "query": query,
                    "useLegacySql": False,
                    "maximumBytesBilled": str(maximum_bytes_billed),
                    "queryParameters": [
                        {
                            "name": "source_as_of_utc",
                            "parameterType": {"type": "TIMESTAMP"},
                            "parameterValue": {"value": source_as_of_utc},
                        }
                    ],
                }
            },
            "statistics": {
                "creationTime": "1787281200000",
                "startTime": "1787281200100",
                "endTime": "1787281200200",
                "query": {
                    "cacheHit": False,
                    "totalBytesProcessed": "123456",
                    "totalBytesBilled": "10485760",
                },
            },
        }

    def test_parent_targets_are_unique_validated_and_sorted(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            context = self.load_context(Path(temporary))
            self.assertEqual(
                context["targets"], ["PDT000000001.1", "PDT000000002.1"]
            )
            self.assertEqual(context["rows"], {"isolates": 2, "ast": 2})
            self.assertEqual(
                context["target_list_sha256"],
                MODULE.sha256_text("PDT000000001.1\nPDT000000002.1\n"),
            )

    def test_duplicate_and_invalid_frozen_targets_fail_closed(self) -> None:
        for option, message in (
            ({"duplicate_target": True}, "not unique"),
            ({"invalid_target": True}, "Invalid or missing frozen target_acc"),
        ):
            with self.subTest(option=option), tempfile.TemporaryDirectory() as temporary:
                with self.assertRaisesRegex(MODULE.EnrichmentError, message):
                    self.load_context(Path(temporary), **option)

    def test_frozen_parent_requires_exact_ecoli_name_and_taxonomy(self) -> None:
        for old, new, message in (
            ("Escherichia coli", "Escherichia fergusonii", "exact scientific_name"),
            (",562,", ",999,", "exact species_taxid"),
        ):
            with self.subTest(new=new), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                isolates, _, _ = self.make_parent(root)
                isolates.write_text(
                    isolates.read_text(encoding="utf-8").replace(old, new, 1),
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(MODULE.EnrichmentError, message):
                    MODULE.read_frozen_isolates(isolates)

    def test_source_timestamp_is_required_and_offsets_normalize_to_utc(self) -> None:
        with self.assertRaises(SystemExit) as caught:
            MODULE.parse_args(["--output-dir", "unused"])
        self.assertEqual(caught.exception.code, 2)
        self.assertEqual(
            MODULE.normalize_utc_timestamp("2026-08-20T23:00:00-04:00"),
            "2026-08-21T03:00:00Z",
        )

    def test_parent_hash_change_fails_before_query_planning(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            isolates, ast, manifest = self.make_parent(root)
            ast.write_text(ast.read_text(encoding="utf-8") + "\n", encoding="utf-8")
            with self.assertRaisesRegex(MODULE.EnrichmentError, "AST output hash mismatch"):
                MODULE.load_parent_context(
                    isolates_path=isolates,
                    ast_path=ast,
                    acquisition_manifest_path=manifest,
                )

    def test_parent_supporting_artifacts_and_append_only_amendment_are_authenticated(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            isolates, ast, manifest_path = self.make_parent(root)
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["selection"] = {
                "antibiotics_exact_case_insensitive": [
                    "ceftriaxone",
                    "ciprofloxacin",
                    "gentamicin",
                ]
            }
            manifest["queries"]["parameters"] = {
                "antibiotics_csv": "ceftriaxone|ciprofloxacin|gentamicin"
            }
            manifest["retrieval"]["client"]["executable"] = (
                r"C:\Users\private-user\GoogleCloudSDK\bin\bq.CMD"
            )
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            source_metadata = manifest_path.parent / manifest["retrieval"][
                "source_table_metadata"
            ]["file"]
            amendment_path = manifest_path.parent / MODULE.PARENT_AMENDMENT_NAME
            amendment = {
                "amendment_status": "complete",
                "original_manifest": {
                    "file": manifest_path.name,
                    "sha256": MODULE.sha256_file(manifest_path),
                },
                "corrections": [
                    {
                        "json_pointer": "/queries/parameters/antibiotics_csv",
                        "corrected_value": "ceftriaxone,ciprofloxacin,gentamicin",
                    },
                    {
                        "json_pointer": "/retrieval/client/executable",
                        "corrected_value": "bq.CMD",
                    },
                ],
                "unchanged_artifacts": {
                    "ast_sha256": MODULE.sha256_file(ast),
                    "isolates_sha256": MODULE.sha256_file(isolates),
                    "source_table_metadata_sha256": MODULE.sha256_file(source_metadata),
                },
            }
            amendment_path.write_text(json.dumps(amendment), encoding="utf-8")
            with self.assertRaisesRegex(MODULE.EnrichmentError, "no verified append-only"):
                MODULE.load_parent_context(
                    isolates_path=isolates,
                    ast_path=ast,
                    acquisition_manifest_path=manifest_path,
                )
            context = MODULE.load_parent_context(
                isolates_path=isolates,
                ast_path=ast,
                acquisition_manifest_path=manifest_path,
                amendment_path=amendment_path,
            )
            self.assertEqual(
                context["manifest_record"]["acquisition_amendment"]["sha256"],
                MODULE.sha256_file(amendment_path),
            )
            tampered_cases = []
            wrong_original = json.loads(json.dumps(amendment))
            wrong_original["original_manifest"]["sha256"] = "0" * 64
            tampered_cases.append((wrong_original, "does not authenticate"))
            missing_unchanged = json.loads(json.dumps(amendment))
            del missing_unchanged["unchanged_artifacts"]["ast_sha256"]
            tampered_cases.append((missing_unchanged, "disagrees with verified ast_sha256"))
            wrong_correction = json.loads(json.dumps(amendment))
            wrong_correction["corrections"][0]["corrected_value"] = "wrong"
            tampered_cases.append((wrong_correction, "lacks the required correction"))
            for tampered, message in tampered_cases:
                with self.subTest(message=message):
                    amendment_path.write_text(json.dumps(tampered), encoding="utf-8")
                    with self.assertRaisesRegex(MODULE.EnrichmentError, message):
                        MODULE.load_parent_context(
                            isolates_path=isolates,
                            ast_path=ast,
                            acquisition_manifest_path=manifest_path,
                            amendment_path=amendment_path,
                        )
            amendment_path.write_text(json.dumps(amendment), encoding="utf-8")
            query_path = manifest_path.parent / manifest["queries"]["isolates"]["file"]
            query_path.write_text("SELECT 'tampered'\n", encoding="utf-8")
            with self.assertRaisesRegex(MODULE.EnrichmentError, "isolates query hash mismatch"):
                MODULE.load_parent_context(
                    isolates_path=isolates,
                    ast_path=ast,
                    acquisition_manifest_path=manifest_path,
                    amendment_path=amendment_path,
                )

    def test_sql_is_deterministic_literal_only_and_uses_documented_aliases(self) -> None:
        targets = ["PDT000000002.1", "PDT000000001.1"]
        parent_hash = "a" * 64
        first = MODULE.build_query(targets, parent_hash)
        second = MODULE.build_query(list(reversed(targets)), parent_hash)
        self.assertEqual(first, second)
        self.assertLess(first.index("PDT000000001.1"), first.index("PDT000000002.1"))
        self.assertIn(f"Parent isolate CSV SHA-256: {parent_hash}", first)
        self.assertIn("FOR SYSTEM_TIME AS OF @source_as_of_utc", first)
        self.assertIn("i.Run AS sra_run_accessions_raw", first)
        self.assertIn("i.collected_by AS collected_by_raw", first)
        self.assertIn("i.sra_center AS sra_submission_center_raw", first)
        self.assertIn("i.host_disease AS host_disease_raw", first)
        self.assertNotIn("pdbrowser.ast", first)
        self.assertNotIn("clinical_indication", "\n".join(MODULE.OUTPUT_COLUMNS))
        self.assertNotIn(" AS site", first)

    def test_offline_dry_run_needs_no_project_bq_or_credentials(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, mock.patch.dict(
            os.environ, {"GOOGLE_CLOUD_PROJECT": "", "PATH": ""}, clear=False
        ):
            root = Path(temporary)
            isolates, ast, parent_manifest = self.make_parent(root)
            output_dir = root / "enrichment-v1"
            with mock.patch.object(MODULE.shutil, "which") as which:
                code = MODULE.main(
                    [
                        "--dry-run",
                        "--output-dir",
                        str(output_dir),
                        "--parent-isolates",
                        str(isolates),
                        "--parent-ast",
                        str(ast),
                        "--parent-manifest",
                        str(parent_manifest),
                        "--source-as-of-utc",
                        "2026-08-21T03:00:00Z",
                    ]
                )
                which.assert_not_called()
            self.assertEqual(code, 0)
            manifest = json.loads(
                (output_dir / MODULE.MANIFEST_NAME).read_text(encoding="utf-8")
            )
            query = (output_dir / MODULE.QUERY_NAME).read_text(encoding="utf-8")
            self.assertEqual(manifest["status"], "planned")
            self.assertIsNone(manifest["retrieval"])
            self.assertFalse(manifest["selection"]["ast_table_reselected"])
            self.assertEqual(manifest["selection"]["target_count"], 2)
            self.assertEqual(manifest["query"]["sha256"], MODULE.sha256_text(query))
            self.assertEqual(
                manifest["parent_snapshot"]["ast"]["sha256"], MODULE.sha256_file(ast)
            )
            self.assertIn("clinical_indication", manifest["field_semantics"]["intentionally_unpopulated_canonical_fields"])
            self.assertEqual(
                set(manifest["field_semantics"]["fields"]), set(MODULE.OUTPUT_COLUMNS)
            )
            self.assertFalse((output_dir / MODULE.OUTPUT_NAME).exists())
            self.assertFalse((output_dir / MODULE.SOURCE_METADATA_NAME).exists())

    def test_dry_run_refuses_to_overwrite_or_share_parent_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            isolates, ast, parent_manifest = self.make_parent(root)
            common = [
                "--dry-run",
                "--parent-isolates",
                str(isolates),
                "--parent-ast",
                str(ast),
                "--parent-manifest",
                str(parent_manifest),
                "--source-as-of-utc",
                "2026-08-21T03:00:00Z",
            ]
            self.assertEqual(
                MODULE.main([*common, "--output-dir", str(isolates.parent)]), 2
            )
            self.assertEqual(
                MODULE.main(
                    [*common, "--output-dir", str(isolates.parent / "nested-overlay")]
                ),
                2,
            )
            output_dir = root / "enrichment-v1"
            self.assertEqual(MODULE.main([*common, "--output-dir", str(output_dir)]), 0)
            self.assertEqual(MODULE.main([*common, "--output-dir", str(output_dir)]), 2)

    def test_valid_output_records_coverage_and_recovered_assembly(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            context = self.load_context(root)
            output = root / "output.csv"
            self.write_enrichment(
                output,
                [
                    self.enrichment_row(
                        "PDT000000002.1", "SAMN00000002", "GCA_000000002.1"
                    ),
                    self.enrichment_row(
                        "PDT000000001.1", "SAMN00000001", "GCA_000000001.1"
                    ),
                ],
            )
            result = MODULE.validate_enrichment_output(output, context)
            self.assertEqual(result["rows"], 2)
            self.assertEqual(result["reconciliation"]["biosample_exact"], 2)
            self.assertEqual(result["reconciliation"]["assembly"]["same"], 1)
            self.assertEqual(result["reconciliation"]["assembly"]["recovered"], 1)
            self.assertEqual(result["coverage"]["collected_by_raw"]["missing"], 2)
            self.assertEqual(result["coverage"]["sra_run_accessions_raw"]["present"], 2)

    def test_output_identity_and_match_conflicts_fail_closed(self) -> None:
        cases = (
            ("source_match_count", "2", "Source-match conflict"),
            ("source_biosample_accession", "SAMN99999999", "BioSample identity conflict"),
            ("source_assembly_accession", "GCA_000000999.1", "Assembly identity conflict"),
        )
        for field, value, message in cases:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                context = self.load_context(root)
                rows = [
                    self.enrichment_row(
                        "PDT000000001.1", "SAMN00000001", "GCA_000000001.1"
                    ),
                    self.enrichment_row(
                        "PDT000000002.1", "SAMN00000002", "GCA_000000002.1"
                    ),
                ]
                rows[0][field] = value
                output = root / "output.csv"
                self.write_enrichment(output, rows)
                with self.assertRaisesRegex(MODULE.EnrichmentError, message):
                    MODULE.validate_enrichment_output(output, context)

    def test_blank_or_non_ecoli_source_identity_fails_closed(self) -> None:
        cases = (
            ("source_scientific_name", "", "Scientific-name identity conflict"),
            ("source_scientific_name", "Escherichia fergusonii", "Scientific-name identity conflict"),
            ("source_species_taxonomy_id", "", "Species-taxonomy identity conflict"),
            ("source_species_taxonomy_id", "999", "Species-taxonomy identity conflict"),
        )
        for field, value, message in cases:
            with self.subTest(field=field, value=value), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                context = self.load_context(root)
                rows = [
                    self.enrichment_row(
                        "PDT000000001.1", "SAMN00000001", "GCA_000000001.1"
                    ),
                    self.enrichment_row(
                        "PDT000000002.1", "SAMN00000002", "GCA_000000002.1"
                    ),
                ]
                rows[0][field] = value
                output = root / "output.csv"
                self.write_enrichment(output, rows)
                with self.assertRaisesRegex(MODULE.EnrichmentError, message):
                    MODULE.validate_enrichment_output(output, context)

    def test_source_schema_validation_detects_drift(self) -> None:
        metadata = self.source_metadata()
        MODULE.validate_source_metadata(metadata)
        metadata["schema"]["fields"] = [
            field
            for field in metadata["schema"]["fields"]
            if field["name"] != "collected_by"
        ]
        with self.assertRaisesRegex(MODULE.EnrichmentError, "lacks enrichment fields"):
            MODULE.validate_source_metadata(metadata)

    def test_mocked_live_acquisition_installs_verified_atomic_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            context = self.load_context(root)
            query = MODULE.build_query(context["targets"], context["hashes"]["isolates"])
            manifest = MODULE.base_manifest(
                context=context,
                query=query,
                generated_at="2026-08-21T03:00:00Z",
                source_as_of_utc="2026-08-21T03:00:00Z",
                maximum_bytes_billed=MODULE.DEFAULT_MAXIMUM_BYTES_BILLED,
            )
            metadata_json = json.dumps(self.source_metadata())
            fixed_uuid = "00000000000000000000000000000001"
            job_id = f"amr_ncbi_isolate_enrichment_{fixed_uuid}"
            job_json = json.dumps(
                self.job_metadata(
                    query=query,
                    job_id=job_id,
                    project="unbilled-sandbox-project",
                    source_as_of_utc="2026-08-21T03:00:00Z",
                    maximum_bytes_billed=MODULE.DEFAULT_MAXIMUM_BYTES_BILLED,
                )
            )
            captures = [
                subprocess.CompletedProcess([], 0, metadata_json, ""),
                subprocess.CompletedProcess([], 0, "BigQuery CLI 2.1.37\n", ""),
                subprocess.CompletedProcess([], 0, job_json, ""),
                subprocess.CompletedProcess([], 0, metadata_json, ""),
            ]

            def fake_query(command: list[str], destination: Path, sql: str) -> None:
                self.assertEqual(sql, query)
                self.assertFalse(any("SELECT" in argument for argument in command))
                self.assertIn(f"--job_id={job_id}", command)
                self.assertIn(
                    f"--maximum_bytes_billed={MODULE.DEFAULT_MAXIMUM_BYTES_BILLED}",
                    command,
                )
                self.write_enrichment(
                    destination,
                    [
                        self.enrichment_row(
                            "PDT000000001.1", "SAMN00000001", "GCA_000000001.1"
                        ),
                        self.enrichment_row(
                            "PDT000000002.1", "SAMN00000002", "GCA_000000002.1"
                        ),
                    ],
                )

            private_bq = r"C:\Users\private-user\GoogleCloudSDK\bin\bq.CMD"
            output_dir = root / "enrichment-v1"
            with mock.patch.object(MODULE.shutil, "which", return_value=private_bq), mock.patch.object(
                MODULE.uuid, "uuid4", return_value=MODULE.uuid.UUID(hex=fixed_uuid)
            ), mock.patch.object(MODULE, "run_capture", side_effect=captures), mock.patch.object(
                MODULE, "run_query_to_file", side_effect=fake_query
            ):
                result = MODULE.acquire(
                    output_dir=output_dir,
                    project="unbilled-sandbox-project",
                    source_as_of_utc="2026-08-21T03:00:00Z",
                    query=query,
                    manifest=manifest,
                    context=context,
                    maximum_bytes_billed=MODULE.DEFAULT_MAXIMUM_BYTES_BILLED,
                )

            self.assertEqual(result["status"], "complete")
            self.assertEqual(result["retrieval"]["client"]["basename"], "bq.CMD")
            self.assertNotIn("private-user", json.dumps(result))
            self.assertEqual(result["retrieval"]["output"]["rows"], 2)
            output_path = output_dir / MODULE.OUTPUT_NAME
            metadata_path = output_dir / MODULE.SOURCE_METADATA_NAME
            job_metadata_path = output_dir / MODULE.JOB_METADATA_NAME
            self.assertEqual(
                result["retrieval"]["output"]["sha256"], MODULE.sha256_file(output_path)
            )
            self.assertEqual(
                result["retrieval"]["source_table_metadata"]["sha256"],
                MODULE.sha256_file(metadata_path),
            )
            self.assertEqual(
                result["retrieval"]["job"]["metadata_sha256"],
                MODULE.sha256_file(job_metadata_path),
            )
            self.assertEqual(result["retrieval"]["job"]["total_bytes_processed"], 123456)
            installed_payload = "".join(
                path.read_text(encoding="utf-8") for path in output_dir.iterdir()
            )
            self.assertNotIn("private-account@example.org", installed_payload)
            self.assertNotIn("private-job-link", installed_payload)
            on_disk = json.loads(
                (output_dir / MODULE.MANIFEST_NAME).read_text(encoding="utf-8")
            )
            self.assertEqual(on_disk["query"]["parameters"]["source_as_of_utc"], "2026-08-21T03:00:00Z")
            self.assertEqual(on_disk["reconciliation"]["assembly"]["recovered"], 1)

    def test_failed_live_acquisition_leaves_no_partial_output_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            context = self.load_context(root)
            query = MODULE.build_query(context["targets"], context["hashes"]["isolates"])
            manifest = MODULE.base_manifest(
                context=context,
                query=query,
                generated_at="2026-08-21T03:00:00Z",
                source_as_of_utc="2026-08-21T03:00:00Z",
                maximum_bytes_billed=MODULE.DEFAULT_MAXIMUM_BYTES_BILLED,
            )
            metadata_json = json.dumps(self.source_metadata())
            captures = [
                subprocess.CompletedProcess([], 0, metadata_json, ""),
                subprocess.CompletedProcess([], 0, "BigQuery CLI 2.1.37\n", ""),
            ]
            output_dir = root / "failed-enrichment"
            with mock.patch.object(MODULE.shutil, "which", return_value="bq"), mock.patch.object(
                MODULE, "run_capture", side_effect=captures
            ), mock.patch.object(
                MODULE,
                "run_query_to_file",
                side_effect=MODULE.EnrichmentError("simulated query failure"),
            ):
                with self.assertRaisesRegex(MODULE.EnrichmentError, "simulated query failure"):
                    MODULE.acquire(
                        output_dir=output_dir,
                        project="unbilled-sandbox-project",
                        source_as_of_utc="2026-08-21T03:00:00Z",
                        query=query,
                        manifest=manifest,
                        context=context,
                        maximum_bytes_billed=MODULE.DEFAULT_MAXIMUM_BYTES_BILLED,
                    )
            self.assertFalse(output_dir.exists())

    def test_bq_command_keeps_sql_on_stdin_and_is_sandbox_compatible(self) -> None:
        command = MODULE.bq_query_command(
            bq="bq.cmd",
            project="unbilled-sandbox-project",
            source_as_of_utc="2026-08-21T03:00:00Z",
            job_id="amr_test_job_1",
            maximum_bytes_billed=123456789,
        )
        self.assertIn("--project_id=unbilled-sandbox-project", command)
        self.assertIn("--job_id=amr_test_job_1", command)
        self.assertIn("--maximum_bytes_billed=123456789", command)
        self.assertIn(
            "--parameter=source_as_of_utc:TIMESTAMP:2026-08-21T03:00:00Z",
            command,
        )
        self.assertFalse(any("SELECT" in argument for argument in command))
        with tempfile.TemporaryDirectory() as temporary, mock.patch.object(
            MODULE.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "", "")
        ) as run:
            destination = Path(temporary) / "output.csv"
            MODULE.run_query_to_file(command, destination, "SELECT 1\n")
            self.assertEqual(run.call_args.kwargs["input"], "SELECT 1\n")
            self.assertNotIn("SELECT 1", run.call_args.args[0])

    def test_client_identity_never_records_private_directories(self) -> None:
        self.assertEqual(
            MODULE.executable_identity(r"C:\Users\private-user\GoogleCloudSDK\bin\bq.CMD"),
            "bq.CMD",
        )
        self.assertEqual(MODULE.executable_identity("/opt/google-cloud-sdk/bin/bq"), "bq")


if __name__ == "__main__":
    unittest.main()
