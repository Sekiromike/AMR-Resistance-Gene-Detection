from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).parents[1]
MODULE_PATH = ROOT / "scripts" / "normalize_ncbi_ast.py"
SPEC = importlib.util.spec_from_file_location("normalize_ncbi_ast", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


LIVE_AST_COLUMNS = (
    "id",
    "biosample_acc",
    "taxgroup_name",
    "scientific_name",
    "epi_type",
    "isolation_source",
    "geo_loc_name",
    "target_acc",
    "antibiotic",
    "phenotype",
    "measurement_sign",
    "mic",
    "mic_secondary",
    "disk_diffusion",
    "disk_diffusion_secondary",
    "standard",
    "reagent",
    "platform",
    "vendor",
    "host",
    "collection_date",
    "creation_date",
    "bioproject_acc",
    "checksum",
)

LIVE_ISOLATE_COLUMNS = (
    "target_acc",
    "biosample_acc",
    "asm_acc",
    "scientific_name",
    "species_taxid",
    "taxid",
    "strain",
    "bioproject_acc",
    "collection_date",
    "geo_loc_name",
    "host",
    "isolation_source",
    "source_type",
    "target_creation_date",
    "erd_group",
    "minsame",
    "mindiff",
)


def live_ast_frame(**overrides: object) -> pd.DataFrame:
    row: dict[str, object] = {
        "id": "008_PDT000000001.1",
        "biosample_acc": "SAMN00000001",
        "taxgroup_name": "E.coli and Shigella",
        "scientific_name": "Escherichia coli",
        "epi_type": "clinical",
        "isolation_source": "blood",
        "geo_loc_name": "USA: Massachusetts",
        "target_acc": "PDT000000001.1",
        "antibiotic": "ciprofloxacin",
        "phenotype": "resistant",
        "measurement_sign": "==",
        "mic": "4.0",
        "mic_secondary": None,
        "disk_diffusion": None,
        "disk_diffusion_secondary": None,
        "standard": "CLSI",
        "reagent": "GM-NEG",
        "platform": "Phoenix",
        "vendor": "Becton Dickinson",
        "host": "Homo sapiens",
        "collection_date": "2024-02-03",
        "creation_date": "2024-02-10T12:30:00Z",
        "bioproject_acc": "PRJNA000001",
        "checksum": "0123456789abcdef0123456789abcdef",
    }
    unknown = set(overrides).difference(row)
    if unknown:
        raise AssertionError(f"Unknown live AST fixture fields: {sorted(unknown)}")
    row.update(overrides)
    return pd.DataFrame([row], columns=LIVE_AST_COLUMNS, dtype="string")


def live_isolates_frame(**overrides: object) -> pd.DataFrame:
    row: dict[str, object] = {
        "target_acc": "PDT000000001.1",
        "biosample_acc": "SAMN00000001",
        "asm_acc": "GCA_000000001.1",
        "scientific_name": "Escherichia coli",
        "species_taxid": "562",
        "taxid": "562",
        "strain": "EC-1",
        "bioproject_acc": "PRJNA000001",
        "collection_date": "2024-02-03",
        "geo_loc_name": "USA: Massachusetts",
        "host": "Homo sapiens",
        "isolation_source": "blood",
        "source_type": "human",
        "target_creation_date": "2024-02-10T16:30:00Z",
        "erd_group": "PDS000000001.1",
        "minsame": "0",
        "mindiff": "25",
    }
    unknown = set(overrides).difference(row)
    if unknown:
        raise AssertionError(f"Unknown live isolate fixture fields: {sorted(unknown)}")
    row.update(overrides)
    return pd.DataFrame([row], columns=LIVE_ISOLATE_COLUMNS, dtype="string")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_acquisition_bundle(root: Path) -> tuple[Path, Path, Path]:
    ast_path = root / "ncbi_ast.csv"
    isolates_path = root / "ncbi_isolates.csv"
    ast_query_path = root / "ncbi_ast_query.sql"
    isolates_query_path = root / "ncbi_isolates_query.sql"
    metadata_path = root / "ncbi_source_tables.json"
    manifest_path = root / MODULE.DEFAULT_ACQUISITION_MANIFEST
    live_ast_frame().to_csv(ast_path, index=False)
    live_isolates_frame().to_csv(isolates_path, index=False)
    ast_query_path.write_text("SELECT * FROM ast\n", encoding="utf-8")
    isolates_query_path.write_text("SELECT * FROM isolates\n", encoding="utf-8")
    metadata_path.write_text('{"tables": {}}\n', encoding="utf-8")
    manifest = {
        "status": "complete",
        "queries": {
            "ast": {"file": ast_query_path.name, "sha256": sha256_file(ast_query_path)},
            "isolates": {
                "file": isolates_query_path.name,
                "sha256": sha256_file(isolates_query_path),
            },
        },
        "retrieval": {
            "outputs": {
                "ast": {"rows": 1, "sha256": sha256_file(ast_path)},
                "isolates": {"rows": 1, "sha256": sha256_file(isolates_path)},
            },
            "source_table_metadata": {
                "file": metadata_path.name,
                "sha256": sha256_file(metadata_path),
            },
        },
    }
    manifest_path.write_text(json.dumps(manifest) + "\n", encoding="utf-8")
    return ast_path, isolates_path, manifest_path


class NcbiAstNormalizationTests(unittest.TestCase):
    def test_exact_live_header_preserves_native_provenance_and_source_geography(self) -> None:
        ast = live_ast_frame()
        self.assertEqual(tuple(ast.columns), LIVE_AST_COLUMNS)
        self.assertEqual(len(ast.columns), 24)

        included, excluded = MODULE.normalize_ast(ast, live_isolates_frame())

        self.assertTrue(excluded.empty)
        self.assertEqual(len(included), 1)
        row = included.iloc[0]
        self.assertEqual(row["source_ast_record_id"], "008_PDT000000001.1")
        self.assertEqual(row["source_record_checksum"], "0123456789abcdef0123456789abcdef")
        self.assertEqual(row["taxgroup_name"], "E.coli and Shigella")
        self.assertEqual(row["epi_type"], "clinical")
        self.assertEqual(row["geo_loc_name"], "USA: Massachusetts")
        self.assertTrue(pd.isna(row["country"]))
        self.assertEqual(row["ast_record_creation_date"], "2024-02-10T12:30:00Z")
        self.assertEqual(row["target_creation_date"], "2024-02-10T16:30:00Z")
        self.assertEqual(row["submitted_ast_category"], "R")
        self.assertEqual(row["raw_submitted_ast_category"], "resistant")
        self.assertEqual(row["ast_measurement_type"], "MIC")
        self.assertEqual(row["ast_value"], 4.0)
        self.assertEqual(row["ast_unit"], "mg/L")
        self.assertEqual(row["measurement_sign"], "=")
        self.assertEqual(row["raw_measurement_sign"], "==")

    def test_explicit_country_is_preserved_without_overwriting_geo_loc_name(self) -> None:
        ast = live_ast_frame()
        ast["country"] = pd.Series(["United States"], dtype="string")

        included, excluded = MODULE.normalize_ast(ast)

        self.assertTrue(excluded.empty)
        self.assertEqual(included.loc[0, "country"], "United States")
        self.assertEqual(included.loc[0, "geo_loc_name"], "USA: Massachusetts")

    def test_missing_breakpoint_metadata_is_flagged_without_structural_exclusion(self) -> None:
        included, excluded = MODULE.normalize_ast(live_ast_frame())

        self.assertTrue(excluded.empty)
        self.assertEqual(len(included), 1)
        self.assertFalse(included.loc[0, "breakpoint_eligible"])
        self.assertEqual(
            included.loc[0, "breakpoint_ineligibility_reasons"],
            "missing_ast_method;missing_ast_testing_date",
        )

    def test_zone_requires_method_disk_content_and_testing_date_for_breakpoints(self) -> None:
        ast = live_ast_frame(mic=None, disk_diffusion="24.0")

        included, excluded = MODULE.normalize_ast(ast)

        self.assertTrue(excluded.empty)
        self.assertEqual(included.loc[0, "ast_measurement_type"], "ZONE")
        self.assertEqual(included.loc[0, "ast_unit"], "mm")
        self.assertFalse(included.loc[0, "breakpoint_eligible"])
        self.assertEqual(
            included.loc[0, "breakpoint_ineligibility_reasons"],
            "missing_ast_method;missing_disk_content;missing_ast_testing_date",
        )

    def test_complete_mic_metadata_is_breakpoint_eligible(self) -> None:
        ast = live_ast_frame()
        ast["ast_method"] = pd.Series(["broth microdilution"], dtype="string")
        ast["ast_testing_date"] = pd.Series(["2024-02-04"], dtype="string")

        included, excluded = MODULE.normalize_ast(ast)

        self.assertTrue(excluded.empty)
        self.assertTrue(included.loc[0, "breakpoint_eligible"])
        reason = included.loc[0, "breakpoint_ineligibility_reasons"]
        self.assertTrue(pd.isna(reason) or str(reason) == "")

    def test_malformed_assay_metadata_never_passes_breakpoint_gate(self) -> None:
        cases = (
            (
                {"ast_method": "unknown-proxy", "ast_testing_date": "2024-02-04"},
                "unsupported_ast_method_or_measurement_mismatch",
            ),
            (
                {"ast_method": "broth microdilution", "ast_testing_date": "not-a-date"},
                "invalid_ast_testing_date",
            ),
            (
                {
                    "mic": None,
                    "disk_diffusion": "24",
                    "ast_method": "disk diffusion",
                    "ast_testing_date": "2024-02-04",
                    "disk_content": "garbage",
                },
                "invalid_disk_content",
            ),
        )
        for additions, expected_reason in cases:
            with self.subTest(expected_reason=expected_reason):
                ast = live_ast_frame(
                    **{key: value for key, value in additions.items() if key in LIVE_AST_COLUMNS}
                )
                for key, value in additions.items():
                    if key not in LIVE_AST_COLUMNS:
                        ast[key] = pd.Series([value], dtype="string")
                included, excluded = MODULE.normalize_ast(ast)
                self.assertTrue(excluded.empty)
                self.assertFalse(included.loc[0, "breakpoint_eligible"])
                self.assertIn(
                    expected_reason,
                    included.loc[0, "breakpoint_ineligibility_reasons"].split(";"),
                )

    def test_rejects_non_ecoli_and_missing_measurements_as_structural_failures(self) -> None:
        non_ecoli = live_ast_frame(scientific_name="Klebsiella pneumoniae")
        missing_measurement = live_ast_frame(
            id="009_PDT000000002.1",
            target_acc="PDT000000002.1",
            biosample_acc="SAMN00000002",
            checksum="fedcba9876543210fedcba9876543210",
        )
        missing_measurement["mic"] = pd.NA
        ast = pd.DataFrame(
            [non_ecoli.iloc[0].to_dict(), missing_measurement.iloc[0].to_dict()],
            columns=LIVE_AST_COLUMNS,
            dtype="string",
        )

        included, excluded = MODULE.normalize_ast(ast)

        self.assertTrue(included.empty)
        self.assertEqual(
            set(excluded["exclusion_reason"]),
            {"not_escherichia_coli", "missing_quantitative_ast"},
        )

    def test_missing_comparator_is_not_invented(self) -> None:
        included, excluded = MODULE.normalize_ast(live_ast_frame(measurement_sign=None))

        self.assertTrue(included.empty)
        self.assertEqual(excluded.loc[0, "exclusion_reason"], "missing_measurement_sign")

    def test_blank_required_identifiers_are_structural_exclusions(self) -> None:
        cases = (
            ("antibiotic", "missing_antibiotic"),
            ("target_acc", "missing_isolate_id"),
            ("biosample_acc", "missing_biosample_accession"),
        )
        for field, expected_reason in cases:
            with self.subTest(field=field):
                included, excluded = MODULE.normalize_ast(live_ast_frame(**{field: "   "}))
                self.assertTrue(included.empty)
                self.assertEqual(len(excluded), 1)
                self.assertEqual(excluded.loc[0, "exclusion_reason"], expected_reason)

    def test_aliases_coalesce_blanks_but_reject_conflicting_values(self) -> None:
        ast = live_ast_frame()
        ast["isolate_id"] = pd.Series(["   "], dtype="string")
        included, excluded = MODULE.normalize_ast(ast)
        self.assertTrue(excluded.empty)
        self.assertEqual(included.loc[0, "isolate_id"], "PDT000000001.1")

        ast["isolate_id"] = pd.Series(["PDT999999999.1"], dtype="string")
        with self.assertRaisesRegex(MODULE.NormalizationError, r"(?i)conflicting aliases"):
            MODULE.normalize_ast(ast)

    def test_isolate_export_supplies_assembly_without_replacing_identity(self) -> None:
        included, excluded = MODULE.normalize_ast(live_ast_frame(), live_isolates_frame())

        self.assertTrue(excluded.empty)
        self.assertEqual(included.loc[0, "assembly_accession"], "GCA_000000001.1")
        self.assertEqual(included.loc[0, "biosample_accession"], "SAMN00000001")

    def test_conflicting_duplicate_isolate_rows_raise_normalization_error(self) -> None:
        isolates = pd.concat(
            [
                live_isolates_frame(),
                live_isolates_frame(asm_acc="GCA_000000002.1"),
            ],
            ignore_index=True,
        )

        with self.assertRaisesRegex(
            MODULE.NormalizationError,
            r"(?i)(conflict|ambiguous).*(isolate|metadata)|(?:isolate|metadata).*(conflict|ambiguous)",
        ):
            MODULE.normalize_ast(live_ast_frame(), isolates)

    def test_ast_and_isolate_biosample_conflict_raises_normalization_error(self) -> None:
        isolates = live_isolates_frame(biosample_acc="SAMN99999999")

        with self.assertRaisesRegex(MODULE.NormalizationError, r"(?i)biosample"):
            MODULE.normalize_ast(live_ast_frame(), isolates)

    def test_ast_and_isolate_assembly_conflict_raises_normalization_error(self) -> None:
        ast = live_ast_frame()
        ast["assembly_accession"] = pd.Series(["GCA_000000009.1"], dtype="string")

        with self.assertRaisesRegex(MODULE.NormalizationError, r"(?i)assembly"):
            MODULE.normalize_ast(ast, live_isolates_frame())

    def test_unmatched_or_wrong_species_isolate_metadata_fails_closed(self) -> None:
        with self.assertRaisesRegex(MODULE.NormalizationError, r"(?i)no exact.*match"):
            MODULE.normalize_ast(
                live_ast_frame(),
                live_isolates_frame(target_acc="PDT999999999.1"),
            )
        with self.assertRaisesRegex(MODULE.NormalizationError, r"(?i)scientific_name"):
            MODULE.normalize_ast(
                live_ast_frame(),
                live_isolates_frame(scientific_name="Klebsiella pneumoniae"),
            )

    def test_cli_manifest_authenticates_inputs_and_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ast_path = root / "ncbi_ast.csv"
            isolates_path = root / "ncbi_isolates.csv"
            normalized_path = root / "normalized.csv"
            exclusions_path = root / "exclusions.csv"
            manifest_path = root / "manifest.json"
            live_ast_frame().to_csv(ast_path, index=False)
            live_isolates_frame().to_csv(isolates_path, index=False)

            result = subprocess.run(
                [
                    sys.executable,
                    str(MODULE_PATH),
                    "--ast",
                    str(ast_path),
                    "--isolates",
                    str(isolates_path),
                    "--output",
                    str(normalized_path),
                    "--exclusions",
                    str(exclusions_path),
                    "--manifest",
                    str(manifest_path),
                    "--allow-unverified-input",
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertRegex(str(manifest["mapping_version"]), r"^\d+\.\d+\.\d+$")
            self.assertEqual(
                manifest["scientific_status"],
                "SMOKE_TEST_ONLY_UNVERIFIED_INPUT",
            )
            self.assertEqual(
                manifest["normalizer"]["script_sha256"],
                sha256_file(MODULE_PATH),
            )
            self.assertEqual(manifest["source"]["columns"], list(LIVE_AST_COLUMNS))
            self.assertEqual(
                manifest["isolates_source"]["columns"],
                list(LIVE_ISOLATE_COLUMNS),
            )

            for entry, path in (
                (manifest["source"], ast_path),
                (manifest["isolates_source"], isolates_path),
                (manifest["outputs"]["normalized"], normalized_path),
                (manifest["outputs"]["exclusions"], exclusions_path),
            ):
                self.assertEqual(entry["sha256"], sha256_file(path))
                self.assertEqual(entry["bytes"], path.stat().st_size)

    def test_cli_requires_acquisition_manifest_without_explicit_test_override(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ast_path = root / "ncbi_ast.csv"
            live_ast_frame().to_csv(ast_path, index=False)
            result = subprocess.run(
                [sys.executable, str(MODULE_PATH), "--ast", str(ast_path)],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("No completed acquisition manifest", result.stderr)

    def test_acquisition_provenance_verifies_every_referenced_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ast_path, isolates_path, manifest_path = write_acquisition_bundle(root)
            manifest_hash = sha256_file(manifest_path)
            amendment_path = root / MODULE.DEFAULT_ACQUISITION_AMENDMENT
            amendment_path.write_text(
                json.dumps(
                    {
                        "amendment_status": "complete",
                        "original_manifest": {"sha256": manifest_hash},
                        "unchanged_artifacts": {
                            "ast_sha256": sha256_file(ast_path),
                            "isolates_sha256": sha256_file(isolates_path),
                            "source_table_metadata_sha256": sha256_file(
                                root / "ncbi_source_tables.json"
                            ),
                        },
                    }
                )
                + "\n",
                encoding="utf-8",
            )

            provenance = MODULE.verify_acquisition_provenance(
                manifest_path,
                ast_path,
                isolates_path,
                ast_rows=1,
                isolate_rows=1,
                amendment_path=amendment_path,
            )

            self.assertEqual(provenance["sha256"], manifest_hash)
            self.assertEqual(
                set(provenance["verified_artifacts"]),
                {
                    "output_ast",
                    "output_isolates",
                    "query_ast",
                    "query_isolates",
                    "source_table_metadata",
                },
            )
            self.assertEqual(provenance["amendment"]["sha256"], sha256_file(amendment_path))

    def test_acquisition_provenance_rejects_tampered_source(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ast_path, isolates_path, manifest_path = write_acquisition_bundle(root)
            ast_path.write_text("tampered\n", encoding="utf-8")

            with self.assertRaisesRegex(MODULE.NormalizationError, r"(?i)hash mismatch"):
                MODULE.verify_acquisition_provenance(
                    manifest_path,
                    ast_path,
                    isolates_path,
                    ast_rows=1,
                    isolate_rows=1,
                )


if __name__ == "__main__":
    unittest.main()


class CountryDerivationTests(unittest.TestCase):
    """Amendment 006: country from the INSDC geo_loc_name prefix, never guessed."""

    VOCAB = {"USA", "Canada", "United Kingdom", "Japan", "USSR"}

    def _frame(self, **columns):
        n = len(next(iter(columns.values())))
        base = {"country": [pd.NA] * n, "geo_loc_name": [pd.NA] * n}
        base.update(columns)
        return pd.DataFrame(base, dtype="string")

    def test_prefix_in_vocabulary_fills_country(self) -> None:
        out = MODULE.derive_country(
            self._frame(geo_loc_name=["USA: California", "United Kingdom:London", "Canada"]), self.VOCAB
        )
        self.assertEqual(out["country"].tolist(), ["USA", "United Kingdom", "Canada"])
        self.assertTrue((out["country_source"] == "insdc_geo_loc_name_prefix").all())

    def test_missingness_token_and_unknown_prefix_stay_absent(self) -> None:
        out = MODULE.derive_country(
            self._frame(geo_loc_name=["not collected", "Atlantis: somewhere", pd.NA]), self.VOCAB
        )
        self.assertTrue(out["country"].isna().all())
        self.assertTrue((out["country_source"] == "absent").all())

    def test_explicit_country_always_wins(self) -> None:
        out = MODULE.derive_country(self._frame(country=["Canada"], geo_loc_name=["USA: Texas"]), self.VOCAB)
        self.assertEqual(out.loc[0, "country"], "Canada")
        self.assertEqual(out.loc[0, "country_source"], "explicit")

    def test_historical_names_are_valid_vocabulary(self) -> None:
        out = MODULE.derive_country(self._frame(geo_loc_name=["USSR: Moscow"]), self.VOCAB)
        self.assertEqual(out.loc[0, "country"], "USSR")

    def test_without_a_vocabulary_nothing_is_derived(self) -> None:
        out = MODULE.derive_country(self._frame(geo_loc_name=["USA: California"]), None)
        self.assertTrue(out["country"].isna().all())

    def test_frozen_vocabulary_loads_and_rejects_bad_files(self) -> None:
        names = MODULE.load_country_vocabulary(Path(__file__).parents[1] / "config" / "insdc_geo_loc_name_vocabulary.tsv")
        self.assertIn("USA", names)
        self.assertIn("Czechoslovakia", names)  # historical
        self.assertEqual(len(names), 297)
        with tempfile.TemporaryDirectory() as temporary:
            bad = Path(temporary) / "bad.tsv"
            bad.write_text("name\tstatus\nUSA\tcurrent\nUSA\tcurrent\n", encoding="utf-8")
            with self.assertRaises(MODULE.NormalizationError):
                MODULE.load_country_vocabulary(bad)
