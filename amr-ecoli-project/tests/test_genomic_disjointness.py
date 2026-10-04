from __future__ import annotations

import csv
import hashlib
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
SPEC = importlib.util.spec_from_file_location(
    "build_genomic_disjointness", SCRIPTS / "build_genomic_disjointness.py"
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_result(root: Path, row: dict[str, str]) -> Path:
    result = root / row["cohort"] / row["isolate_id"]
    result.mkdir(parents=True)
    assembly = b">contig\nACGT\n"
    provenance = {
        "completed_at_utc": "2026-08-28T00:00:00Z",
        "isolate_id": row["isolate_id"],
        "assembly_sha256": hashlib.sha256(assembly).hexdigest(),
    }
    if row["source_type"] == "raw_reads":
        provenance.update(
            {
                "cohort": row["cohort"],
                "run_accession": row["sequence_accession"],
                "expected_sra_md5": row["expected_read_md5"],
                "observed_sra_md5": row["expected_read_md5"],
            }
        )
    else:
        provenance["assembly_accession"] = row["sequence_accession"]
    files = {
        "assembly.fna": assembly,
        "provenance.txt": "".join(f"{key}={value}\n" for key, value in provenance.items()).encode(),
    }
    if row["source_type"] == "raw_reads":
        files.update({"read-files.sha256": b"reads preserved by hash\n", "vdb-validate.txt": b"ok\n"})
    for name, content in files.items():
        (result / name).write_bytes(content)
    (result / "artifact.sha256").write_text(
        "".join(f"{hashlib.sha256(content).hexdigest()}  {name}\n" for name, content in files.items()),
        encoding="ascii",
    )
    (result / "COMPLETE").touch()
    return result


class PlanTests(unittest.TestCase):
    def test_plan_uses_assembly_or_one_run_and_never_needs_phenotypes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            development = root / "development.csv"
            external = root / "external.csv"
            write_csv(
                development,
                [
                    {"target_acc": "PDT1", "source_biosample_accession": "SAMN1", "source_assembly_accession": "GCF_000000001.1", "source_species_taxonomy_id": "562", "sra_run_accessions_raw": ""},
                    {"target_acc": "PDT2", "source_biosample_accession": "SAMN2", "source_assembly_accession": "", "source_species_taxonomy_id": "562", "sra_run_accessions_raw": "SRR2"},
                ],
            )
            write_csv(
                external,
                [{"isolate_id": "JB1", "supplement_species": "Escherichia coli", "biosample_accession": "SAMD1", "primary_short_read_run": "DRR1", "primary_short_read_md5": "a" * 32, "sequence_link_status": "linked_one_short_read_run"}],
            )
            rows = MODULE.build_sequence_requests(
                development, external, expected_development=2, expected_external=1
            )
            self.assertEqual([row["source_type"] for row in rows], ["registered_assembly", "raw_reads", "raw_reads"])
            self.assertEqual(rows[-1]["sequence_id"], "external::JB1")

    def test_plan_rejects_ambiguous_development_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            development = root / "development.csv"
            external = root / "external.csv"
            write_csv(development, [{"target_acc": "PDT1", "source_biosample_accession": "SAMN1", "source_assembly_accession": "", "source_species_taxonomy_id": "562", "sra_run_accessions_raw": "SRR1;SRR2"}])
            write_csv(external, [{"isolate_id": "JB1", "supplement_species": "Escherichia coli", "biosample_accession": "SAMD1", "primary_short_read_run": "DRR1", "primary_short_read_md5": "b" * 32, "sequence_link_status": "linked_one_short_read_run"}])
            with self.assertRaisesRegex(MODULE.DisjointnessError, "has 2 SRA runs"):
                MODULE.build_sequence_requests(development, external, expected_development=1, expected_external=1)

    def test_written_plan_is_lf_only_and_freezes_broad_search_floor(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            development = root / "development.csv"
            external = root / "external.csv"
            write_csv(development, [{"target_acc": "PDT1", "source_biosample_accession": "SAMN1", "source_assembly_accession": "GCF_000000001.1", "source_species_taxonomy_id": "562", "sra_run_accessions_raw": ""}])
            write_csv(external, [{"isolate_id": "JB1", "supplement_species": "Escherichia coli", "biosample_accession": "SAMD1", "primary_short_read_run": "DRR1", "primary_short_read_md5": "c" * 32, "sequence_link_status": "linked_one_short_read_run"}])
            requests = root / "requests.tsv"
            payload = MODULE.write_plan(
                development,
                external,
                requests,
                root / "assemblies.txt",
                root / "runs.txt",
                root / "contract.json",
                expected_development=1,
                expected_external=1,
                search_min_ani=95.0,
                search_min_af=0.5,
                near_ani=99.9,
                near_af=0.9,
                duplicate_ani=99.99,
                duplicate_af=0.95,
                skani_version="0.3.1",
            )
            self.assertNotIn(b"\r\n", requests.read_bytes())
            self.assertEqual(payload["comparison"]["search_min_reciprocal_aligned_fraction"], 0.5)


class EvaluationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.requests = [
            {"sequence_id": "development::D1", "cohort": "development"},
            {"sequence_id": "development::D2", "cohort": "development"},
            {"sequence_id": "external::E1", "cohort": "external"},
        ]

    def test_duplicate_and_near_neighbor_classes_use_minimum_reciprocal_af(self) -> None:
        comparisons = [
            {"development_sequence_id": "development::D1", "external_sequence_id": "external::E1", "ani_percent": "99.995", "aligned_fraction_development": "0.97", "aligned_fraction_external": "0.96"},
            {"development_sequence_id": "development::D2", "external_sequence_id": "external::E1", "ani_percent": "99.95", "aligned_fraction_development": "0.94", "aligned_fraction_external": "0.91"},
            {"development_sequence_id": "development::D2", "external_sequence_id": "external::E1", "ani_percent": "99.95", "aligned_fraction_development": "0.94", "aligned_fraction_external": "0.89"},
        ]
        rows = MODULE.classify_comparisons(comparisons, self.requests, near_ani=99.9, near_af=0.90, duplicate_ani=99.99, duplicate_af=0.95)
        self.assertEqual([row["collision_class"] for row in rows], ["PUTATIVE_GENOMIC_DUPLICATE", "NEAR_NEIGHBOR"])
        self.assertEqual(len({row["component_id"] for row in rows}), 1)

    def test_evaluate_fails_closed_when_external_query_is_not_completed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            requests = root / "requests.tsv"
            comparisons = root / "comparisons.tsv"
            completed = root / "completed.txt"
            contract = root / "contract.json"
            write_tsv(requests, [{**row, "isolate_id": "", "biosample_accession": "", "source_type": "", "sequence_accession": "", "expected_read_md5": ""} for row in self.requests], MODULE.REQUEST_FIELDS)
            write_tsv(comparisons, [], ["development_sequence_id", "external_sequence_id", "ani_percent", "aligned_fraction_development", "aligned_fraction_external"])
            completed.write_text("", encoding="utf-8")
            contract.write_text(json.dumps({"comparison": {"near_neighbor_min_ani_percent": 99.9, "near_neighbor_min_reciprocal_aligned_fraction": 0.9, "putative_duplicate_min_ani_percent": 99.99, "putative_duplicate_min_reciprocal_aligned_fraction": 0.95}}), encoding="utf-8")
            with self.assertRaisesRegex(MODULE.DisjointnessError, "completion mismatch"):
                MODULE.evaluate(requests, contract, comparisons, completed, root / "collisions.tsv", root / "manifest.json")


class AssemblyAuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self.requests = [
            {
                "sequence_id": "development::D1",
                "cohort": "development",
                "isolate_id": "D1",
                "biosample_accession": "SAMN1",
                "source_type": "registered_assembly",
                "sequence_accession": "GCA_1.1",
                "expected_read_md5": "",
            },
            {
                "sequence_id": "external::E1",
                "cohort": "external",
                "isolate_id": "E1",
                "biosample_accession": "SAMD1",
                "source_type": "raw_reads",
                "sequence_accession": "DRR1",
                "expected_read_md5": "a" * 32,
            },
        ]

    def test_audit_verifies_every_declared_artifact_and_retains_partial_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            requests = root / "requests.tsv"
            assemblies = root / "assemblies"
            write_tsv(requests, self.requests, MODULE.REQUEST_FIELDS)
            for row in self.requests:
                write_result(assemblies, row)
            (assemblies / "external" / "E1.partial-old-job").mkdir()
            payload = MODULE.audit_assemblies(
                requests, assemblies, root / "audit.json", workers=2
            )
            self.assertEqual(payload["status"], "complete_and_verified_with_retained_partials")
            self.assertEqual(payload["counts"]["checksum_verified"], 2)
            self.assertEqual(payload["counts"]["retained_partial_directories"], 1)

    def test_audit_fails_closed_for_corruption_and_missing_completion(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            requests = root / "requests.tsv"
            assemblies = root / "assemblies"
            write_tsv(requests, self.requests, MODULE.REQUEST_FIELDS)
            result = write_result(assemblies, self.requests[0])
            (result / "assembly.fna").write_text(">contig\nCORRUPTED\n", encoding="ascii")
            payload = MODULE.audit_assemblies(
                requests, assemblies, root / "audit.json", workers=1
            )
            self.assertEqual(payload["status"], "incomplete_or_invalid")
            self.assertEqual(payload["counts"]["failed_or_incomplete"], 2)
            issues = {row["sequence_id"]: row["issues"] for row in payload["failures"]}
            self.assertIn("CHECKSUM_MISMATCH:assembly.fna", issues["development::D1"])
            self.assertEqual(issues["external::E1"], ["MISSING_COMPLETE"])

    def test_audit_rejects_valid_hashes_with_wrong_frozen_accession(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            requests = root / "requests.tsv"
            assemblies = root / "assemblies"
            write_tsv(requests, self.requests[:1], MODULE.REQUEST_FIELDS)
            result = write_result(assemblies, self.requests[0])
            provenance = (result / "provenance.txt").read_text(encoding="utf-8")
            provenance = provenance.replace("assembly_accession=GCA_1.1", "assembly_accession=GCA_9.9")
            (result / "provenance.txt").write_text(provenance, encoding="utf-8")
            lines = []
            for name in ("assembly.fna", "provenance.txt"):
                lines.append(f"{hashlib.sha256((result / name).read_bytes()).hexdigest()}  {name}\n")
            (result / "artifact.sha256").write_text("".join(lines), encoding="ascii")
            payload = MODULE.audit_assemblies(
                requests, assemblies, root / "audit.json", workers=1
            )
            self.assertIn("PROVENANCE_MISMATCH:assembly_accession", payload["failures"][0]["issues"])


class ComparisonPreparationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.requests = [
            {
                "sequence_id": "development::D1",
                "cohort": "development",
                "isolate_id": "D1",
                "biosample_accession": "SAMN1",
                "source_type": "registered_assembly",
                "sequence_accession": "GCA_1.1",
                "expected_read_md5": "",
            },
            {
                "sequence_id": "external::E1",
                "cohort": "external",
                "isolate_id": "E1",
                "biosample_accession": "SAMD1",
                "source_type": "raw_reads",
                "sequence_accession": "DRR1",
                "expected_read_md5": "a" * 32,
            },
        ]

    def _fixture(self, root: Path, status: str = "complete_and_verified") -> tuple[Path, Path, Path]:
        requests = root / "requests.tsv"
        assemblies = root / "assemblies"
        audit = root / "audit.json"
        write_tsv(requests, self.requests, MODULE.REQUEST_FIELDS)
        for row in self.requests:
            write_result(assemblies, row)
        audit.write_text(
            json.dumps(
                {
                    "status": status,
                    "inputs": {
                        "requests": {"sha256": MODULE.sha256_file(requests)},
                        "assemblies_root": str(assemblies),
                    },
                    "counts": {
                        "requests": 2,
                        "checksum_verified": 2,
                        "failed_or_incomplete": 0,
                        "extra_complete": 0,
                    },
                }
            ),
            encoding="utf-8",
        )
        return requests, assemblies, audit

    def test_prepare_requires_verified_audit_and_writes_complete_frozen_lists(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            requests, assemblies, audit = self._fixture(root)
            payload = MODULE.prepare_comparison_inputs(
                requests,
                audit,
                assemblies,
                root / "development.txt",
                root / "external.txt",
                root / "paths.tsv",
                root / "completed.txt",
                root / "manifest.json",
            )
            self.assertEqual(payload["counts"]["total_assemblies"], 2)
            self.assertEqual(payload["counts"]["reauthenticated_assemblies"], 2)
            self.assertEqual((root / "completed.txt").read_text(), "external::E1\n")
            self.assertIn("development/D1/assembly.fna", (root / "development.txt").read_text())

    def test_prepare_rejects_incomplete_audit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            requests, assemblies, audit = self._fixture(root, status="incomplete_or_invalid")
            with self.assertRaisesRegex(MODULE.DisjointnessError, "not complete and verified"):
                MODULE.prepare_comparison_inputs(
                    requests,
                    audit,
                    assemblies,
                    root / "development.txt",
                    root / "external.txt",
                    root / "paths.tsv",
                    root / "completed.txt",
                    root / "manifest.json",
                )

    def test_prepare_rejects_assembly_changed_after_successful_audit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            requests, assemblies, audit = self._fixture(root)
            changed = assemblies / "external" / "E1" / "assembly.fna"
            changed.write_text(">contig\nCHANGED\n", encoding="ascii")
            with self.assertRaisesRegex(MODULE.DisjointnessError, "integrity changed after audit"):
                MODULE.prepare_comparison_inputs(
                    requests,
                    audit,
                    assemblies,
                    root / "development.txt",
                    root / "external.txt",
                    root / "paths.tsv",
                    root / "completed.txt",
                    root / "manifest.json",
                    workers=2,
                )

    def test_normalize_maps_paths_and_converts_percent_af_to_fractions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, assemblies, _ = self._fixture(root)
            path_map = root / "paths.tsv"
            write_tsv(
                path_map,
                [
                    {
                        "sequence_id": "development::D1",
                        "cohort": "development",
                        "isolate_id": "D1",
                        "assembly_path": (assemblies / "development" / "D1" / "assembly.fna").as_posix(),
                    },
                    {
                        "sequence_id": "external::E1",
                        "cohort": "external",
                        "isolate_id": "E1",
                        "assembly_path": (assemblies / "external" / "E1" / "assembly.fna").as_posix(),
                    },
                ],
                MODULE.PATH_MAP_FIELDS,
            )
            raw = root / "raw.tsv"
            write_tsv(
                raw,
                [
                    {
                        "Ref_file": (assemblies / "development" / "D1" / "assembly.fna").as_posix(),
                        "Query_file": (assemblies / "external" / "E1" / "assembly.fna").as_posix(),
                        "ANI": "99.99",
                        "Align_fraction_ref": "98.21",
                        "Align_fraction_query": "99.68",
                    }
                ],
                MODULE.SKANI_RAW_FIELDS,
            )
            output = root / "normalized.tsv"
            payload = MODULE.normalize_skani_search(raw, path_map, output, root / "normalize.json")
            rows = MODULE._read_tsv(output)
            self.assertEqual(payload["counts"]["candidate_comparisons"], 1)
            self.assertEqual(rows[0]["development_sequence_id"], "development::D1")
            self.assertEqual(rows[0]["external_sequence_id"], "external::E1")
            self.assertEqual(rows[0]["aligned_fraction_development"], "0.982100")
            self.assertEqual(rows[0]["aligned_fraction_external"], "0.996800")

    def test_normalize_accepts_header_only_zero_hits_but_rejects_missing_schema(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, assemblies, _ = self._fixture(root)
            path_map = root / "paths.tsv"
            write_tsv(
                path_map,
                [
                    {
                        "sequence_id": "development::D1",
                        "cohort": "development",
                        "isolate_id": "D1",
                        "assembly_path": (assemblies / "development" / "D1" / "assembly.fna").as_posix(),
                    },
                    {
                        "sequence_id": "external::E1",
                        "cohort": "external",
                        "isolate_id": "E1",
                        "assembly_path": (assemblies / "external" / "E1" / "assembly.fna").as_posix(),
                    },
                ],
                MODULE.PATH_MAP_FIELDS,
            )
            raw = root / "raw.tsv"
            write_tsv(raw, [], MODULE.SKANI_RAW_FIELDS)
            output = root / "normalized.tsv"
            payload = MODULE.normalize_skani_search(raw, path_map, output, root / "normalize.json")
            self.assertEqual(payload["counts"]["candidate_comparisons"], 0)
            self.assertEqual(MODULE._read_tsv(output), [])

            raw.write_text("Ref_file\tQuery_file\n", encoding="utf-8")
            with self.assertRaisesRegex(MODULE.DisjointnessError, "missing required columns"):
                MODULE.normalize_skani_search(raw, path_map, output, root / "normalize.json")


if __name__ == "__main__":
    unittest.main()


CONTRACT_STUB = {
    "comparison": {
        "near_neighbor_min_ani_percent": 99.9,
        "near_neighbor_min_reciprocal_aligned_fraction": 0.9,
        "putative_duplicate_min_ani_percent": 99.99,
        "putative_duplicate_min_reciprocal_aligned_fraction": 0.95,
    }
}


def _cluster_fixture(root: Path, comparisons: list[dict[str, str]]) -> dict:
    """Run cluster-development over a five-isolate development request table."""
    requests = root / "requests.tsv"
    with requests.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sequence_id", "cohort"], delimiter="\t")
        writer.writeheader()
        for index in range(1, 6):
            writer.writerow({"sequence_id": f"development::D{index}", "cohort": "development"})
    contract = root / "contract.json"
    contract.write_text(json.dumps(CONTRACT_STUB), encoding="utf-8")
    comparisons_path = root / "comparisons.tsv"
    fields = [
        "development_sequence_id",
        "external_sequence_id",
        "ani_percent",
        "aligned_fraction_development",
        "aligned_fraction_external",
    ]
    with comparisons_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(comparisons)
    return MODULE.cluster_development(
        requests, contract, comparisons_path, root / "clusters.tsv", root / "manifest.json"
    )


def _edge(left: str, right: str, ani: float, af: float) -> dict[str, str]:
    return {
        "development_sequence_id": f"development::{left}",
        "external_sequence_id": f"development::{right}",
        "ani_percent": f"{ani}",
        "aligned_fraction_development": f"{af}",
        "aligned_fraction_external": f"{af}",
    }


class DevelopmentClusteringTests(unittest.TestCase):
    def test_every_isolate_gets_a_cluster_including_singletons(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            payload = _cluster_fixture(root, [_edge("D1", "D2", 99.95, 0.95)])
            self.assertEqual(payload["counts"]["isolates"], 5)
            # D1+D2 merge; D3, D4, D5 remain singletons.
            self.assertEqual(payload["counts"]["clusters"], 4)
            self.assertEqual(payload["counts"]["singleton_clusters"], 3)
            self.assertEqual(payload["counts"]["largest_cluster_size"], 2)
            rows = list(csv.DictReader((root / "clusters.tsv").open(encoding="utf-8"), delimiter="\t"))
            self.assertEqual(len(rows), 5)
            self.assertEqual({row["sequence_id"] for row in rows}, {f"development::D{i}" for i in range(1, 6)})

    def test_edges_below_the_frozen_threshold_do_not_merge(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            payload = _cluster_fixture(
                root,
                [
                    _edge("D1", "D2", 99.89, 0.95),  # ANI below 99.9
                    _edge("D3", "D4", 99.95, 0.89),  # aligned fraction below 0.9
                ],
            )
            self.assertEqual(payload["counts"]["linking_edges"], 0)
            self.assertEqual(payload["counts"]["clusters"], 5)

    def test_clusters_are_transitive(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            payload = _cluster_fixture(
                root,
                [_edge("D1", "D2", 99.95, 0.95), _edge("D2", "D3", 99.95, 0.95)],
            )
            self.assertEqual(payload["counts"]["largest_cluster_size"], 3)
            rows = list(csv.DictReader((root / "clusters.tsv").open(encoding="utf-8"), delimiter="\t"))
            assignment = {row["sequence_id"]: row["genomic_cluster"] for row in rows}
            self.assertEqual(assignment["development::D1"], assignment["development::D3"])

    def test_self_comparisons_are_ignored(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            payload = _cluster_fixture(root, [_edge("D1", "D1", 100.0, 1.0)])
            self.assertEqual(payload["counts"]["linking_edges"], 0)
            self.assertEqual(payload["counts"]["clusters"], 5)

    def test_non_development_sequence_is_rejected(self) -> None:
        """Guards against clustering a cross-cohort comparison table by mistake."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaises(MODULE.DisjointnessError):
                _cluster_fixture(
                    root,
                    [
                        {
                            "development_sequence_id": "development::D1",
                            "external_sequence_id": "external::JB-1",
                            "ani_percent": "99.95",
                            "aligned_fraction_development": "0.95",
                            "aligned_fraction_external": "0.95",
                        }
                    ],
                )

    def test_thresholds_are_taken_from_the_frozen_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            payload = _cluster_fixture(root, [_edge("D1", "D2", 99.95, 0.95)])
            self.assertEqual(payload["thresholds"]["near_neighbor_min_ani_percent"], 99.9)
            self.assertEqual(payload["thresholds"]["near_neighbor_min_reciprocal_aligned_fraction"], 0.9)
            self.assertFalse(payload["scientific_boundary"]["phenotypes_read"])


class NormalizeQueryCohortTests(unittest.TestCase):
    """A development-versus-development search must normalize, not fail closed."""

    def _run(self, query_cohort: str, query_is_development: bool) -> dict:
        root = Path(self._tmp.name)
        path_map = root / "path-map.tsv"
        rows = [
            {
                "sequence_id": "development::D1",
                "cohort": "development",
                "isolate_id": "D1",
                "assembly_path": "a/D1.fna",
            },
            {
                "sequence_id": ("development::D2" if query_is_development else "external::E1"),
                "cohort": ("development" if query_is_development else "external"),
                "isolate_id": ("D2" if query_is_development else "E1"),
                "assembly_path": "a/Q.fna",
            },
        ]
        with path_map.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MODULE.PATH_MAP_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        raw = root / "raw.tsv"
        with raw.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MODULE.SKANI_RAW_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerow(
                {
                    "Ref_file": "a/D1.fna",
                    "Query_file": "a/Q.fna",
                    "ANI": "99.95",
                    "Align_fraction_ref": "95.0",
                    "Align_fraction_query": "96.0",
                }
            )
        return MODULE.normalize_skani_search(
            raw,
            path_map,
            root / "out.tsv",
            root / "manifest.json",
            query_cohort=query_cohort,
        )

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)

    def test_development_query_cohort_is_accepted(self) -> None:
        payload = self._run("development", query_is_development=True)
        self.assertEqual(payload["query_cohort"], "development")

    def test_external_query_cohort_remains_the_default_behaviour(self) -> None:
        payload = self._run("external", query_is_development=False)
        self.assertEqual(payload["query_cohort"], "external")

    def test_cohort_mismatch_still_fails_closed(self) -> None:
        with self.assertRaises(MODULE.DisjointnessError):
            self._run("development", query_is_development=False)
        with self.assertRaises(MODULE.DisjointnessError):
            self._run("external", query_is_development=True)

    def test_unsupported_query_cohort_is_rejected(self) -> None:
        with self.assertRaises(MODULE.DisjointnessError):
            self._run("nonsense", query_is_development=True)


def _membership_fixture(root: Path, collisions: list[dict[str, str]], rule: str) -> dict:
    requests = root / "requests.tsv"
    with requests.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sequence_id", "cohort", "isolate_id"], delimiter="\t"
        )
        writer.writeheader()
        for index in range(1, 5):
            writer.writerow(
                {
                    "sequence_id": f"external::E{index}",
                    "cohort": "external",
                    "isolate_id": f"E{index}",
                }
            )
        writer.writerow(
            {"sequence_id": "development::D1", "cohort": "development", "isolate_id": "D1"}
        )
    collisions_path = root / "collisions.tsv"
    with collisions_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MODULE.COLLISION_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(collisions)
    return MODULE.assign_external_membership(
        requests,
        collisions_path,
        root / "membership.tsv",
        root / "manifest.json",
        rule=rule,
    )


def _collision(external: str, collision_class: str) -> dict[str, str]:
    return {
        "component_id": "c1",
        "development_sequence_id": "development::D1",
        "external_sequence_id": f"external::{external}",
        "ani_percent": "99.995" if collision_class == "PUTATIVE_GENOMIC_DUPLICATE" else "99.93",
        "aligned_fraction_development": "0.97",
        "aligned_fraction_external": "0.97",
        "minimum_reciprocal_aligned_fraction": "0.97",
        "collision_class": collision_class,
        "resolution_status": "UNRESOLVED",
    }


class ExternalMembershipTests(unittest.TestCase):
    """Amendment 003: the duplicate rule excludes clones, not lineage members."""

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        # E1 duplicate, E2 near-neighbour only, E3 both, E4 no collision.
        self.collisions = [
            _collision("E1", "PUTATIVE_GENOMIC_DUPLICATE"),
            _collision("E2", "NEAR_NEIGHBOR"),
            _collision("E3", "PUTATIVE_GENOMIC_DUPLICATE"),
            _collision("E3", "NEAR_NEIGHBOR"),
        ]

    def test_duplicate_rule_retains_near_neighbour_only_isolates(self) -> None:
        payload = _membership_fixture(self.root, self.collisions, "duplicate")
        self.assertEqual(payload["counts"]["external_locked"], 2)  # E2 and E4
        self.assertEqual(payload["counts"]["external_excluded"], 2)  # E1 and E3
        rows = list(csv.DictReader((self.root / "membership.tsv").open(encoding="utf-8"), delimiter="\t"))
        by_id = {row["sequence_id"]: row for row in rows}
        self.assertEqual(by_id["external::E2"]["membership"], "external_locked")
        self.assertEqual(by_id["external::E1"]["membership"], "external_excluded")

    def test_near_neighbour_rule_is_stricter(self) -> None:
        payload = _membership_fixture(self.root, self.collisions, "near_neighbor")
        self.assertEqual(payload["counts"]["external_locked"], 1)  # only E4
        self.assertEqual(payload["counts"]["external_excluded"], 3)

    def test_lineage_novel_flag_is_independent_of_the_rule(self) -> None:
        for rule in ("duplicate", "near_neighbor"):
            payload = _membership_fixture(self.root, self.collisions, rule)
            self.assertEqual(payload["counts"]["lineage_novel"], 1, rule)

    def test_partner_counts_are_recorded(self) -> None:
        _membership_fixture(self.root, self.collisions, "duplicate")
        rows = list(csv.DictReader((self.root / "membership.tsv").open(encoding="utf-8"), delimiter="\t"))
        by_id = {row["sequence_id"]: row for row in rows}
        self.assertEqual(by_id["external::E3"]["n_duplicate_partners"], "1")
        self.assertEqual(by_id["external::E3"]["n_near_neighbor_partners"], "1")
        self.assertEqual(by_id["external::E4"]["n_duplicate_partners"], "0")

    def test_unknown_rule_is_rejected(self) -> None:
        with self.assertRaises(MODULE.DisjointnessError):
            _membership_fixture(self.root, self.collisions, "anything")

    def test_unknown_collision_class_fails_closed(self) -> None:
        bad = [_collision("E1", "SOMETHING_NEW")]
        with self.assertRaises(MODULE.DisjointnessError):
            _membership_fixture(self.root, bad, "duplicate")

    def test_no_phenotype_is_read(self) -> None:
        payload = _membership_fixture(self.root, self.collisions, "duplicate")
        self.assertFalse(payload["scientific_boundary"]["phenotypes_read"])
        self.assertEqual(payload["exclusion_rule"], "duplicate")


class ExternalClusteringTests(unittest.TestCase):
    """The 1,888 locked externals need clusters of their own (blocker 2)."""

    def _run(self, cohort: str) -> dict:
        root = Path(self._tmp.name)
        requests = root / "requests.tsv"
        with requests.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sequence_id", "cohort"], delimiter="\t")
            writer.writeheader()
            for index in range(1, 4):
                writer.writerow({"sequence_id": f"external::E{index}", "cohort": "external"})
            writer.writerow({"sequence_id": "development::D1", "cohort": "development"})
        contract = root / "contract.json"
        contract.write_text(json.dumps(CONTRACT_STUB), encoding="utf-8")
        comparisons = root / "comparisons.tsv"
        fields = [
            "development_sequence_id",
            "external_sequence_id",
            "ani_percent",
            "aligned_fraction_development",
            "aligned_fraction_external",
        ]
        with comparisons.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerow(
                {
                    "development_sequence_id": "external::E1",
                    "external_sequence_id": "external::E2",
                    "ani_percent": "99.95",
                    "aligned_fraction_development": "0.95",
                    "aligned_fraction_external": "0.95",
                }
            )
        return MODULE.cluster_development(
            requests, contract, comparisons, root / "c.tsv", root / "m.json", cohort=cohort
        )

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)

    def test_external_cohort_clusters_and_is_namespaced(self) -> None:
        payload = self._run("external")
        self.assertEqual(payload["cohort"], "external")
        self.assertEqual(payload["counts"]["isolates"], 3)
        self.assertEqual(payload["counts"]["clusters"], 2)
        rows = list(csv.DictReader((Path(self._tmp.name) / "c.tsv").open(encoding="utf-8"), delimiter="\t"))
        # Namespacing keeps external clusters from colliding with development ones.
        self.assertTrue(all(r["genomic_cluster"].startswith("external-cluster-") for r in rows))

    def test_development_cohort_rejects_external_sequences(self) -> None:
        with self.assertRaises(MODULE.DisjointnessError):
            self._run("development")

    def test_unknown_cohort_is_rejected(self) -> None:
        with self.assertRaises(MODULE.DisjointnessError):
            self._run("nonsense")


def _reference_normalize(raw_rows, path_map, reference_cohort, query_cohort):
    """Frozen copy of the pre-streaming in-memory normalizer, used as an oracle."""
    normalized = []
    seen = set()
    for row in raw_rows:
        reference = path_map[Path(row["Ref_file"]).as_posix()]
        query = path_map[Path(row["Query_file"]).as_posix()]
        assert reference["cohort"] == reference_cohort and query["cohort"] == query_cohort
        pair = (reference["sequence_id"], query["sequence_id"])
        assert pair not in seen
        seen.add(pair)
        ani = float(row["ANI"])
        af_reference = float(row["Align_fraction_ref"])
        af_query = float(row["Align_fraction_query"])
        normalized.append(
            {
                "development_sequence_id": reference["sequence_id"],
                "external_sequence_id": query["sequence_id"],
                "ani_percent": f"{ani:.6f}",
                "aligned_fraction_development": f"{af_reference / 100:.6f}",
                "aligned_fraction_external": f"{af_query / 100:.6f}",
            }
        )
    normalized.sort(key=lambda row: (row["development_sequence_id"], row["external_sequence_id"]))
    lines = ["\t".join(MODULE.NORMALIZED_COMPARISON_FIELDS)]
    for row in normalized:
        lines.append("\t".join(row[field] for field in MODULE.NORMALIZED_COMPARISON_FIELDS))
    return "\n".join(lines) + "\n"


class StreamingNormalizerEquivalenceTests(unittest.TestCase):
    """The streaming normalizer must be byte-identical to the in-memory one."""

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def _fixture(self, *, reference_cohort, query_cohort, n=37, seed=7, extra_rows=()):
        import random

        rng = random.Random(seed)
        # Mixed-length IDs so string order differs from insertion and numeric order.
        ids = [f"{reference_cohort}::ID{rng.randint(1, 10**rng.randint(1, 6))}_{i}" for i in range(n)]
        query_ids = ids if reference_cohort == query_cohort else [
            f"{query_cohort}::Q{rng.randint(1, 10**rng.randint(1, 6))}_{i}" for i in range(n)
        ]
        map_rows = []
        path_map = {}
        for cohort, pool in ((reference_cohort, ids), (query_cohort, query_ids)):
            for sequence_id in pool:
                path = f"a/{cohort}/{sequence_id.split('::')[1]}.fna"
                if path in path_map:
                    continue
                row = {"sequence_id": sequence_id, "cohort": cohort, "isolate_id": sequence_id, "assembly_path": path}
                map_rows.append(row)
                path_map[path] = row
        ref_paths = [p for p, r in path_map.items() if r["cohort"] == reference_cohort]
        query_paths = [p for p, r in path_map.items() if r["cohort"] == query_cohort]
        pairs = [(r, q) for r in ref_paths for q in query_paths]
        rng.shuffle(pairs)
        pairs = pairs[: max(1, len(pairs) // 2)]
        # Awkward float spellings exercise the float -> .6f round trip.
        ani_values = ["100.00", "99.9", "99.999999", "95.123456789", "97", "99.905"]
        af_values = ["99.89", "100", "50.5", "90.0000001", "0", "95.123456"]
        raw_rows = [
            {
                "Ref_file": r,
                "Query_file": q,
                "ANI": rng.choice(ani_values),
                "Align_fraction_ref": rng.choice(af_values),
                "Align_fraction_query": rng.choice(af_values),
                "Ref_name": "x",
                "Query_name": "y",
            }
            for r, q in pairs
        ]
        raw_rows.extend(extra_rows)
        with (self.root / "map.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MODULE.PATH_MAP_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(map_rows)
        fields = MODULE.SKANI_RAW_FIELDS + ["Ref_name", "Query_name"]
        with (self.root / "raw.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(raw_rows)
        return raw_rows, path_map

    def _run(self, reference_cohort, query_cohort, chunk_rows):
        return MODULE.normalize_skani_search(
            self.root / "raw.tsv",
            self.root / "map.tsv",
            self.root / "out.tsv",
            self.root / "manifest.json",
            reference_cohort=reference_cohort,
            query_cohort=query_cohort,
            chunk_rows=chunk_rows,
        )

    def test_byte_identical_to_reference_across_modes_and_chunk_sizes(self) -> None:
        for reference_cohort, query_cohort in (
            ("development", "external"),
            ("development", "development"),
            ("external", "external"),
        ):
            for chunk_rows in (1, 3, 10_000):
                with self.subTest(ref=reference_cohort, query=query_cohort, chunk=chunk_rows):
                    raw_rows, path_map = self._fixture(
                        reference_cohort=reference_cohort, query_cohort=query_cohort
                    )
                    payload = self._run(reference_cohort, query_cohort, chunk_rows)
                    expected = _reference_normalize(raw_rows, path_map, reference_cohort, query_cohort)
                    observed = (self.root / "out.tsv").read_text(encoding="utf-8")
                    self.assertEqual(observed, expected)
                    self.assertEqual(payload["counts"]["candidate_comparisons"], len(raw_rows))

    def test_duplicate_pair_is_still_rejected(self) -> None:
        raw_rows, _ = self._fixture(reference_cohort="development", query_cohort="external")
        self._fixture(
            reference_cohort="development", query_cohort="external", extra_rows=[dict(raw_rows[0])]
        )
        with self.assertRaisesRegex(MODULE.DisjointnessError, "Duplicate skani comparison pair"):
            self._run("development", "external", 4)

    def test_empty_search_writes_header_only(self) -> None:
        self._fixture(reference_cohort="development", query_cohort="external")
        fields = MODULE.SKANI_RAW_FIELDS
        (self.root / "raw.tsv").write_text("\t".join(fields) + "\n", encoding="utf-8")
        payload = self._run("development", "external", 5)
        self.assertEqual(payload["counts"]["candidate_comparisons"], 0)
        self.assertEqual(
            (self.root / "out.tsv").read_text(encoding="utf-8"),
            "\t".join(MODULE.NORMALIZED_COMPARISON_FIELDS) + "\n",
        )

    def test_truncated_row_fails_closed(self) -> None:
        self._fixture(reference_cohort="development", query_cohort="external")
        with (self.root / "raw.tsv").open("a", encoding="utf-8") as handle:
            handle.write("a/only-one-column\n")
        with self.assertRaisesRegex(MODULE.DisjointnessError, "truncated"):
            self._run("development", "external", 5)

    def test_out_of_range_value_fails_closed(self) -> None:
        raw_rows, _ = self._fixture(reference_cohort="development", query_cohort="external")
        bad = dict(raw_rows[0])
        bad["ANI"] = "100.5"
        bad["Query_file"] = raw_rows[0]["Query_file"]
        self._fixture(
            reference_cohort="development", query_cohort="external", extra_rows=[bad]
        )
        with self.assertRaises(MODULE.DisjointnessError):
            self._run("development", "external", 5)

    def test_distinct_query_and_reference_counts_are_recorded(self) -> None:
        self._fixture(reference_cohort="development", query_cohort="development", n=5)
        payload = self._run("development", "development", 2)
        self.assertGreater(payload["counts"]["distinct_queries"], 0)
        self.assertGreater(payload["counts"]["distinct_references"], 0)
        self.assertEqual(payload["implementation"], "streaming")


class IterativeClusteringTests(unittest.TestCase):
    def test_long_chain_does_not_hit_the_recursion_limit(self) -> None:
        """A recursive find overflows on single-linkage chains this long."""
        import sys

        n = sys.getrecursionlimit() * 3
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (root / "requests.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["sequence_id", "cohort"], delimiter="\t")
                writer.writeheader()
                for index in range(n):
                    writer.writerow({"sequence_id": f"development::D{index:06d}", "cohort": "development"})
            (root / "contract.json").write_text(json.dumps(CONTRACT_STUB), encoding="utf-8")
            fields = [
                "development_sequence_id",
                "external_sequence_id",
                "ani_percent",
                "aligned_fraction_development",
                "aligned_fraction_external",
            ]
            with (root / "comparisons.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                # Link in descending order, which builds the deepest parent chain.
                for index in range(n - 1, 0, -1):
                    writer.writerow(
                        {
                            "development_sequence_id": f"development::D{index:06d}",
                            "external_sequence_id": f"development::D{index - 1:06d}",
                            "ani_percent": "99.95",
                            "aligned_fraction_development": "0.95",
                            "aligned_fraction_external": "0.95",
                        }
                    )
            payload = MODULE.cluster_development(
                root / "requests.tsv",
                root / "contract.json",
                root / "comparisons.tsv",
                root / "clusters.tsv",
                root / "manifest.json",
            )
            self.assertEqual(payload["counts"]["clusters"], 1)
            self.assertEqual(payload["counts"]["largest_cluster_size"], n)
            self.assertEqual(payload["counts"]["candidate_comparisons"], n - 1)


class DuplicatePairExtractionTests(unittest.TestCase):
    def test_keeps_only_frozen_duplicate_pairs_and_drops_self(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "contract.json").write_text(json.dumps(CONTRACT_STUB), encoding="utf-8")
            rows = [
                ("development::A", "development::A", "100.000000", "1.000000", "1.000000"),  # self
                ("development::A", "development::B", "99.995000", "0.960000", "0.970000"),  # keep
                ("development::A", "development::C", "99.995000", "0.940000", "0.990000"),  # AF too low
                ("development::B", "development::C", "99.980000", "0.990000", "0.990000"),  # ANI too low
                ("development::B", "development::A", "99.990000", "0.950000", "0.950000"),  # boundary keep
            ]
            with (root / "cmp.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(MODULE.NORMALIZED_COMPARISON_FIELDS)
                writer.writerows(rows)
            payload = MODULE.extract_duplicate_pairs(
                root / "contract.json", root / "cmp.tsv", root / "pairs.tsv", root / "m.json"
            )
            self.assertEqual(payload["counts"]["comparisons_scanned"], 5)
            self.assertEqual(payload["counts"]["duplicate_pair_rows"], 2)
            kept = list(csv.DictReader((root / "pairs.tsv").open(encoding="utf-8"), delimiter="\t"))
            self.assertEqual(
                {(r["development_sequence_id"], r["external_sequence_id"]) for r in kept},
                {("development::A", "development::B"), ("development::B", "development::A")},
            )
            self.assertEqual(payload["thresholds"]["putative_duplicate_min_ani_percent"], 99.99)
