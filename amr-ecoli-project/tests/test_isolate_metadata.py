from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

MODULE_PATH = Path(__file__).parents[1] / "scripts" / "build_isolate_metadata.py"
SPEC = importlib.util.spec_from_file_location("build_isolate_metadata", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

PAIR_FIELDS = ["development_sequence_id", "external_sequence_id", "ani_percent",
               "aligned_fraction_development", "aligned_fraction_external"]


def write(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class Fixture:
    """Development isolates D1..Dn with chosen lineages, plus external E1."""

    def __init__(self, root: Path, lineages: dict[str, str], *, pairs=(), shared_biosample=None,
                 genome_overrides=None, membership_state="external_locked", with_ext_cluster=True,
                 specimen_missing=()):
        self.root = root
        isolates = list(lineages) + ["E1"]
        requests, clusters, genomes, specimens = [], [], [], []
        for isolate in isolates:
            cohort = "external" if isolate.startswith("E") else "development"
            biosample = f"SAMN{isolate}"
            if shared_biosample and isolate in shared_biosample:
                biosample = "SAMN_SHARED"
            requests.append({"sequence_id": f"{cohort}::{isolate}", "cohort": cohort, "isolate_id": isolate,
                             "biosample_accession": biosample, "sequence_accession": f"GCA_{isolate}.1"})
            if cohort == "development" or with_ext_cluster:
                clusters.append({"sequence_id": f"{cohort}::{isolate}", "genomic_cluster": f"{cohort}-cluster-{isolate}"})
            genome = {"isolate_id": isolate, "assembly_accession": f"GCA_{isolate}.1", "species_method": "ANI",
                      "species_result": "Escherichia coli", "genome_qc_status": "PASS",
                      "lineage_group": lineages.get(isolate, "ST999"), "eligible_for_modeling": "true",
                      "genome_source": "registered_assembly", "genome_source_accession": f"GCA_{isolate}.1",
                      "assembly_sha256": "e" * 64, "exclusion_reasons": ""}
            if genome_overrides:
                genome.update(genome_overrides)
            genomes.append(genome)
            if isolate not in specimen_missing:
                specimens.append({"isolate_id": isolate, "specimen_source": "urine",
                                  "intended_use_population": "human_clinical"})
        write(root / "requests.tsv", list(requests[0]), requests)
        write(root / "clusters.tsv", ["sequence_id", "genomic_cluster"], clusters)
        write(root / "membership.tsv", ["sequence_id", "membership"],
              [{"sequence_id": "external::E1", "membership": membership_state}])
        write(root / "genomes.tsv", list(genomes[0]), genomes)
        write(root / "specimens.tsv", ["isolate_id", "specimen_source", "intended_use_population"], specimens)
        write(root / "pairs.tsv", PAIR_FIELDS, [
            {"development_sequence_id": f"development::{a}", "external_sequence_id": f"development::{b}",
             "ani_percent": str(ani), "aligned_fraction_development": "0.99", "aligned_fraction_external": "0.99"}
            for a, b, ani in pairs
        ])

    def build(self, **kwargs):
        r = self.root
        return MODULE.build_isolate_metadata(
            r / "requests.tsv", r / "clusters.tsv", r / "membership.tsv", r / "genomes.tsv",
            r / "specimens.tsv", r / "pairs.tsv", r / "metadata.csv", r / "exclusions.csv",
            r / "ledger.csv", r / "manifest.json", **kwargs)

    def rows(self, name="metadata.csv"):
        return {row.get("isolate_id") or row.get("isolate_a"): row
                for row in csv.DictReader((self.root / name).open(encoding="utf-8"))}


def lineages(n: int, *, st: str = "ST1") -> dict[str, str]:
    return {f"D{i:02d}": st for i in range(1, n + 1)}


class BaseTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def test_happy_path_assigns_split_cluster_and_singleton_dedup(self) -> None:
        f = Fixture(self.root, {"D01": "ST1", "D02": "ST2"})
        payload = f.build()
        self.assertEqual(payload["counts"]["development"], 2)
        self.assertEqual(payload["counts"]["external"], 1)
        rows = f.rows()
        self.assertEqual(rows["D01"]["deduplication_group"], "dedup-D01")
        self.assertEqual(rows["E1"]["evaluation_split"], "external")

    def test_missing_attribute_row_excludes_because_population_is_absent(self) -> None:
        f = Fixture(self.root, {"D01": "ST1", "D02": "ST1"}, specimen_missing=("D02",))
        payload = f.build()
        self.assertEqual(payload["exclusion_counts"].get("required_field_absent"), 1)
        self.assertEqual(f.rows()["D01"]["specimen_source"], "urine")

    def test_blank_specimen_source_is_kept_blank_not_excluded(self) -> None:
        f = Fixture(self.root, {"D01": "ST1", "D02": "ST1"})
        rows = list(csv.DictReader((self.root / "specimens.tsv").open(encoding="utf-8"), delimiter="\t"))
        rows[1]["specimen_source"] = ""
        write(self.root / "specimens.tsv", ["isolate_id", "specimen_source", "intended_use_population"], rows)
        payload = f.build()
        self.assertEqual(payload["counts"]["development"], 2)
        self.assertEqual(payload["counts"]["specimen_source_absent"], 1)
        self.assertEqual(f.rows()["D02"]["specimen_source"], "")

    def test_absent_lineage_is_not_substituted(self) -> None:
        f = Fixture(self.root, {"D01": "ST1"}, genome_overrides={"lineage_group": ""})
        payload = f.build()
        self.assertEqual(payload["counts"]["metadata_rows"], 0)

    def test_membership_excluded_external_is_dropped(self) -> None:
        f = Fixture(self.root, {"D01": "ST1"}, membership_state="external_excluded")
        payload = f.build()
        self.assertIn("external_excluded_by_membership", payload["exclusion_counts"])

    def test_missing_external_cluster_is_reported_not_invented(self) -> None:
        f = Fixture(self.root, {"D01": "ST1"}, with_ext_cluster=False)
        payload = f.build()
        self.assertIn("genomic_cluster_absent", payload["exclusion_counts"])


class RuleITests(BaseTests):
    def test_shared_biosample_merges_dedup_groups(self) -> None:
        f = Fixture(self.root, {"D01": "ST1", "D02": "ST2"}, shared_biosample={"D01", "D02"})
        payload = f.build()
        self.assertEqual(payload["rule_i"]["identifier_merges"], 1)
        rows = f.rows()
        self.assertEqual(rows["D01"]["deduplication_group"], rows["D02"]["deduplication_group"])

    def test_identifier_shared_across_the_split_is_refused(self) -> None:
        f = Fixture(self.root, {"D01": "ST1"}, shared_biosample={"D01", "E1"})
        with self.assertRaisesRegex(MODULE.MetadataError, "Deduplication groups cross"):
            f.build()


class RuleGTests(BaseTests):
    def test_within_lineage_duplicates_need_no_action(self) -> None:
        f = Fixture(self.root, lineages(4), pairs=[("D01", "D02", 99.995)])
        payload = f.build()
        self.assertEqual(payload["rule_g"]["within_lineage"], 1)
        self.assertEqual(payload["rule_g"]["cross_lineage"], 0)
        self.assertEqual(f.rows()["D01"]["deduplication_group"], "dedup-D01")

    def test_symmetric_pair_rows_count_once(self) -> None:
        f = Fixture(self.root, lineages(4), pairs=[("D01", "D02", 99.995), ("D02", "D01", 99.996)])
        self.assertEqual(f.build()["rule_g"]["genomic_duplicate_pairs"], 1)

    def test_cross_lineage_pair_merges_under_the_cap(self) -> None:
        lin = {f"D{i:02d}": f"ST{i}" for i in range(1, 41)}  # 40 singleton lineages, cap = 4
        many_within = [(f"D{i:02d}", f"D{i:02d}", 100.0) for i in range(1, 2)]  # self rows ignored
        f = Fixture(self.root, lin, pairs=[("D01", "D02", 99.995)] + many_within)
        payload = f.build(reopen_cross_lineage_fraction=1.0)
        self.assertEqual(payload["rule_g"]["merged"], 1)
        rows = f.rows()
        self.assertEqual(rows["D01"]["deduplication_group"], rows["D02"]["deduplication_group"])

    def test_cross_lineage_pair_over_the_cap_excludes_the_smaller_lineage(self) -> None:
        lin = {**lineages(30, st="ST_BIG"), "D90": "ST_SMALL"}  # cap = 3; merged would be 31
        f = Fixture(self.root, lin, pairs=[("D01", "D90", 99.995)])
        payload = f.build(reopen_cross_lineage_fraction=1.0)
        self.assertEqual(payload["rule_g"]["excluded"], 1)
        self.assertNotIn("D90", f.rows())
        self.assertIn("D01", f.rows())
        ledger = list(csv.DictReader((self.root / "ledger.csv").open(encoding="utf-8")))
        self.assertEqual(ledger[0]["action"], "excluded:D90")

    def test_equal_lineage_sizes_exclude_the_larger_identifier(self) -> None:
        lin = {**lineages(10, st="STA"), **{f"D{i:02d}": "STB" for i in range(11, 21)}}  # cap 2
        f = Fixture(self.root, lin, pairs=[("D03", "D15", 99.995)])
        f.build(reopen_cross_lineage_fraction=1.0)
        self.assertNotIn("D15", f.rows())
        self.assertIn("D03", f.rows())

    def test_too_many_cross_lineage_pairs_reopens_amendment_004(self) -> None:
        lin = {f"D{i:02d}": f"ST{i}" for i in range(1, 11)}
        f = Fixture(self.root, lin, pairs=[("D01", "D02", 99.995)])
        with self.assertRaisesRegex(MODULE.MetadataError, "amendment 004 must be reopened"):
            f.build()

    def test_merges_cannot_chain_lineages_past_the_cap(self) -> None:
        """The anti-chaining guarantee: a run of cross-lineage pairs stops at the cap."""
        lin = {f"D{i:02d}": f"ST{i}" for i in range(1, 51)}  # 50 singletons, cap = 5
        chain = [(f"D{i:02d}", f"D{i + 1:02d}", 99.999 - i * 1e-4) for i in range(1, 20)]
        f = Fixture(self.root, lin, pairs=chain)
        payload = f.build(reopen_cross_lineage_fraction=1.0)
        self.assertLessEqual(payload["rule_g"]["largest_cv_group"], 5)
        self.assertGreater(payload["rule_g"]["excluded"], 0)

    def test_no_phenotype_is_read(self) -> None:
        f = Fixture(self.root, {"D01": "ST1"})
        self.assertFalse(f.build()["scientific_boundary"]["phenotypes_read"])


class GenomeEligibilityTests(BaseTests):
    def test_ineligible_genome_is_excluded_with_its_reasons(self) -> None:
        f = Fixture(self.root, {"D01": "ST1", "D02": "ST1"},
                    genome_overrides=None)
        # Mark only D02 ineligible by rewriting its genome row.
        rows = list(csv.DictReader((self.root / "genomes.tsv").open(encoding="utf-8"), delimiter="\t"))
        for row in rows:
            if row["isolate_id"] == "D02":
                row["eligible_for_modeling"] = "false"
                row["exclusion_reasons"] = "SPECIES_NOT_CONFIRMED"
        write(self.root / "genomes.tsv", list(rows[0]), rows)
        payload = f.build()
        self.assertEqual(payload["exclusion_counts"].get("genome_not_eligible"), 1)
        self.assertNotIn("D02", f.rows())
        self.assertEqual(f.rows()["D01"]["genome_sha256"], "e" * 64)


class SplitClusterFileTests(BaseTests):
    """Development and external clusterings arrive as separate files."""

    def _split_clusters(self) -> tuple[Path, Path]:
        rows = list(csv.DictReader((self.root / "clusters.tsv").open(encoding="utf-8"), delimiter="\t"))
        dev = [r for r in rows if r["sequence_id"].startswith("development::")]
        ext = [r for r in rows if r["sequence_id"].startswith("external::")]
        write(self.root / "dev-clusters.tsv", ["sequence_id", "genomic_cluster"], dev)
        write(self.root / "ext-clusters.tsv", ["sequence_id", "genomic_cluster"], ext)
        return self.root / "dev-clusters.tsv", self.root / "ext-clusters.tsv"

    def _build(self, clusters):
        r = self.root
        return MODULE.build_isolate_metadata(
            r / "requests.tsv", clusters, r / "membership.tsv", r / "genomes.tsv",
            r / "specimens.tsv", r / "pairs.tsv", r / "metadata.csv", r / "exclusions.csv",
            r / "ledger.csv", r / "manifest.json")

    def test_two_cluster_files_combine_and_are_both_hashed(self) -> None:
        Fixture(self.root, {"D01": "ST1", "D02": "ST2"})
        payload = self._build(list(self._split_clusters()))
        self.assertEqual(payload["counts"]["development"], 2)
        self.assertEqual(payload["counts"]["external"], 1)
        self.assertIn("clusters_0", payload["inputs"])
        self.assertIn("clusters_1", payload["inputs"])

    def test_development_file_alone_leaves_external_without_cluster(self) -> None:
        Fixture(self.root, {"D01": "ST1", "D02": "ST2"})
        dev, _ = self._split_clusters()
        payload = self._build([dev])
        self.assertEqual(payload["counts"]["external"], 0)

    def test_sequence_in_both_files_is_refused(self) -> None:
        Fixture(self.root, {"D01": "ST1", "D02": "ST2"})
        dev, _ = self._split_clusters()
        with self.assertRaises(MODULE.MetadataError):
            self._build([dev, dev])
