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
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location("finalize_genome_analyses", ROOT / "scripts" / "finalize_genome_analyses.py")
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

ACHTMAN_LOCI = ["adk", "fumC", "gyrB", "icd", "mdh", "purA", "recA"]
PASTEUR_LOCI = ["dinB", "icdA", "pabB", "polB", "putP", "trpA", "trpB", "uidA"]
CONFIG = {"genomics": {"assembly_qc": {"min_total_length": 4000000, "max_total_length": 6500000,
                                       "max_contigs": 500, "min_n50": 20000, "max_ns_per_100kb": 1000.0},
                       "mlst": {"scheme": "ecoli_achtman_4", "loci": ACHTMAN_LOCI}}}


def mlst_text(st: str = "131", status: str = "PERFECT", scheme: str = "ecoli_achtman_4",
              loci: list[str] | None = None) -> str:
    alleles = ";".join(f"{locus}(1)" for locus in (loci or ACHTMAN_LOCI))
    return f"FILE\tSCHEME\tST\tSTATUS\tSCORE\tALLELES\nassembly.fna\t{scheme}\t{st}\t{status}\t100\t{alleles}\n"


def make_genome(root: Path, cohort: str, isolate: str, *, length=5_000_000, contigs=150, n50=120_000,
                ani="99.1", af="88.0", skani_rows=True, st="131", status="PERFECT", amr_exit=0,
                complete=True, scheme="ecoli_achtman_4", loci=None) -> None:
    d = root / cohort / isolate
    for sub in ("quast", "mlst", "amrfinder", "species"):
        (d / sub).mkdir(parents=True, exist_ok=True)
    (d / "quast" / "report.tsv").write_text(
        f"Assembly\tassembly\n# contigs\t{contigs}\nTotal length\t{length}\nN50\t{n50}\nGC (%)\t50.7\n# N's per 100 kbp\t0.00\n",
        encoding="utf-8")
    (d / "mlst" / "mlst.tsv").write_text(mlst_text(st, status, scheme, loci), encoding="utf-8")
    (d / "amrfinder" / "amrfinder.tsv").write_text("Name\tGene symbol\nx\tblaCTX-M-15\nx\tgyrA_S83L\n", encoding="utf-8")
    header = "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
    body = f"ref.fna\tassembly.fna\t{ani}\t85.0\t{af}\tref\tq\n" if skani_rows else ""
    (d / "species" / "skani.tsv").write_text(header + body, encoding="utf-8")
    for sub, code in (("quast", 0), ("mlst", 0), ("amrfinder", amr_exit), ("species", 0)):
        (d / sub / "exit_code.txt").write_text(f"{code}\n", encoding="utf-8")
    (d / "provenance.txt").write_text(f"assembly_sha256={'a' * 64}\n", encoding="utf-8")
    lines = []
    for f in sorted(p for p in d.rglob("*") if p.is_file()):
        lines.append(f"{hashlib.sha256(f.read_bytes()).hexdigest()}  ./{f.relative_to(d).as_posix()}")
    (d / "artifact.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if complete:
        (d / "COMPLETE").touch()


class FinalizeTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.analysis = self.root / "analysis"
        (self.root / "config.json").write_text(json.dumps(CONFIG), encoding="utf-8")

    def _requests(self, rows):
        path = self.root / "requests.tsv"
        with path.open("w", encoding="utf-8", newline="") as h:
            w = csv.DictWriter(h, fieldnames=["sequence_id", "cohort", "isolate_id", "biosample_accession",
                                              "source_type", "sequence_accession"], delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        return path

    def _run(self, rows):
        payload = MODULE.finalize(self._requests(rows), self.analysis, self.root / "config.json",
                                  self.root / "genomes.tsv", self.root / "m.json")
        with (self.root / "genomes.tsv").open(encoding="utf-8") as handle:
            out = {r["isolate_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
        return payload, out

    @staticmethod
    def req(isolate, cohort="development", source="registered_assembly", acc="GCA_000000001.1"):
        return {"sequence_id": f"{cohort}::{isolate}", "cohort": cohort, "isolate_id": isolate,
                "biosample_accession": "SAMN1", "source_type": source, "sequence_accession": acc}

    def _run_with_mlst_root(self, rows, mlst_root):
        return MODULE.finalize(self._requests(rows), self.analysis, self.root / "config.json",
                               self.root / "genomes.tsv", self.root / "m.json", mlst_root=mlst_root)

    def test_pasteur_scheme_stops_the_run(self) -> None:
        # Regression: array 64927454 ran "--scheme ecoli", which is Pasteur in mlst >= 2.23.
        make_genome(self.analysis, "development", "D1", st="43", scheme="ecoli", loci=PASTEUR_LOCI)
        with self.assertRaisesRegex(MODULE.FinalizeError, "ecoli_achtman_4"):
            self._run([self.req("D1")])

    def test_right_scheme_name_with_wrong_loci_stops_the_run(self) -> None:
        make_genome(self.analysis, "development", "D1", loci=ACHTMAN_LOCI[:6])
        with self.assertRaises(MODULE.FinalizeError):
            self._run([self.req("D1")])

    def test_mlst_root_overrides_the_analysis_directory(self) -> None:
        make_genome(self.analysis, "development", "D1", st="43", scheme="ecoli", loci=PASTEUR_LOCI)
        typed = self.root / "mlst" / "development" / "D1"
        typed.mkdir(parents=True)
        (typed / "mlst.tsv").write_text(mlst_text("131"), encoding="utf-8")
        (typed / "exit_code.txt").write_text("0\n", encoding="utf-8")
        (typed / "artifact.sha256").write_text("".join(
            f"{hashlib.sha256((typed / n).read_bytes()).hexdigest()}  {n}\n" for n in ("exit_code.txt", "mlst.tsv")),
            encoding="utf-8")
        (typed / "COMPLETE").touch()
        payload = self._run_with_mlst_root([self.req("D1")], self.root / "mlst")
        with (self.root / "genomes.tsv").open(encoding="utf-8") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(row["lineage_group"], "SLV:ST131")
        self.assertEqual(payload["rules"]["mlst_scheme"], "ecoli_achtman_4")
        # A tampered re-typing output fails closed.
        (typed / "mlst.tsv").write_text(mlst_text("10"), encoding="utf-8")
        with self.assertRaisesRegex(MODULE.FinalizeError, "checksum"):
            self._run_with_mlst_root([self.req("D1")], self.root / "mlst")

    def test_missing_mlst_config_stops_the_run(self) -> None:
        make_genome(self.analysis, "development", "D1")
        (self.root / "config.json").write_text(json.dumps(
            {"genomics": {"assembly_qc": CONFIG["genomics"]["assembly_qc"]}}), encoding="utf-8")
        with self.assertRaisesRegex(MODULE.FinalizeError, "genomics.mlst"):
            self._run([self.req("D1")])

    def test_passing_registered_genome(self) -> None:
        make_genome(self.analysis, "development", "D1")
        _, out = self._run([self.req("D1")])
        row = out["D1"]
        self.assertEqual(row["eligible_for_modeling"], "true")
        self.assertEqual(row["species_result"], "Escherichia coli")
        self.assertEqual(row["genome_qc_status"], "PASS")
        self.assertEqual(row["lineage_group"], "SLV:ST131")
        self.assertEqual(row["assembly_accession"], "GCA_000000001.1")
        self.assertEqual(row["amrfinder_hit_count"], "2")

    def test_raw_read_genome_has_no_invented_assembly_accession(self) -> None:
        make_genome(self.analysis, "external", "JB1")
        _, out = self._run([self.req("JB1", cohort="external", source="raw_reads", acc="DRR000001")])
        self.assertEqual(out["JB1"]["assembly_accession"], "")
        self.assertEqual(out["JB1"]["genome_source_accession"], "DRR000001")
        self.assertEqual(out["JB1"]["genome_source"], "raw_reads")

    def test_species_boundaries(self) -> None:
        make_genome(self.analysis, "development", "LOW_ANI", ani="94.9")
        make_genome(self.analysis, "development", "LOW_AF", af="49.9")
        make_genome(self.analysis, "development", "NO_HIT", skani_rows=False)
        make_genome(self.analysis, "development", "EDGE", ani="95.0", af="50.0")
        _, out = self._run([self.req(i) for i in ("LOW_ANI", "LOW_AF", "NO_HIT", "EDGE")])
        for isolate in ("LOW_ANI", "LOW_AF", "NO_HIT"):
            self.assertEqual(out[isolate]["species_result"], "not_confirmed", isolate)
            self.assertIn("SPECIES_NOT_CONFIRMED", out[isolate]["exclusion_reasons"])
        self.assertEqual(out["EDGE"]["species_result"], "Escherichia coli")

    def test_assembly_qc_thresholds(self) -> None:
        make_genome(self.analysis, "development", "SHORT", length=3_900_000)
        make_genome(self.analysis, "development", "FRAGMENTED", contigs=501)
        make_genome(self.analysis, "development", "LOWN50", n50=19_999)
        _, out = self._run([self.req(i) for i in ("SHORT", "FRAGMENTED", "LOWN50")])
        self.assertIn("ASSEMBLY_LENGTH_BELOW_MIN", out["SHORT"]["exclusion_reasons"])
        self.assertIn("ASSEMBLY_CONTIG_COUNT_ABOVE_MAX", out["FRAGMENTED"]["exclusion_reasons"])
        self.assertIn("ASSEMBLY_N50_BELOW_MIN", out["LOWN50"]["exclusion_reasons"])
        self.assertTrue(all(out[i]["genome_qc_status"] == "FAIL" for i in out))

    def test_unresolved_mlst_leaves_lineage_blank_as_a_warning(self) -> None:
        make_genome(self.analysis, "development", "NOVEL", st="-", status="NOVEL")
        make_genome(self.analysis, "development", "OKAY", st="73", status="OKAY")
        _, out = self._run([self.req("NOVEL"), self.req("OKAY")])
        self.assertEqual(out["NOVEL"]["lineage_group"], "")
        self.assertIn("LINEAGE_UNRESOLVED:NOVEL", out["NOVEL"]["warning_reasons"])
        self.assertEqual(out["NOVEL"]["eligible_for_modeling"], "true")  # handled downstream
        self.assertEqual(out["OKAY"]["lineage_group"], "SLV:ST73")

    def test_incomplete_or_tampered_analysis_is_excluded(self) -> None:
        make_genome(self.analysis, "development", "MISSING", complete=False)
        make_genome(self.analysis, "development", "TAMPERED")
        (self.analysis / "development" / "TAMPERED" / "mlst" / "mlst.tsv").write_text("changed\n", encoding="utf-8")
        _, out = self._run([self.req("MISSING"), self.req("TAMPERED")])
        self.assertIn("ANALYSIS_INCOMPLETE", out["MISSING"]["exclusion_reasons"])
        self.assertIn("ANALYSIS_CHECKSUM_MISMATCH", out["TAMPERED"]["exclusion_reasons"])

    def test_amrfinder_failure_excludes(self) -> None:
        make_genome(self.analysis, "development", "AMRFAIL", amr_exit=1)
        _, out = self._run([self.req("AMRFAIL")])
        self.assertIn("AMRFINDER_EXECUTION_FAILED", out["AMRFAIL"]["exclusion_reasons"])

    def test_manifest_counts_by_cohort(self) -> None:
        make_genome(self.analysis, "development", "D1")
        make_genome(self.analysis, "external", "E1", ani="80")
        payload, _ = self._run([self.req("D1"), self.req("E1", cohort="external", source="raw_reads", acc="DRR1")])
        self.assertEqual(payload["counts"]["development"]["eligible"], 1)
        self.assertEqual(payload["counts"]["external"]["species_confirmed"], 0)
        self.assertFalse(payload["scientific_boundary"]["phenotypes_read"])


def typed_row(isolate: str, st: str, alleles: list[int], cohort: str = "development", eligible: str = "true"):
    return {"isolate_id": isolate, "cohort": cohort, "eligible_for_modeling": eligible, "lineage_group": st,
            "mlst_alleles": ";".join(f"{locus}({value})" for locus, value in zip(ACHTMAN_LOCI, alleles))}


ST131 = [53, 40, 47, 13, 36, 28, 29]
ST9126 = [53, 40, 47, 13, 36, 28, 999]      # SLV of ST131
ST9999 = [53, 40, 47, 13, 36, 777, 999]     # SLV of ST9126 (DLV of ST131)
ST73 = [36, 24, 9, 13, 17, 11, 25]
ST95 = [37, 38, 19, 37, 17, 11, 26]         # DLV of ST73


class SlvGroupingTests(unittest.TestCase):
    def test_single_locus_variants_chain_into_one_named_group(self) -> None:
        rows = [typed_row("A", "ST131", ST131), typed_row("B", "ST131", ST131),
                typed_row("C", "ST9126", ST9126), typed_row("D", "ST9999", ST9999),
                typed_row("E", "ST73", ST73), typed_row("F", "ST95", ST95)]
        summary = MODULE.assign_slv_groups(rows)
        groups = {row["isolate_id"]: row["lineage_group"] for row in rows}
        self.assertEqual({groups[i] for i in "ABCD"}, {"SLV:ST131"})  # named by most frequent ST
        self.assertEqual(groups["E"], "SLV:ST73")
        self.assertEqual(groups["F"], "SLV:ST95")                     # two loci apart: separate
        self.assertEqual(summary["development_lineage_groups"], 3)

    def test_ties_name_by_lowest_st_and_unresolved_stays_blank(self) -> None:
        rows = [typed_row("A", "ST9126", ST9126), typed_row("B", "ST131", ST131),
                {"isolate_id": "U", "cohort": "development", "eligible_for_modeling": "true",
                 "lineage_group": "", "mlst_alleles": ""}]
        MODULE.assign_slv_groups(rows)
        self.assertEqual(rows[0]["lineage_group"], "SLV:ST131")
        self.assertEqual(rows[2]["lineage_group"], "")

    def test_external_and_ineligible_genomes_never_define_or_merge_groups(self) -> None:
        # E1 would bridge ST73 and ST95 if it were allowed to define groups.
        bridge = [36, 24, 9, 13, 17, 11, 26]  # SLV of ST73 and of ST95
        rows = [typed_row("A", "ST73", ST73), typed_row("B", "ST95", [36, 24, 9, 13, 17, 12, 26]),
                typed_row("E1", "ST5000", bridge, cohort="external"),
                typed_row("E2", "ST9126", ST9126, cohort="external"),
                typed_row("E3", "ST131", ST131, cohort="external"),
                typed_row("X", "ST4000", bridge, eligible="false")]
        MODULE.assign_slv_groups(rows)
        groups = {row["isolate_id"]: (row["lineage_group"], row["lineage_basis"]) for row in rows}
        self.assertNotEqual(groups["A"][0], groups["B"][0])
        self.assertEqual(groups["E1"], ("SLV:ST5000", "outside_development_groups"))
        self.assertEqual(groups["X"], ("SLV:ST4000", "outside_development_groups"))
        self.assertEqual(groups["E2"], ("SLV:ST9126", "outside_development_groups"))  # no ST131 in development
        self.assertEqual(groups["E3"], ("SLV:ST131", "outside_development_groups"))

    def test_external_st_joins_its_development_group_or_single_slv_neighbour(self) -> None:
        rows = [typed_row("A", "ST131", ST131), typed_row("E1", "ST131", ST131, cohort="external"),
                typed_row("E2", "ST9126", ST9126, cohort="external")]
        MODULE.assign_slv_groups(rows)
        self.assertEqual((rows[1]["lineage_group"], rows[1]["lineage_basis"]), ("SLV:ST131", "development_slv_graph"))
        self.assertEqual((rows[2]["lineage_group"], rows[2]["lineage_basis"]),
                         ("SLV:ST131", "slv_of_development_group"))

    def test_same_st_with_different_profiles_stops_the_run(self) -> None:
        rows = [typed_row("A", "ST131", ST131), typed_row("B", "ST131", ST9126)]
        with self.assertRaises(MODULE.FinalizeError):
            MODULE.assign_slv_groups(rows)
