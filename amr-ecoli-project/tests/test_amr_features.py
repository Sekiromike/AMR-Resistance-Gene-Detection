from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from scripts import build_amr_features as MODULE

HEADER = ["Name", "Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol", "Element name",
          "Scope", "Type", "Subtype", "Class", "Subclass", "Method"]


def hit(symbol, type_="AMR", subtype="AMR", klass="BETA-LACTAM", method="EXACTX"):
    return ["-", "NA", "c1", "1", "100", "+", symbol, symbol, "core", type_, subtype, klass, "CEPHALOSPORIN", method]


def write_genome(root: Path, cohort: str, isolate: str, hits: list[list[str]]) -> Path:
    directory = root / cohort / isolate
    (directory / "amrfinder").mkdir(parents=True)
    report = directory / "amrfinder" / "amrfinder.tsv"
    report.write_text("\n".join("\t".join(r) for r in [HEADER, *hits]) + "\n", encoding="utf-8")
    digest = hashlib.sha256(report.read_bytes()).hexdigest()
    (directory / "artifact.sha256").write_text(f"{digest}  ./amrfinder/amrfinder.tsv\n", encoding="utf-8")
    return report


def write_manifest(root: Path, rows: list[tuple[str, str, str]]) -> Path:
    path = root / "genomes.tsv"
    lines = ["isolate_id\tcohort\teligible_for_modeling"] + ["\t".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class AmrFeatureTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def test_only_complete_amr_calls_are_features(self) -> None:
        write_genome(self.root, "development", "D1", [
            hit("blaCTX-M-15"),
            hit("gyrA_S83L", subtype="POINT", klass="QUINOLONE", method="POINTX"),
            hit("aac(3)-IId", method="PARTIAL_CONTIG_ENDX"),
            hit("blaTEM-1", method="PARTIALX"),
            hit("tet(A)", method="INTERNAL_STOP"),
            hit("fdeC", type_="VIRULENCE", subtype="VIRULENCE"),
            hit("merA", type_="STRESS", subtype="METAL"),
            hit("blaCTX-M-15", method="BLASTX"),
        ])
        write_genome(self.root, "external", "E1", [])
        manifest = write_manifest(self.root, [("D1", "development", "true"), ("E1", "external", "true"),
                                              ("X1", "development", "false")])
        rows, summary = MODULE.extract(manifest, self.root)
        self.assertEqual([r["element_symbol"] for r in rows], ["aac(3)-IId", "blaCTX-M-15", "gyrA_S83L"])
        self.assertEqual(summary["genomes"], 2)
        self.assertEqual(summary["other_counts"]["genomes_without_amr_elements"], 1)
        self.assertEqual(summary["other_counts"]["excluded_type_virulence"], 1)
        self.assertEqual(summary["other_counts"]["duplicate_calls_collapsed"], 1)

    def test_altered_report_is_refused(self) -> None:
        report = write_genome(self.root, "development", "D1", [hit("blaCTX-M-15")])
        report.write_text(report.read_text(encoding="utf-8") + "\t".join(hit("blaTEM-1")) + "\n", encoding="utf-8")
        manifest = write_manifest(self.root, [("D1", "development", "true")])
        with self.assertRaisesRegex(MODULE.FeatureError, "changed"):
            MODULE.extract(manifest, self.root)

    def test_missing_report_is_refused(self) -> None:
        manifest = write_manifest(self.root, [("D9", "development", "true")])
        with self.assertRaisesRegex(MODULE.FeatureError, "Missing"):
            MODULE.extract(manifest, self.root)
