from __future__ import annotations

import csv
import random
import tempfile
import unittest
from pathlib import Path

from scripts import extract_resistance_loci as MODULE

random.seed(4)
CONTIG = "".join(random.choice("ACGT") for _ in range(5000))
HEADER = ["Name", "Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol", "Element name",
          "Scope", "Type", "Subtype", "Class", "Subclass", "Method"]


def hit(**kw):
    base = {"contig": "c1", "pident": 99.0, "qstart": 1, "qend": 1000, "sstart": 2001, "send": 3000,
            "bitscore": 1800.0, "qlen": 1000}
    base.update(kw)
    return base


class PanelTests(unittest.TestCase):
    def test_forward_hit_is_extended_to_the_full_reference_extent(self) -> None:
        sequence, status = MODULE.panel_unit(hit(qstart=11, qend=1000, sstart=2011), {"c1": CONTIG})
        self.assertEqual(status, "found")
        self.assertEqual(sequence, CONTIG[2000:3000])

    def test_reverse_hit_is_oriented_to_the_reference(self) -> None:
        sequence, status = MODULE.panel_unit(hit(sstart=3000, send=2001), {"c1": CONTIG})
        self.assertEqual(status, "found")
        self.assertEqual(sequence, MODULE.reverse_complement(CONTIG[2000:3000]))

    def test_thresholds_and_contig_ends_leave_the_unit_missing(self) -> None:
        self.assertIsNone(MODULE.panel_unit(hit(pident=79.9), {"c1": CONTIG})[0])
        self.assertIsNone(MODULE.panel_unit(hit(qend=899), {"c1": CONTIG})[0])
        seq, status = MODULE.panel_unit(hit(qstart=51, sstart=10, send=959), {"c1": CONTIG})
        self.assertIsNone(seq)
        self.assertEqual(status, "truncated_at_contig_end")

    def test_best_hit_per_locus_by_bitscore(self) -> None:
        tsv = "\n".join([
            "gyrA\tc1\t99\t1000\t1\t1000\t1\t1000\t0\t1500\t1000",
            "gyrA\tc2\t99\t1000\t1\t1000\t1\t1000\t0\t1900\t1000",
        ])
        self.assertEqual(MODULE.best_panel_hits(tsv)["gyrA"]["contig"], "c2")

    def test_reference_panel_takes_upstream_on_the_gene_strand(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            gff = Path(temporary) / "ref.gff"
            lines = ["##gff-version 3"]
            for i, name in enumerate(MODULE.PANEL):
                start = 400 + i * 300
                strand = "-" if name == "gyrB" else "+"
                lines.append(f"c1\tRefSeq\tgene\t{start}\t{start + 99}\t.\t{strand}\t.\tID=g{i};Name={name};gene={name}")
            gff.write_text("\n".join(lines) + "\n", encoding="utf-8")
            units = dict(MODULE.panel_from_reference({"c1": CONTIG}, gff))
        self.assertEqual(units["gyrA"], CONTIG[99:499])            # 300 upstream + 100 gene
        gyrb_start = 400 + 300
        self.assertEqual(units["gyrB"], MODULE.reverse_complement(CONTIG[gyrb_start - 1:gyrb_start + 99 + 300]))


class AcquiredTests(unittest.TestCase):
    def test_complete_acquired_genes_with_flanks_and_strand(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            report = Path(temporary) / "amrfinder.tsv"
            rows = [
                ["-", "NA", "c1", "2001", "2900", "-", "blaCTX-M-15", "x", "core", "AMR", "AMR", "BETA-LACTAM", "C", "EXACTX"],
                ["-", "NA", "c1", "100", "400", "+", "blaTEM-1", "x", "core", "AMR", "AMR", "BETA-LACTAM", "C", "PARTIALX"],
                ["-", "NA", "c1", "3000", "3500", "+", "gyrA_S83L", "x", "core", "AMR", "POINT", "QUINOLONE", "Q", "POINTX"],
                ["-", "NA", "c1", "100", "600", "+", "blaCTX-M-15", "x", "core", "AMR", "AMR", "BETA-LACTAM", "C", "BLASTX"],
            ]
            report.write_text("\n".join("\t".join(r) for r in [HEADER, *rows]) + "\n", encoding="utf-8")
            units = MODULE.acquired_units(report, {"c1": CONTIG})
        self.assertEqual([u[0] for u in units], ["blaCTX-M-15", "blaCTX-M-15#2"])
        self.assertEqual(units[0][1], MODULE.reverse_complement(CONTIG[1500:3400]))
        self.assertEqual(units[1][1], CONTIG[0:1100])               # left flank clipped at contig start


class CollectTests(unittest.TestCase):
    def test_identical_sequences_are_embedded_once_and_hashes_are_checked(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            table = Path(temporary) / "units.tsv"
            rows = []
            for isolate in ("A", "B"):
                seq = CONTIG[:900]
                rows.append({"isolate_id": isolate, "cohort": "development", "unit": "panel:gyrA", "status": "found",
                             "length": "900", "sequence_sha256": MODULE.sha256_text(seq), "sequence": seq})
            rows.append({"isolate_id": "B", "cohort": "development", "unit": "panel:parC", "status": "not_found",
                         "length": "0", "sequence_sha256": "", "sequence": ""})
            with table.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=MODULE.UNIT_FIELDS, delimiter="\t")
                writer.writeheader(); writer.writerows(rows)
            unique, mapping = MODULE.collect([table])
            self.assertEqual(len(unique), 1)
            self.assertEqual(len(mapping), 3)
            rows[0]["sequence"] = CONTIG[1:901]
            with table.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=MODULE.UNIT_FIELDS, delimiter="\t")
                writer.writeheader(); writer.writerows(rows)
            with self.assertRaises(MODULE.LocusError):
                MODULE.collect([table])
