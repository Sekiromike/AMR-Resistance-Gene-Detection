"""End-to-end check of the post-genomics chain on synthetic data.

finalize_genome_analyses -> assess_lineage_adequacy -> build_isolate_metadata
-> construct_research_cohort (NCBI + JARBS sources, mic_regression)
-> validate_research_cohort.

Covers a registered development assembly, a development raw-read assembly, a
locked external raw-read assembly and an external isolate excluded by the
frozen membership, so every stage sees both genome sources and both cohorts.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "tests"))


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FINALIZE = load("finalize_genome_analyses")
ADEQUACY = load("assess_lineage_adequacy")
METADATA = load("build_isolate_metadata")
CONSTRUCT = load("construct_research_cohort")
VALIDATE = load("validate_research_cohort")
JARBS = load("normalize_jarbs_ast")

from test_finalize_genome_analyses import CONFIG, make_genome  # noqa: E402

REQUEST_FIELDS = ["sequence_id", "cohort", "isolate_id", "biosample_accession", "source_type",
                  "sequence_accession"]
GENOMES = [
    # isolate, cohort, source, accession, biosample, ST
    ("PDT1", "development", "registered_assembly", "GCA_000000001.1", "SAMN00000001", "131"),
    ("PDT2", "development", "raw_reads", "SRR0000002", "SAMN00000002", "73"),
    ("JB1", "external", "raw_reads", "DRR0000001", "SAMD00000001", "10"),
    ("JB2", "external", "raw_reads", "DRR0000002", "SAMD00000002", "69"),
]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]], delimiter: str = "\t") -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


def ncbi_row(record: str, isolate: str, biosample: str, accession: str, drug: str,
             sign: str, value: str) -> dict[str, str]:
    return {
        "source_dataset": "NCBI Pathogen Detection AST", "source_release": "snapshot",
        "source_fingerprint": "a" * 64, "source_row_number": "2", "source_ast_record_id": record,
        "isolate_id": isolate, "biosample_accession": biosample, "assembly_accession": accession,
        "scientific_name": "Escherichia coli", "antibiotic": drug, "ast_measurement_type": "MIC",
        "raw_measurement_sign": sign, "measurement_sign": sign, "raw_ast_value": value,
        "ast_value": value, "ast_unit": "mg/L", "ast_method": "broth microdilution",
        "ast_platform": "", "ast_testing_date": "", "mic_eligible": "true",
        "mic_ineligibility_reasons": "", "collection_date": "2020", "geo_loc_name": "USA",
        "country": "USA", "country_source": "geo_loc_name_prefix",
        "bioproject_accession": "PRJNA1", "isolation_source": "urine",
    }


class PostGenomicsPipelineTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def test_chain_produces_a_valid_mic_cohort_from_both_sources(self) -> None:
        r = self.root
        analysis = r / "analysis"
        for isolate, cohort, _, _, _, st in GENOMES:
            make_genome(analysis, cohort, isolate, st=st)
        requests = write_tsv(r / "requests.tsv", REQUEST_FIELDS, [
            {"sequence_id": f"{cohort}::{isolate}", "cohort": cohort, "isolate_id": isolate,
             "biosample_accession": biosample, "source_type": source, "sequence_accession": accession}
            for isolate, cohort, source, accession, biosample, _ in GENOMES
        ])
        (r / "config.json").write_text(json.dumps(CONFIG), encoding="utf-8")

        # 1. Genome manifest from per-genome analysis outputs.
        FINALIZE.finalize(requests, analysis, r / "config.json", r / "genomes.tsv", r / "genomes.json")
        with (r / "genomes.tsv").open(encoding="utf-8") as handle:
            genomes = {row["isolate_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(genomes["PDT2"]["assembly_accession"], "")
        self.assertTrue(all(row["eligible_for_modeling"] == "true" for row in genomes.values()))

        # 2. Lineage adequacy (thresholds relaxed only because the fixture is tiny).
        attributes = write_tsv(r / "attributes.tsv",
                               ["isolate_id", "cohort", "specimen_source", "intended_use_population"], [
            # PDT2's specimen source is absent in source: kept, never inferred.
            {"isolate_id": i, "cohort": c, "specimen_source": "" if i == "PDT2" else "urine",
             "intended_use_population": "human_clinical"}
            for i, c, *_ in GENOMES
        ])
        adequacy = ADEQUACY.assess(r / "genomes.tsv", attributes, {
            "min_typing_coverage": 0.95, "max_largest_lineage_fraction": 1.0,
            "min_inverse_simpson": 1, "folds": 5})
        self.assertEqual(adequacy["outcome"], "keep_st")
        self.assertEqual(adequacy["eligible_development_genomes"], 2)

        # 3. Isolate metadata from split cluster files, membership and pairs.
        dev_clusters = write_tsv(r / "dev-clusters.tsv", ["sequence_id", "genomic_cluster"], [
            {"sequence_id": "development::PDT1", "genomic_cluster": "development-cluster-000001"},
            {"sequence_id": "development::PDT2", "genomic_cluster": "development-cluster-000002"},
        ])
        ext_clusters = write_tsv(r / "ext-clusters.tsv", ["sequence_id", "genomic_cluster"], [
            {"sequence_id": "external::JB1", "genomic_cluster": "external-cluster-000001"},
            {"sequence_id": "external::JB2", "genomic_cluster": "external-cluster-000002"},
        ])
        membership = write_tsv(r / "membership.tsv", ["sequence_id", "membership"], [
            {"sequence_id": "external::JB1", "membership": "external_locked"},
            {"sequence_id": "external::JB2", "membership": "external_excluded"},
        ])
        pairs = write_tsv(r / "pairs.tsv", ["development_sequence_id", "external_sequence_id", "ani_percent",
                                             "aligned_fraction_development", "aligned_fraction_external"], [])
        METADATA.build_isolate_metadata(
            requests, [dev_clusters, ext_clusters], membership, r / "genomes.tsv", attributes, pairs,
            r / "metadata.csv", r / "metadata-exclusions.csv", r / "ledger.csv", r / "metadata.json")
        with (r / "metadata.csv").open(encoding="utf-8") as handle:
            metadata = {row["isolate_id"]: row for row in csv.DictReader(handle)}
        self.assertEqual(set(metadata), {"PDT1", "PDT2", "JB1"})
        self.assertEqual(metadata["JB1"]["evaluation_split"], "external")
        self.assertEqual(metadata["PDT2"]["genome_source"], "raw_reads")

        # 4. AST: NCBI development rows plus JARBS rows from the real normalizer.
        ncbi = [
            ncbi_row("n1", "PDT1", "SAMN00000001", "GCA_000000001.1", "ciprofloxacin", "=", "0.25"),
            ncbi_row("n2", "PDT2", "SAMN00000002", "", "ciprofloxacin", ">", "4"),
        ]
        write_tsv(r / "ncbi.csv", list(ncbi[0]), ncbi, delimiter=",")
        jarbs_table = pd.DataFrame([
            {"isolate_id": isolate, "genome_species": "Escherichia coli", "biosample_accession": biosample,
             "collection_date_raw": "2019", "geo_loc_name_raw": "Japan", "isolation_source_raw": "urine",
             "ceftriaxone_mic_raw": "", "ciprofloxacin_mic_raw": "0.5", "gentamicin_mic_raw": ""}
            for isolate, cohort, _, _, biosample, _ in GENOMES if cohort == "external"
        ])
        jarbs_table.to_csv(r / "jarbs.csv", index=False)
        vocabulary = ROOT / "config" / "insdc_geo_loc_name_vocabulary.tsv"
        frame, _ = JARBS.normalize_jarbs(r / "jarbs.csv", vocabulary, "test")
        frame.to_csv(r / "jarbs-normalized.csv", index=False, lineterminator="\n")

        payload = CONSTRUCT.run_construction(
            [r / "ncbi.csv", r / "jarbs-normalized.csv"], r / "metadata.csv", r / "cohort.csv",
            r / "cohort-exclusions.csv", r / "cohort.json", r / "genomes.tsv", profile="mic_regression")
        self.assertEqual(len(payload["sources"]["ast_sources"]), 2)
        cohort = pd.read_csv(r / "cohort.csv", dtype=str, keep_default_na=False)
        included = set(zip(cohort["isolate_id"], cohort["antibiotic"]))
        self.assertEqual(included, {("PDT1", "ciprofloxacin"), ("PDT2", "ciprofloxacin"),
                                    ("JB1", "ciprofloxacin")})
        self.assertEqual(cohort.loc[cohort["isolate_id"] == "PDT2", "specimen_source"].tolist(), [""])
        # Raw-read genomes (no registered assembly) enter the cohort.
        self.assertEqual(set(cohort.loc[cohort["genome_source"] == "raw_reads", "isolate_id"]), {"PDT2", "JB1"})
        exclusions = pd.read_csv(r / "cohort-exclusions.csv", dtype=str, keep_default_na=False)
        self.assertTrue(set(exclusions["isolate_id"]) <= {"JB2"})

        # 5. Validation under the primary endpoint; external censoring stays sealed.
        report = VALIDATE.validate_cohort(pd.read_csv(r / "cohort.csv"), min_per_category=1,
                                          profile="mic_regression")
        self.assertEqual(report["failures"], [])
        for record in report["censoring"]:
            if record["evaluation_split"] == "external":
                self.assertEqual(record["fraction_censored"], "sealed")


if __name__ == "__main__":
    unittest.main()
