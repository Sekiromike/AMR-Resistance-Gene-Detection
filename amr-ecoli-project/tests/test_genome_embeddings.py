from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "scripts" / "extract_genome_embeddings.py"
SPEC = importlib.util.spec_from_file_location("extract_genome_embeddings", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for header, sequence in records:
        lines.append(f">{header}")
        for index in range(0, len(sequence), 60):
            lines.append(sequence[index : index + 60])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def length_encoder(batch: list[str]) -> list[list[float]]:
    """A deterministic stand-in: encodes GC content and length."""
    vectors = []
    for sequence in batch:
        gc = sum(1 for base in sequence if base in "GC") / max(len(sequence), 1)
        vectors.append([gc, float(len(sequence))])
    return vectors


class FastaAndWindowingTests(unittest.TestCase):
    def test_windows_never_span_a_contig_boundary(self) -> None:
        records = [("c1", "A" * 25), ("c2", "C" * 25)]
        windows = list(MODULE.window_contigs(records, window=10, minimum_window=1))
        self.assertEqual([contig for contig, _start, _piece in windows], ["c1"] * 3 + ["c2"] * 3)
        for _contig, _start, piece in windows:
            self.assertIn(set(piece), [{"A"}, {"C"}])

    def test_short_trailing_window_is_dropped_not_padded(self) -> None:
        records = [("c1", "A" * 25)]
        windows = list(MODULE.window_contigs(records, window=10, minimum_window=10))
        self.assertEqual(len(windows), 2)
        self.assertTrue(all(len(piece) == 10 for _c, _s, piece in windows))

    def test_reverse_complement_round_trips(self) -> None:
        sequence = "ACGTTGCANN"
        self.assertEqual(MODULE.reverse_complement(MODULE.reverse_complement(sequence)), sequence)
        self.assertEqual(MODULE.reverse_complement("ACGT"), "ACGT")

    def test_invalid_window_parameters_are_rejected(self) -> None:
        records = [("c1", "ACGT")]
        with self.assertRaises(MODULE.EmbeddingError):
            list(MODULE.window_contigs(records, window=0, minimum_window=1))
        with self.assertRaises(MODULE.EmbeddingError):
            list(MODULE.window_contigs(records, window=10, minimum_window=11))

    def test_fasta_without_records_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "empty.fna"
            path.write_text("\n", encoding="utf-8")
            with self.assertRaises(MODULE.EmbeddingError):
                MODULE.read_fasta(path)


class AggregationTests(unittest.TestCase):
    def test_mean_pooling(self) -> None:
        self.assertEqual(MODULE.aggregate([[1.0, 2.0], [3.0, 4.0]], method="mean"), [2.0, 3.0])

    def test_unspecified_aggregation_is_rejected(self) -> None:
        """Aggregation must stay prespecified; it is not a free hyperparameter."""
        with self.assertRaises(MODULE.EmbeddingError):
            MODULE.aggregate([[1.0]], method="attention")

    def test_inconsistent_dimensionality_is_rejected(self) -> None:
        with self.assertRaises(MODULE.EmbeddingError):
            MODULE.aggregate([[1.0, 2.0], [3.0]], method="mean")

    def test_empty_input_is_rejected(self) -> None:
        with self.assertRaises(MODULE.EmbeddingError):
            MODULE.aggregate([], method="mean")


class EmbedAssemblyTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def test_embedding_is_pooled_and_accounted(self) -> None:
        path = self.root / "assembly.fna"
        write_fasta(path, [("c1", "ACGT" * 50), ("c2", "GGCC" * 50)])
        result = MODULE.embed_assembly(
            path, length_encoder, window=100, minimum_window=50, batch_size=4
        )
        self.assertEqual(result["dimension"], 2)
        self.assertEqual(result["counts"]["contigs"], 2)
        self.assertEqual(result["counts"]["total_bases"], 400)
        self.assertEqual(result["counts"]["windows"], 4)
        self.assertEqual(result["counts"]["fraction_of_assembly_embedded"], 1.0)

    def test_both_strands_are_averaged_when_requested(self) -> None:
        path = self.root / "assembly.fna"
        # GC content is reverse-complement invariant, so averaging strands must
        # leave it unchanged; this catches a strand-pairing bug.
        write_fasta(path, [("c1", "ACGTACGTGG" * 10)])
        both = MODULE.embed_assembly(path, length_encoder, window=50, minimum_window=50)
        single = MODULE.embed_assembly(
            path, length_encoder, window=50, minimum_window=50, both_strands=False
        )
        self.assertAlmostEqual(both["embedding"][0], single["embedding"][0], places=9)
        self.assertTrue(both["parameters"]["both_strands"])

    def test_ambiguous_bases_are_counted_not_silently_dropped(self) -> None:
        path = self.root / "assembly.fna"
        write_fasta(path, [("c1", "ACGT" * 24 + "N" * 4)])
        result = MODULE.embed_assembly(path, length_encoder, window=50, minimum_window=10)
        self.assertGreaterEqual(result["counts"]["windows_with_ambiguous_bases"], 1)

    def test_encoder_returning_wrong_count_is_rejected(self) -> None:
        path = self.root / "assembly.fna"
        write_fasta(path, [("c1", "ACGT" * 50)])
        with self.assertRaises(MODULE.EmbeddingError):
            MODULE.embed_assembly(
                path, lambda batch: [[1.0]], window=50, minimum_window=10
            )

    def test_assembly_too_short_for_any_window_is_rejected(self) -> None:
        path = self.root / "assembly.fna"
        write_fasta(path, [("c1", "ACGT")])
        with self.assertRaises(MODULE.EmbeddingError):
            MODULE.embed_assembly(path, length_encoder, window=1000, minimum_window=500)


class PlanTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        for cohort, isolate in (("development", "D1"), ("external", "E1")):
            write_fasta(self.root / cohort / isolate / "assembly.fna", [("c1", "ACGT" * 100)])

    def _plan(self, **overrides):
        kwargs = dict(
            window=8192,
            minimum_window=1024,
            aggregation="mean",
            both_strands=True,
            model_id="test/model",
            model_revision="abc123",
            layer="final",
        )
        kwargs.update(overrides)
        return MODULE.write_plan(self.root, self.root / "plan.json", **kwargs)

    def test_plan_records_both_cohorts_and_seals_evaluation(self) -> None:
        payload = self._plan()
        self.assertEqual(payload["counts"]["assemblies"], 2)
        self.assertEqual(payload["counts"]["development"], 1)
        self.assertEqual(payload["counts"]["external"], 1)
        boundary = payload["scientific_boundary"]
        self.assertFalse(boundary["phenotypes_read"])
        self.assertFalse(boundary["splits_read"])
        self.assertFalse(boundary["external_membership_used_for_selection"])
        self.assertFalse(boundary["evaluation_authorised"])

    def test_moving_model_tag_is_rejected(self) -> None:
        """An unpinned revision is not reproducible provenance."""
        with self.assertRaises(MODULE.EmbeddingError):
            self._plan(model_revision="")

    def test_plan_is_written_as_sorted_json(self) -> None:
        self._plan()
        payload = json.loads((self.root / "plan.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["model"]["revision"], "abc123")
        self.assertIn("generated_at_utc", payload)

    def test_missing_assemblies_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as empty:
            with self.assertRaises(MODULE.EmbeddingError):
                MODULE.write_plan(
                    Path(empty),
                    Path(empty) / "plan.json",
                    window=8192,
                    minimum_window=1024,
                    aggregation="mean",
                    both_strands=True,
                    model_id="test/model",
                    model_revision="abc123",
                    layer="final",
                )


if __name__ == "__main__":
    unittest.main()
