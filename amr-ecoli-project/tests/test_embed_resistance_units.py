from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

import numpy as np

from scripts import embed_resistance_units as MODULE


def fake_embed(sequence: str) -> dict[str, np.ndarray]:
    """Deterministic, strand-sensitive stand-in for Evo 2."""
    counts = np.array([sequence.count(b) for b in "ACGT"], dtype=np.float32)
    first = np.array([float("ACGT".index(sequence[0]))], dtype=np.float32)
    return {layer: np.concatenate([counts, first * (i + 1)]) for i, layer in enumerate(MODULE.LAYERS)}


def fasta(path: Path, sequences: list[str]) -> None:
    path.write_text("".join(f">{hashlib.sha256(s.encode()).hexdigest()}\n{s}\n" for s in sequences), encoding="utf-8")


class EmbedTests(unittest.TestCase):
    def test_strands_are_averaged_and_rc_consistency_recorded(self) -> None:
        sequence = "AACGTTTG"
        out = MODULE.embed_shard([("h", sequence)], fake_embed)
        f, r = fake_embed(sequence), fake_embed(MODULE.reverse_complement(sequence))
        layer = "blocks_6_mlp_l3"
        np.testing.assert_allclose(out[layer][0], ((f["blocks.6.mlp.l3"] + r["blocks.6.mlp.l3"]) / 2).astype(np.float16))
        self.assertEqual(out[layer].dtype, np.float16)
        self.assertLess(out[f"{layer}__rc_cosine"][0], 1.0)
        # A reverse-complement palindrome is embedded identically on both strands.
        palindrome = MODULE.embed_shard([("p", "ACGT")], fake_embed)
        self.assertAlmostEqual(float(palindrome[f"{layer}__rc_cosine"][0]), 1.0, places=5)

    def test_unique_fasta_is_hash_checked_and_length_capped(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "u.fna"
            fasta(path, ["ACGT" * 10, "GGCC" * 5])
            records = MODULE.read_unique(path)
            self.assertEqual(records, sorted(records))
            path.write_text(">deadbeef\nACGT\n", encoding="utf-8")
            with self.assertRaises(MODULE.EmbeddingError):
                MODULE.read_unique(path)
            fasta(path, ["A" * (MODULE.MAX_UNIT_LENGTH + 1)])
            with self.assertRaises(MODULE.EmbeddingError):
                MODULE.read_unique(path)

    def test_shards_are_written_atomically_with_checksums(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            shard = Path(temporary) / "shard-00000.npz"
            arrays = MODULE.embed_shard([("h1", "AACG"), ("h2", "TTGA")], fake_embed)
            digest = MODULE.write_shard(shard, arrays)
            self.assertTrue(shard.exists())
            self.assertFalse(shard.with_suffix(".partial.npz").exists())
            self.assertEqual(hashlib.sha256(shard.read_bytes()).hexdigest(), digest)
            loaded = np.load(shard)
            self.assertEqual(list(loaded["hashes"]), ["h1", "h2"])

    def test_sharding_covers_every_sequence_once(self) -> None:
        records = [(f"{i:04d}", "ACGT") for i in range(1234)]
        parts = MODULE.shards(records)
        self.assertEqual(sum(len(p) for p in parts), 1234)
        self.assertEqual(len(parts), 3)
