from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

import numpy as np

from scripts import embedding_features as MODULE

LOCI = ["gyrA", "parC"]


def setup(n_isolates=60, dim=24, seed=0):
    rng = np.random.default_rng(seed)
    hashes = [f"h{i:03d}" for i in range(40)]
    vectors = rng.normal(size=(40, dim)).astype(np.float32)
    index = {h: i for i, h in enumerate(hashes)}
    mapping = {}
    for i in range(n_isolates):
        units = {"panel:gyrA": [hashes[i % 10]], "acquired": [hashes[20 + (i % 7)], hashes[30 + (i % 5)]]}
        if i % 6 != 0:  # parC missing for every sixth isolate
            units["panel:parC"] = [hashes[10 + (i % 8)]]
        if i % 9 == 0:
            del units["acquired"]
        mapping[f"I{i}"] = units
    return mapping, index, vectors


class FeaturizerTests(unittest.TestCase):
    def test_shape_missing_indicators_and_counts(self) -> None:
        mapping, index, vectors = setup()
        isolates = [f"I{i}" for i in range(60)]
        f = MODULE.EmbeddingFeaturizer(mapping, index, vectors, LOCI).fit(isolates, seed=1)
        X = f.transform(isolates)
        self.assertEqual(X.shape, (60, f.width()))
        parc_missing = X[:, 2 * (MODULE.PANEL_COMPONENTS + 1) - 1]
        self.assertEqual(parc_missing.tolist(), [1.0 if i % 6 == 0 else 0.0 for i in range(60)])
        self.assertTrue(np.all(X[parc_missing == 1, MODULE.PANEL_COMPONENTS + 1:2 * MODULE.PANEL_COMPONENTS + 1] == 0))
        self.assertEqual(X[:, -1].tolist(), [0.0 if i % 9 == 0 else 2.0 for i in range(60)])

    def test_pca_is_fitted_only_on_training_isolates(self) -> None:
        mapping, index, vectors = setup()
        train, test = [f"I{i}" for i in range(40)], [f"I{i}" for i in range(40, 60)]
        a = MODULE.EmbeddingFeaturizer(mapping, index, vectors, LOCI).fit(train, seed=1)
        before = a.transform(test)
        # Changing test isolates' vectors must not change the fitted PCA (only their projections).
        mean_before = a.pcas["gyrA"].mean_.copy()
        a.transform(train + test)
        np.testing.assert_array_equal(a.pcas["gyrA"].mean_, mean_before)
        b = MODULE.EmbeddingFeaturizer(mapping, index, vectors, LOCI).fit(train + test, seed=1)
        self.assertFalse(np.allclose(before, b.transform(test)))

    def test_transform_requires_fit_and_embeddings_must_exist(self) -> None:
        mapping, index, vectors = setup()
        with self.assertRaises(MODULE.EmbeddingFeatureError):
            MODULE.EmbeddingFeaturizer(mapping, index, vectors, LOCI).transform(["I1"])
        mapping["I1"]["panel:gyrA"] = ["unknown"]
        with self.assertRaises(MODULE.EmbeddingFeatureError):
            MODULE.EmbeddingFeaturizer(mapping, index, vectors, LOCI).fit(["I1", "I2", "I3"], seed=1)

    def test_shards_are_checksum_verified(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            shard = root / "shard-00000.npz"
            np.savez(shard, hashes=np.array(["a", "b"]), blocks_6_mlp_l3=np.ones((2, 3), dtype=np.float16))
            shard.with_suffix(".sha256").write_text(hashlib.sha256(shard.read_bytes()).hexdigest() + "  x\n")
            index, vectors = MODULE.load_layer(root, "blocks.6.mlp.l3")
            self.assertEqual(index, {"a": 0, "b": 1})
            self.assertEqual(vectors.dtype, np.float32)
            shard.with_suffix(".sha256").write_text("0" * 64 + "  x\n")
            with self.assertRaises(MODULE.EmbeddingFeatureError):
                MODULE.load_layer(root, "blocks.6.mlp.l3")


class AblationPartsTests(unittest.TestCase):
    def test_panel_only_and_acquired_only_widths(self) -> None:
        mapping, index, vectors = setup()
        isolates = [f"I{i}" for i in range(60)]
        panel = MODULE.EmbeddingFeaturizer(mapping, index, vectors, LOCI, ("panel",)).fit(isolates, seed=1)
        acquired = MODULE.EmbeddingFeaturizer(mapping, index, vectors, LOCI, ("acquired",)).fit(isolates, seed=1)
        self.assertEqual(panel.transform(isolates).shape[1], len(LOCI) * (MODULE.PANEL_COMPONENTS + 1))
        self.assertEqual(acquired.transform(isolates).shape[1], MODULE.ACQUIRED_COMPONENTS + 1)
        with self.assertRaises(MODULE.EmbeddingFeatureError):
            MODULE.EmbeddingFeaturizer(mapping, index, vectors, LOCI, ("genome",))
