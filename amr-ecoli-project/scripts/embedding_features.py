"""Isolate-level features from Evo 2 resistance-unit embeddings, fitted in-fold.

Specification (docs/FOUNDATION_MODEL_PROTOCOL.md revision 0.2, feature step),
fixed before any embedding is evaluated:

- Chromosomal panel: each of the 12 loci keeps its own strand-averaged
  embedding, reduced by PCA to PANEL_COMPONENTS components. A missing locus
  contributes zeros plus a missing indicator.
- Acquired determinants: an isolate carries a variable number, so their
  embeddings are averaged, reduced by PCA to ACQUIRED_COMPONENTS components,
  plus the count of acquired units (zeros when there are none).
- Every PCA is fitted on the rows passed to fit() only, i.e. inside the
  training part of a fold; transform() never refits.
- The layer is a hyperparameter chosen inside training folds.
- parts selects "panel", "acquired" or both; the single-part variants exist
  only for the diagnostic ablation that locates any Evo 2 gain.
"""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA

PANEL_COMPONENTS = 16
ACQUIRED_COMPONENTS = 32


class EmbeddingFeatureError(ValueError):
    """Raised when embeddings are missing or inconsistent."""


def load_layer(shard_root: Path, layer: str) -> tuple[dict[str, int], np.ndarray]:
    """Unique-sequence vectors for one layer, with each shard's checksum verified."""
    key = layer.replace(".", "_")
    index: dict[str, int] = {}
    blocks = []
    for shard in sorted(shard_root.glob("shard-*.npz")):
        recorded = shard.with_suffix(".sha256").read_text(encoding="utf-8").split()[0]
        if hashlib.sha256(shard.read_bytes()).hexdigest() != recorded:
            raise EmbeddingFeatureError(f"Embedding shard changed: {shard}")
        data = np.load(shard)
        for digest in data["hashes"]:
            index[str(digest)] = len(index)
        blocks.append(data[key].astype(np.float32))
    if not blocks:
        raise EmbeddingFeatureError(f"No embedding shards in {shard_root}")
    return index, np.concatenate(blocks)


def load_mapping(path: Path) -> dict[str, dict[str, list[str]]]:
    """isolate -> {"panel:<locus>": [hash], "acquired": [hashes]}."""
    mapping: dict[str, dict[str, list[str]]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row["sequence_sha256"]:
                continue
            unit = row["unit"] if row["unit"].startswith("panel:") else "acquired"
            mapping.setdefault(row["isolate_id"], {}).setdefault(unit, []).append(row["sequence_sha256"])
    return mapping


class EmbeddingFeaturizer:
    def __init__(self, mapping: dict[str, dict[str, list[str]]], index: dict[str, int], vectors: np.ndarray,
                 panel_loci: list[str], parts: tuple[str, ...] = ("panel", "acquired")):
        if not parts or not set(parts) <= {"panel", "acquired"}:
            raise EmbeddingFeatureError(f"parts must be drawn from panel and acquired, got {parts}")
        self.mapping, self.index, self.vectors = mapping, index, vectors
        self.panel_loci = panel_loci if "panel" in parts else []
        self.use_acquired = "acquired" in parts
        self.pcas: dict[str, PCA | None] = {}

    def _panel_block(self, isolates: list[str], locus: str) -> tuple[np.ndarray, np.ndarray]:
        rows = np.full(len(isolates), -1)
        for i, isolate in enumerate(isolates):
            hashes = self.mapping.get(isolate, {}).get(f"panel:{locus}", [])
            if hashes:
                if hashes[0] not in self.index:
                    raise EmbeddingFeatureError(f"No embedding for {isolate} {locus}")
                rows[i] = self.index[hashes[0]]
        present = rows >= 0
        block = np.zeros((len(isolates), self.vectors.shape[1]), dtype=np.float32)
        block[present] = self.vectors[rows[present]]
        return block, present

    def _acquired_block(self, isolates: list[str]) -> tuple[np.ndarray, np.ndarray]:
        block = np.zeros((len(isolates), self.vectors.shape[1]), dtype=np.float32)
        counts = np.zeros(len(isolates))
        for i, isolate in enumerate(isolates):
            hashes = self.mapping.get(isolate, {}).get("acquired", [])
            missing = [h for h in hashes if h not in self.index]
            if missing:
                raise EmbeddingFeatureError(f"No embedding for {len(missing)} acquired units of {isolate}")
            if hashes:
                block[i] = self.vectors[[self.index[h] for h in hashes]].mean(axis=0)
                counts[i] = len(hashes)
        return block, counts

    @staticmethod
    def _fit_pca(block: np.ndarray, present: np.ndarray, components: int, seed: int) -> PCA | None:
        rows = block[present]
        k = min(components, rows.shape[0] - 1, rows.shape[1])
        if k < 1:
            return None
        return PCA(n_components=k, svd_solver="randomized", random_state=seed).fit(rows)

    def fit(self, isolates: list[str], seed: int) -> "EmbeddingFeaturizer":
        self.pcas = {}
        for locus in self.panel_loci:
            block, present = self._panel_block(isolates, locus)
            self.pcas[locus] = self._fit_pca(block, present, PANEL_COMPONENTS, seed)
        if self.use_acquired:
            block, counts = self._acquired_block(isolates)
            self.pcas["acquired"] = self._fit_pca(block, counts > 0, ACQUIRED_COMPONENTS, seed)
        else:
            self.pcas["acquired"] = None
        return self

    def transform(self, isolates: list[str]) -> np.ndarray:
        if not self.pcas:
            raise EmbeddingFeatureError("transform() before fit()")
        parts = []
        for locus in self.panel_loci:
            block, present = self._panel_block(isolates, locus)
            pca = self.pcas[locus]
            width = PANEL_COMPONENTS
            reduced = np.zeros((len(isolates), width), dtype=np.float32)
            if pca is not None and present.any():
                reduced[present, :pca.n_components_] = pca.transform(block[present])
            parts += [reduced, (~present).astype(np.float32)[:, None]]
        if self.use_acquired:
            block, counts = self._acquired_block(isolates)
            pca = self.pcas["acquired"]
            reduced = np.zeros((len(isolates), ACQUIRED_COMPONENTS), dtype=np.float32)
            has = counts > 0
            if pca is not None and has.any():
                reduced[has, :pca.n_components_] = pca.transform(block[has])
            parts += [reduced, counts.astype(np.float32)[:, None]]
        return np.hstack(parts)

    def width(self) -> int:
        return len(self.panel_loci) * (PANEL_COMPONENTS + 1) + (ACQUIRED_COMPONENTS + 1) * self.use_acquired
