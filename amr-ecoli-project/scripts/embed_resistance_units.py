"""Embed deduplicated resistance units with Evo 2 (protocol revision 0.2).

Each unique sequence (from extract_resistance_loci.py collect) is embedded on
both strands; per candidate layer the token embeddings are mean-pooled and the
two strands averaged (the model is not reverse-complement invariant). The
cosine between forward and reverse-complement vectors is recorded per
sequence as the reverse-complement consistency diagnostic.

Work is split into shards of sequences sorted by hash; a shard is written
atomically with its SHA-256, so a preempted job resumes by skipping finished
shards. Genotype only; no phenotype, split or external membership is read.

Run from the repository root (GPU node):
    python scripts/embed_resistance_units.py --unique-fasta ... --output-root ... \\
        --model evo2_20b --shard-index 0 --shard-count 20
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import numpy as np

LAYERS = ("blocks.6.mlp.l3", "blocks.12.mlp.l3", "blocks.17.mlp.l3")
SEQUENCES_PER_SHARD = 500
MAX_UNIT_LENGTH = 20_000
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


class EmbeddingError(ValueError):
    """Raised when the embedding contract cannot be honoured."""


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def read_unique(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name, parts = None, []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts)))
                name, parts = line[1:], []
            elif line:
                parts.append(line)
    if name is not None:
        records.append((name, "".join(parts)))
    for digest, sequence in records:
        if hashlib.sha256(sequence.encode("ascii")).hexdigest() != digest:
            raise EmbeddingError(f"Unique sequence {digest[:12]} does not match its hash")
        if len(sequence) > MAX_UNIT_LENGTH:
            raise EmbeddingError(f"Unit {digest[:12]} is {len(sequence)} bp, above {MAX_UNIT_LENGTH}")
    return sorted(records)


def shards(records: list[tuple[str, str]]) -> list[list[tuple[str, str]]]:
    return [records[i:i + SEQUENCES_PER_SHARD] for i in range(0, len(records), SEQUENCES_PER_SHARD)]


def embed_shard(records: list[tuple[str, str]],
                embed: Callable[[str], dict[str, np.ndarray]]) -> dict[str, np.ndarray]:
    """Strand-averaged vectors per layer plus the per-sequence RC cosine."""
    vectors: dict[str, list[np.ndarray]] = {layer: [] for layer in LAYERS}
    rc_cosine: dict[str, list[float]] = {layer: [] for layer in LAYERS}
    for _, sequence in records:
        forward, reverse = embed(sequence), embed(reverse_complement(sequence))
        for layer in LAYERS:
            f, r = forward[layer].astype(np.float32), reverse[layer].astype(np.float32)
            vectors[layer].append(((f + r) / 2.0).astype(np.float16))
            rc_cosine[layer].append(float(f @ r / (np.linalg.norm(f) * np.linalg.norm(r) + 1e-12)))
    out = {"hashes": np.array([digest for digest, _ in records])}
    for layer in LAYERS:
        key = layer.replace(".", "_")
        out[key] = np.stack(vectors[layer])
        out[f"{key}__rc_cosine"] = np.array(rc_cosine[layer], dtype=np.float32)
    return out


def write_shard(path: Path, arrays: dict[str, np.ndarray]) -> str:
    partial = path.with_suffix(".partial.npz")
    np.savez(partial, **arrays)
    digest = hashlib.sha256(partial.read_bytes()).hexdigest()
    partial.rename(path)
    path.with_suffix(".sha256").write_text(f"{digest}  {path.name}\n", encoding="utf-8")
    return digest


def evo2_embedder(model_name: str) -> tuple[Callable[[str], dict[str, np.ndarray]], dict[str, str]]:
    import torch
    from evo2 import Evo2

    model = Evo2(model_name)

    def embed(sequence: str) -> dict[str, np.ndarray]:
        ids = torch.tensor(model.tokenizer.tokenize(sequence), dtype=torch.int, device="cuda").unsqueeze(0)
        with torch.inference_mode():
            _, layers = model(ids, return_embeddings=True, layer_names=list(LAYERS))
        return {layer: layers[layer].float().mean(dim=1).squeeze(0).cpu().numpy() for layer in LAYERS}

    import importlib.metadata as metadata
    info = {"gpu": torch.cuda.get_device_name(0), "torch": torch.__version__, "evo2": metadata.version("evo2"),
            "model": model_name}
    return embed, info


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--unique-fasta", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--model", default="evo2_20b")
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--max-shards", type=int, help="canary: stop after this many shards")
    parser.add_argument("--canary-sequences", type=int, help="canary: embed only the first N sequences of a shard")
    args = parser.parse_args()
    records = read_unique(args.unique_fasta)
    all_shards = shards(records)
    mine = [i for i in range(len(all_shards)) if i % args.shard_count == args.shard_index]
    if args.max_shards is not None:
        mine = mine[:args.max_shards]
    args.output_root.mkdir(parents=True, exist_ok=True)
    todo = [i for i in mine if not (args.output_root / f"shard-{i:05d}.npz").exists()]
    if not todo:
        print(json.dumps({"shards_assigned": len(mine), "todo": 0}))
        return 0
    import torch
    embed, info = evo2_embedder(args.model)
    started = time.time()
    bases = 0
    for i in todo:
        batch = all_shards[i][:args.canary_sequences] if args.canary_sequences else all_shards[i]
        arrays = embed_shard(batch, embed)
        write_shard(args.output_root / f"shard-{i:05d}.npz", arrays)
        bases += 2 * sum(len(sequence) for _, sequence in batch)
        print(json.dumps({"shard": i, "sequences": len(batch), "elapsed_s": round(time.time() - started, 1),
                          "bases_per_s": round(bases / max(time.time() - started, 1e-9), 1)}), flush=True)
    info.update({"peak_gpu_memory_gb": round(torch.cuda.max_memory_allocated() / 1e9, 2),
                 "shards_done": len(todo), "bases_embedded_both_strands": bases,
                 "wall_seconds": round(time.time() - started, 1), "layers": list(LAYERS),
                 "generated_at_utc": datetime.now(timezone.utc).isoformat()})
    (args.output_root / f"run-{args.shard_index:03d}-{int(started)}.json").write_text(
        json.dumps(info, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(info, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
