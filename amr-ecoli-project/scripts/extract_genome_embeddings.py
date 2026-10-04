"""Extract sequence foundation-model embeddings from frozen genome assemblies.

This module never reads a phenotype, a category, a split, or external cohort
membership.  It consumes assemblies and emits fixed-length representations plus
provenance, under the windowing and aggregation rules prespecified in
docs/FOUNDATION_MODEL_PROTOCOL.md.

Computing a representation is not crediting a model: nothing produced here may
be evaluated against phenotype until the locked baselines in rungs 1-4 of
docs/STUDY_PROTOCOL.md are complete.

Run from the repository root:
    python scripts/extract_genome_embeddings.py plan --assemblies ... --output ...
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator, Sequence


VALID_BASES = frozenset("ACGT")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
AGGREGATIONS = ("mean",)


class EmbeddingError(ValueError):
    """Raised when the embedding contract cannot be verified."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Read a FASTA assembly, preserving contig order and uppercasing bases."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            else:
                if header is None:
                    raise EmbeddingError(f"FASTA content precedes any header: {path}")
                chunks.append(line.upper())
    if header is not None:
        records.append((header, "".join(chunks)))
    if not records:
        raise EmbeddingError(f"Assembly contains no sequence records: {path}")
    return records


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def window_contigs(
    records: Sequence[tuple[str, str]],
    *,
    window: int,
    minimum_window: int,
) -> Iterator[tuple[str, int, str]]:
    """Yield non-overlapping windows in the given contig order.

    Windows never span a contig boundary, because adjacency between contigs of a
    draft assembly is not evidence of adjacency in the genome. A trailing window
    shorter than `minimum_window` is dropped rather than padded, so padding is
    never embedded as if it were sequence.
    """
    if window < 1:
        raise EmbeddingError("Window length must be positive")
    if not 0 < minimum_window <= window:
        raise EmbeddingError("Minimum window must be positive and at most the window length")
    for contig, sequence in records:
        for start in range(0, len(sequence), window):
            piece = sequence[start : start + window]
            if len(piece) < minimum_window:
                continue
            yield contig, start, piece


def aggregate(vectors: Sequence[Sequence[float]], *, method: str) -> list[float]:
    """Pool per-window vectors into one genome-level representation."""
    if method not in AGGREGATIONS:
        raise EmbeddingError(
            f"Unsupported aggregation {method!r}; prespecified options are {list(AGGREGATIONS)}"
        )
    if not vectors:
        raise EmbeddingError("Cannot aggregate an empty set of window vectors")
    width = len(vectors[0])
    if width == 0:
        raise EmbeddingError("Window vectors must be non-empty")
    if any(len(vector) != width for vector in vectors):
        raise EmbeddingError("Window vectors have inconsistent dimensionality")
    return [sum(vector[index] for vector in vectors) / len(vectors) for index in range(width)]


def embed_assembly(
    path: Path,
    encoder: Callable[[list[str]], list[list[float]]],
    *,
    window: int,
    minimum_window: int,
    aggregation: str = "mean",
    both_strands: bool = True,
    batch_size: int = 8,
) -> dict[str, Any]:
    """Embed one assembly into a single pooled vector plus window accounting.

    `encoder` maps a batch of sequences to a batch of vectors. It is injected so
    the windowing, strand, and pooling contract is testable without a GPU or
    model weights.
    """
    if batch_size < 1:
        raise EmbeddingError("Batch size must be positive")
    records = read_fasta(path)
    windows = list(window_contigs(records, window=window, minimum_window=minimum_window))
    if not windows:
        raise EmbeddingError(f"Assembly yielded no windows at window={window}: {path}")

    sequences: list[str] = []
    for _contig, _start, piece in windows:
        sequences.append(piece)
        if both_strands:
            # The architecture is not assumed reverse-complement invariant, so
            # both strands are embedded and averaged. The choice is recorded in
            # the manifest rather than left implicit.
            sequences.append(reverse_complement(piece))

    vectors: list[list[float]] = []
    for index in range(0, len(sequences), batch_size):
        batch = sequences[index : index + batch_size]
        encoded = encoder(batch)
        if len(encoded) != len(batch):
            raise EmbeddingError("Encoder returned a different number of vectors than inputs")
        vectors.extend(list(map(float, vector)) for vector in encoded)

    if both_strands:
        paired: list[list[float]] = []
        for index in range(0, len(vectors), 2):
            forward, reverse = vectors[index], vectors[index + 1]
            if len(forward) != len(reverse):
                raise EmbeddingError("Forward and reverse vectors differ in dimensionality")
            paired.append([(a + b) / 2 for a, b in zip(forward, reverse, strict=True)])
        vectors = paired

    pooled = aggregate(vectors, method=aggregation)
    total_bases = sum(len(sequence) for _header, sequence in records)
    embedded_bases = sum(len(piece) for _contig, _start, piece in windows)
    ambiguous = sum(
        1 for _contig, _start, piece in windows if not set(piece).issubset(VALID_BASES)
    )
    return {
        "embedding": pooled,
        "dimension": len(pooled),
        "counts": {
            "contigs": len(records),
            "total_bases": total_bases,
            "windows": len(windows),
            "embedded_bases": embedded_bases,
            "windows_with_ambiguous_bases": ambiguous,
            "fraction_of_assembly_embedded": round(embedded_bases / total_bases, 6)
            if total_bases
            else 0.0,
        },
        "parameters": {
            "window": window,
            "minimum_window": minimum_window,
            "aggregation": aggregation,
            "both_strands": both_strands,
        },
    }


def write_plan(
    assemblies_root: Path,
    output_path: Path,
    *,
    window: int,
    minimum_window: int,
    aggregation: str,
    both_strands: bool,
    model_id: str,
    model_revision: str,
    layer: str,
) -> dict[str, Any]:
    """Freeze the extraction request before any GPU time is spent."""
    if aggregation not in AGGREGATIONS:
        raise EmbeddingError(f"Unsupported aggregation: {aggregation}")
    if not model_revision:
        raise EmbeddingError(
            "A pinned model revision is required; a moving tag is not reproducible provenance."
        )
    assemblies = sorted(assemblies_root.glob("*/*/assembly.fna"))
    if not assemblies:
        raise EmbeddingError(f"No assemblies found under {assemblies_root}")
    entries = [
        {
            "cohort": path.parent.parent.name,
            "isolate_id": path.parent.name,
            "assembly_path": path.as_posix(),
        }
        for path in assemblies
    ]
    payload = {
        "schema_version": "1.0.0",
        "operation": "plan_genome_embeddings",
        "generated_at_utc": utc_now(),
        "model": {"id": model_id, "revision": model_revision, "layer": layer},
        "parameters": {
            "window": window,
            "minimum_window": minimum_window,
            "aggregation": aggregation,
            "both_strands": both_strands,
        },
        "counts": {
            "assemblies": len(entries),
            "development": sum(1 for row in entries if row["cohort"] == "development"),
            "external": sum(1 for row in entries if row["cohort"] == "external"),
        },
        "assemblies": entries,
        "scientific_boundary": {
            "phenotypes_read": False,
            "splits_read": False,
            "external_membership_used_for_selection": False,
            "evaluation_authorised": False,
            "note": (
                "Phenotype-free feature computation. Evaluation against phenotype is "
                "blocked until locked baselines are complete and the pretraining-overlap "
                "audit is attached. See docs/FOUNDATION_MODEL_PROTOCOL.md."
            ),
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    plan = sub.add_parser("plan")
    plan.add_argument("--assemblies-root", type=Path, required=True)
    plan.add_argument("--output", type=Path, required=True)
    plan.add_argument("--window", type=int, default=8192)
    plan.add_argument("--minimum-window", type=int, default=1024)
    plan.add_argument("--aggregation", choices=AGGREGATIONS, default="mean")
    plan.add_argument("--single-strand", action="store_true")
    plan.add_argument("--model-id", required=True)
    plan.add_argument("--model-revision", required=True)
    plan.add_argument("--layer", required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = write_plan(
        args.assemblies_root,
        args.output,
        window=args.window,
        minimum_window=args.minimum_window,
        aggregation=args.aggregation,
        both_strands=not args.single_strand,
        model_id=args.model_id,
        model_revision=args.model_revision,
        layer=args.layer,
    )
    print(json.dumps({k: v for k, v in payload.items() if k != "assemblies"}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
