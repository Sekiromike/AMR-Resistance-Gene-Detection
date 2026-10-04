"""Nearest-relative baseline (docs/EVALUATION_DESIGN_AMENDMENT.md, addition 4).

For each outer-test isolate and drug, predict from the training isolate
(outer training folds, same population, with that drug's MIC) of highest ANI
in the frozen development-versus-development comparisons (ties: higher
minimum aligned fraction, then isolate ID). The prediction is that relative's
reference-set boundary: its exact dilution, or the finite end of its censored
range. No tuning. A model that cannot beat this has learned relatedness, not
resistance.

The comparison table (~9 GB) is streamed once; for every isolate the best
TOP_K relatives in other folds are kept, separately for all sources and for
human-clinical relatives, using the ANI with the isolate as query. Development
only.

Run from the repository root:
    python scripts/run_neighbor_baseline.py --cohort ... --splits ... \\
        --comparisons ... --output-prefix ... --manifest ...
"""
from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

try:
    from scripts.evaluate_mic_predictions import reference_bounds
except ModuleNotFoundError:  # Direct execution
    from evaluate_mic_predictions import reference_bounds  # type: ignore[no-redef]

TOP_K = 200
POPULATIONS = ("human_clinical", "all_sources")
PREFIX = "development::"


class NeighborError(ValueError):
    """Raised when the baseline cannot be computed as specified."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def boundary_prediction(sign: str, value: float) -> float:
    lo, hi = reference_bounds(sign, value)
    return hi if math.isinf(lo) else lo


def collect_neighbors(comparisons: Path, fold: dict[str, int], clinical: set[str],
                      top_k: int = TOP_K) -> dict[str, dict[str, list[tuple[float, float, str]]]]:
    """Best other-fold relatives per isolate: {population: {isolate: [(ani, af, id), ...]}}."""
    heaps: dict[str, dict[str, list]] = {pop: {} for pop in POPULATIONS}

    def offer(query: str, relative: str, key: tuple[float, float, str]) -> None:
        for pop in POPULATIONS:
            if pop == "human_clinical" and relative not in clinical:
                continue
            heap = heaps[pop].setdefault(query, [])
            item = (key[0], key[1], _reverse(relative))  # min-heap keeps the best top_k
            if len(heap) < top_k:
                heapq.heappush(heap, item)
            elif item > heap[0]:
                heapq.heapreplace(heap, item)

    with comparisons.open("r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        expected = ["development_sequence_id", "external_sequence_id", "ani_percent",
                    "aligned_fraction_development", "aligned_fraction_external"]
        if header[:5] != expected:
            raise NeighborError(f"Unexpected comparison header: {header}")
        # The table lists (almost) every ordered pair. A row offers its second
        # genome as a relative of its first only, so each ordered pair counts
        # once without holding a set of ~53M pairs in memory.
        for line in handle:
            a, b, ani, af_a, af_b = line.rstrip("\n").split("\t")[:5]
            a, b = a.removeprefix(PREFIX), b.removeprefix(PREFIX)
            if a == b or a not in fold or b not in fold or fold[a] == fold[b]:
                continue
            offer(a, b, (float(ani), min(float(af_a), float(af_b)), b))
    return {pop: {q: sorted(((ani, af, _reverse(r)) for ani, af, r in heap), key=lambda t: (-t[0], -t[1], t[2]))
                  for q, heap in per.items()} for pop, per in heaps.items()}


def _reverse(isolate: str) -> str:
    """Order key so that, at equal ANI and AF, the lexicographically smaller ID ranks higher."""
    return "".join(chr(0x10FFFF - ord(c)) for c in isolate)


def predict(cohort: pd.DataFrame, fold: dict[str, int], neighbors: dict[str, dict[str, list]],
            population: str) -> pd.DataFrame:
    development = cohort.loc[cohort["evaluation_split"].astype(str).eq("development")].copy()
    if population == "human_clinical":
        development = development.loc[development["intended_use_population"].astype(str).eq("human_clinical")]
    development["cv_fold"] = development["isolate_id"].map(fold)
    reference = {(r.isolate_id, r.antibiotic): (str(r.measurement_sign), float(r.ast_value))
                 for r in development.itertuples()}
    rows, missing = [], 0
    for row in development.itertuples():
        chosen = None
        for ani, af, relative in neighbors[population].get(row.isolate_id, []):
            if (relative, row.antibiotic) in reference:
                chosen = (relative, ani)
                break
        if chosen is None:
            missing += 1
            continue
        sign, value = reference[(chosen[0], row.antibiotic)]
        rows.append({
            "isolate_id": row.isolate_id, "antibiotic": row.antibiotic, "evaluation_split": "development",
            "lineage_group": row.lineage_group, "intended_use_population": row.intended_use_population,
            "cv_fold": int(row.cv_fold), "measurement_sign": row.measurement_sign, "ast_value": row.ast_value,
            "predicted_log2_mic": boundary_prediction(sign, value),
            "neighbor_isolate_id": chosen[0], "neighbor_ani": chosen[1],
        })
    if missing:
        raise NeighborError(f"{missing} {population} rows have no training relative with that drug in the top {TOP_K}")
    return pd.DataFrame(rows).sort_values(["antibiotic", "isolate_id"]).reset_index(drop=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--comparisons", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    cohort = pd.read_csv(args.cohort, dtype={"isolate_id": str, "lineage_group": str}, low_memory=False)
    splits = pd.read_csv(args.splits, dtype={"isolate_id": str})
    development = splits[splits["evaluation_split"].astype(str).eq("development")]
    fold = dict(zip(development["isolate_id"], development["cv_fold"].astype(int)))
    in_cohort = set(cohort.loc[cohort["evaluation_split"].eq("development"), "isolate_id"])
    fold = {isolate: f for isolate, f in fold.items() if isolate in in_cohort}
    clinical = set(cohort.loc[cohort["intended_use_population"].eq("human_clinical"), "isolate_id"])
    try:
        neighbors = collect_neighbors(args.comparisons, fold, clinical)
        outputs = {}
        for population in POPULATIONS:
            frame = predict(cohort, fold, neighbors, population)
            path = Path(f"{args.output_prefix}_{population}.csv")
            path.parent.mkdir(parents=True, exist_ok=True)
            frame.to_csv(path, index=False)
            outputs[population] = {"path": str(path), "sha256": sha256_file(path), "rows": int(len(frame)),
                                   "median_neighbor_ani": float(frame["neighbor_ani"].median())}
    except NeighborError as exc:
        raise SystemExit(f"Nearest-relative baseline failed: {exc}") from exc
    payload = {
        "schema_version": "1.0.0",
        "operation": "run_neighbor_baseline",
        "model": "nearest other-fold relative by ANI; prediction is its reference-set boundary",
        "generated_at_utc": utc_now(),
        "top_k": TOP_K,
        "inputs": {name: {"path": str(path), "sha256": sha256_file(path)} for name, path in
                   (("cohort", args.cohort), ("splits", args.splits))} | {
                       "comparisons": {"path": str(args.comparisons)}},
        "outputs": outputs,
        "scientific_boundary": {"external_rows_read": False, "tuning": "none"},
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
