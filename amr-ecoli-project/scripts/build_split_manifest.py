"""Create deterministic, genomic-cluster-disjoint development folds.

External membership must already be assigned from study provenance. This script
never chooses or modifies the external cohort; it only assigns development
clusters to cross-validation folds while balancing drug/stratum counts.

Balancing stratum (development rows only): the S/I/R category under the
categorical endpoint; under the MIC endpoint, whether the reference MIC is
exact, left-censored (at or below the panel) or right-censored (above it),
combined with the intended-use population when present, so the primary
human-clinical subset is balanced as well. External rows are reported by
antibiotic count only: their censoring mix approximates resistance prevalence
and the external set is sealed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


CENSORING_CLASS = {"=": "exact", "==": "exact", "<": "left", "<=": "left", ">": "right", ">=": "right"}


def balance_stratum(df: pd.DataFrame) -> pd.Series:
    """Per-row balancing stratum; see module docstring."""
    if "ast_category" in df.columns and df["ast_category"].astype(str).str.strip().ne("").all():
        stratum = df["ast_category"].astype(str)
    elif "measurement_sign" in df.columns:
        signs = df["measurement_sign"].astype(str).str.strip()
        unknown = sorted(set(signs) - set(CENSORING_CLASS))
        if unknown:
            raise ValueError(f"Unknown MIC comparators: {unknown}")
        stratum = signs.map(CENSORING_CLASS)
    else:
        raise ValueError("Split balancing needs ast_category (categorical) or measurement_sign (MIC).")
    if "intended_use_population" in df.columns:
        stratum = df["intended_use_population"].astype(str) + "|" + stratum
    return stratum


def development_cv_groups(df: pd.DataFrame, n_folds: int = 5) -> dict[str, str]:
    """Map each isolate to a cross-validation group.

    Amendment 007: the unit is the union of `lineage_group` (SLV group,
    amendment 004), `genomic_cluster` (99.9% ANI component) and
    `deduplication_group`, so neither lineages, near-neighbours nor duplicates
    straddle a fold boundary. In development the clusters nest inside lineages
    (489 lineage groups become 485 units). Should clusters ever chain lineages
    together, as they did in the external set, a unit larger than one fold's
    share stops the split instead of silently collapsing the folds.
    """
    required = ("isolate_id", "lineage_group", "genomic_cluster", "deduplication_group", "evaluation_split")
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(
            f"Cross-validation grouping requires columns absent from the cohort: {missing}"
        )
    development = df.loc[df["evaluation_split"].eq("development")]
    isolates = sorted(development["isolate_id"].astype(str).unique())
    parent = {isolate: isolate for isolate in isolates}

    def find(isolate: str) -> str:
        while parent[isolate] != isolate:
            parent[isolate] = parent[parent[isolate]]
            isolate = parent[isolate]
        return isolate

    def union(left: str, right: str) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    # Population structure: one fold per lineage, never split a lineage.
    for _, rows in development.groupby("lineage_group", sort=False):
        members = sorted(rows["isolate_id"].astype(str).unique())
        for member in members[1:]:
            union(members[0], member)

    def largest_unit() -> int:
        return max(Counter(find(isolate) for isolate in isolates).values(), default=0)

    # Leakage: duplicates stay together even when assigned different lineages.
    def unite_by(column: str) -> None:
        for _, rows in development.groupby(column, sort=False):
            members = sorted(rows["isolate_id"].astype(str).unique())
            for member in members[1:]:
                union(members[0], member)

    unite_by("deduplication_group")
    before_clusters = largest_unit()
    # Amendment 007: 99.9% ANI near-neighbours stay together too. If that
    # merges units into one larger than both the largest lineage/duplicate
    # unit and one fold's share, clusters have chained: stop, do not collapse.
    unite_by("genomic_cluster")
    after_clusters = largest_unit()
    if after_clusters > before_clusters and after_clusters > len(isolates) / n_folds:
        raise ValueError(
            f"Genomic clusters chained lineages into a unit of {after_clusters} of {len(isolates)} development "
            f"isolates (largest lineage/duplicate unit: {before_clusters}), more than one fold's share "
            "(amendment 007)."
        )

    components: dict[str, list[str]] = {}
    for isolate in isolates:
        components.setdefault(find(isolate), []).append(isolate)
    mapping = {}
    for index, root in enumerate(sorted(components), start=1):
        label = f"cv-group-{index:06d}"
        for isolate in components[root]:
            mapping[isolate] = label
    return mapping


def assign_development_folds(df: pd.DataFrame, n_folds: int = 5, seed: int = 20260815) -> dict[str, int]:
    development = df.loc[df["evaluation_split"].eq("development")].copy()
    if development.empty:
        raise ValueError("No development rows are available for fold assignment.")
    isolate_to_group = development_cv_groups(df, n_folds)
    development["cv_group"] = development["isolate_id"].astype(str).map(isolate_to_group)
    development["_stratum"] = balance_stratum(development)
    groups = sorted(development["cv_group"].unique())
    if len(groups) < n_folds:
        raise ValueError(f"Need at least {n_folds} development cross-validation groups; observed {len(groups)}.")

    dimensions = sorted(
        set(zip(development["antibiotic"].astype(str), development["_stratum"].astype(str)))
    )
    dim_index = {dimension: index for index, dimension in enumerate(dimensions)}
    group_vectors: dict[str, np.ndarray] = {}
    for group, rows in development.groupby("cv_group", sort=False):
        vector = np.zeros(len(dimensions) + 1, dtype=float)
        for key, count in rows.groupby(["antibiotic", "_stratum"]).size().items():
            vector[dim_index[(str(key[0]), str(key[1]))]] = float(count)
        vector[-1] = float(len(rows))
        group_vectors[str(group)] = vector

    total = sum(group_vectors.values(), np.zeros(len(dimensions) + 1, dtype=float))
    target = total / n_folds
    scale = np.maximum(target, 1.0)

    def tie_key(group: str) -> str:
        return hashlib.sha256(f"{seed}:{group}".encode("utf-8")).hexdigest()

    ordered = sorted(groups, key=lambda group: (-group_vectors[group][-1], tie_key(group)))
    fold_totals = np.zeros((n_folds, len(dimensions) + 1), dtype=float)
    assignments: dict[str, int] = {}
    for position, group in enumerate(ordered):
        if position < n_folds:
            selected = position
        else:
            scores = []
            for fold in range(n_folds):
                candidate = fold_totals.copy()
                candidate[fold] += group_vectors[group]
                imbalance = np.square((candidate - target) / scale).sum()
                scores.append((float(imbalance), float(candidate[fold, -1]), fold))
            selected = min(scores)[2]
        assignments[group] = selected
        fold_totals[selected] += group_vectors[group]
    return assignments


def build_manifest(df: pd.DataFrame, n_folds: int = 5, seed: int = 20260815) -> pd.DataFrame:
    required = {
        "isolate_id",
        "antibiotic",
        "evaluation_split",
        "genomic_cluster",
        "lineage_group",
        "deduplication_group",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing split columns: {missing}")
    if not set(df["evaluation_split"].astype(str)).issubset({"development", "external"}):
        raise ValueError("evaluation_split must contain only development or external.")

    cluster_splits = df.groupby("genomic_cluster")["evaluation_split"].nunique()
    if (cluster_splits > 1).any():
        bad = sorted(map(str, cluster_splits[cluster_splits > 1].index))
        raise ValueError(f"Genomic clusters cross development/external boundaries: {bad[:10]}")
    isolate_metadata = df.groupby("isolate_id").agg(
        evaluation_splits=("evaluation_split", "nunique"),
        genomic_clusters=("genomic_cluster", "nunique"),
        lineage_groups=("lineage_group", "nunique"),
    )
    inconsistent = isolate_metadata[(isolate_metadata > 1).any(axis=1)]
    if not inconsistent.empty:
        raise ValueError(f"Isolates have inconsistent split/cluster/lineage metadata: {list(inconsistent.index[:10])}")

    isolate_to_group = development_cv_groups(df, n_folds)
    assignments = assign_development_folds(df, n_folds=n_folds, seed=seed)
    rows = (
        df[["isolate_id", "evaluation_split", "genomic_cluster", "lineage_group"]]
        .drop_duplicates()
        .copy()
    )
    rows["cv_group"] = pd.Series(pd.NA, index=rows.index, dtype="string")
    rows["cv_fold"] = pd.Series(pd.NA, index=rows.index, dtype="Int64")
    development = rows["evaluation_split"].eq("development")
    rows.loc[development, "cv_group"] = (
        rows.loc[development, "isolate_id"].astype(str).map(isolate_to_group)
    )
    rows.loc[development, "cv_fold"] = rows.loc[development, "cv_group"].map(assignments)
    return rows.sort_values(["evaluation_split", "cv_fold", "isolate_id"], na_position="last").reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("data/splits/isolate_splits.csv"))
    parser.add_argument("--manifest", type=Path, default=Path("data/splits/split_manifest.json"))
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260815)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cohort = pd.read_csv(args.cohort)
    split_table = build_manifest(cohort, n_folds=args.folds, seed=args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    split_table.to_csv(args.output, index=False)
    output_hash = sha256_file(args.output)
    merged = cohort.merge(split_table, on=["isolate_id", "evaluation_split", "genomic_cluster", "lineage_group"])
    development = merged.loc[merged["evaluation_split"].eq("development")]
    distribution = (
        development.assign(stratum=balance_stratum(development))
        .groupby(["cv_fold", "antibiotic", "stratum"])
        .size()
        .rename("n")
        .reset_index()
        .to_dict(orient="records")
    )
    external_counts = (
        merged.loc[merged["evaluation_split"].eq("external")]
        .groupby("antibiotic")
        .size()
        .to_dict()
    )
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "algorithm": "deterministic greedy multilabel balance of lineage/deduplication groups over antibiotic x stratum",
        "seed": args.seed,
        "n_folds": args.folds,
        "cohort_sha256": sha256_file(args.cohort),
        "split_sha256": output_hash,
        "n_isolates": int(split_table["isolate_id"].nunique()),
        "n_development_clusters": int(
            split_table.loc[split_table["evaluation_split"].eq("development"), "genomic_cluster"].nunique()
        ),
        "n_development_cv_groups": int(
            split_table.loc[split_table["evaluation_split"].eq("development"), "cv_group"].nunique()
        ),
        "n_external_clusters": int(
            split_table.loc[split_table["evaluation_split"].eq("external"), "genomic_cluster"].nunique()
        ),
        "distribution": distribution,
        "external_rows_by_antibiotic": {str(k): int(v) for k, v in external_counts.items()},
        "external_selection": "PRESPECIFIED_INPUT_NOT_MODIFIED",
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {len(split_table)} isolate assignments to {args.output} (sha256={output_hash})")


if __name__ == "__main__":
    main()
