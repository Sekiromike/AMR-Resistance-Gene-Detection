"""Finalize per-genome analyses into the genome manifest used for the cohort.

Reads the outputs of hpc/unity/genome-analysis-array.sbatch, which ran QUAST,
MLST, AMRFinderPlus and a skani species check on the same checksum-verified
assemblies used for genomic disjointness. Keyed by isolate, so genomes
assembled from raw reads -- which have no NCBI assembly accession -- are
represented honestly rather than forced into an accession-keyed schema.

Rules are fixed before any result is inspected:

- Species: "Escherichia coli" when skani ANI to the frozen E. coli type-strain
  reference is >= 95% and >= 50% of the query genome aligns (the GTDB species
  convention). Anything else is not confirmed. ANI cannot separate Shigella,
  which falls inside the E. coli species boundary; this is a stated limitation.
- Assembly QC: the thresholds frozen in config/study.json genomics.assembly_qc.
- Sequence type: "ST<n>" only when MLST reports an integer sequence type with
  status PERFECT or OKAY. Novel, mixed or missing calls leave lineage_group
  blank (amendment 004: never a placeholder).
- Lineage group (amendment 004, revised 2026-09-27): single-locus-variant
  (SLV) groups. STs of eligible development genomes are linked when their
  Achtman profiles differ at exactly one locus; each connected component is
  one group, named "SLV:" plus its most frequent ST (ties: lowest ST number).
  Only development defines the groups. Any other genome joins the group of its
  ST, else the single development group it is an SLV of, else a group of its
  own ST; it never merges development groups. lineage_basis records which.
- MLST scheme: every report must name the scheme and exactly the loci frozen
  in config/study.json genomics.mlst (Achtman). Any other scheme is a pipeline
  defect and stops the run; it is never recorded as an unresolved lineage.
  MLST may be read from a separate, checksum-verified re-typing tree
  (--mlst-root) instead of the analysis directory.

Run from the repository root:
    python scripts/finalize_genome_analyses.py --requests ... --output ...
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    from scripts.build_genome_manifest import (
        _quast_value,
        _read_exit_code,
        count_amrfinder_hits,
        parse_mlst_report,
        parse_quast_report,
    )
except ModuleNotFoundError:  # Direct execution
    from build_genome_manifest import (  # type: ignore[no-redef]
        _quast_value,
        _read_exit_code,
        count_amrfinder_hits,
        parse_mlst_report,
        parse_quast_report,
    )

SPECIES_NAME = "Escherichia coli"
SPECIES_MIN_ANI = 95.0
SPECIES_MIN_QUERY_AF = 50.0  # percent, as skani reports it
LINEAGE_STATUSES = frozenset({"PERFECT", "OKAY"})
ALLELE_LOCUS = re.compile(r"([A-Za-z0-9_]+)\(")

OUTPUT_FIELDS = [
    "isolate_id",
    "cohort",
    "sequence_id",
    "biosample_accession",
    "genome_source",
    "genome_source_accession",
    "assembly_accession",
    "assembly_sha256",
    "species_method",
    "species_result",
    "species_ani",
    "species_query_aligned_fraction",
    "genome_qc_status",
    "quast_total_length",
    "quast_contigs",
    "quast_n50",
    "quast_gc_percent",
    "quast_ns_per_100kb",
    "mlst_sequence_type",
    "mlst_status",
    "mlst_alleles",
    "lineage_group",
    "lineage_basis",
    "amrfinder_hit_count",
    "eligible_for_modeling",
    "overall_status",
    "exclusion_reasons",
    "warning_reasons",
]


class FinalizeError(ValueError):
    """Raised when inputs are missing or inconsistent."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _verify_artifacts(directory: Path) -> bool:
    manifest = directory / "artifact.sha256"
    if not manifest.is_file():
        return False
    for line in manifest.read_text(encoding="utf-8").splitlines():
        parts = line.split(maxsplit=1)
        if len(parts) != 2:
            return False
        relative = parts[1].lstrip("* ").removeprefix("./")
        target = directory / relative
        if not target.is_file() or sha256_file(target) != parts[0].lower():
            return False
    return True


def parse_species(path: Path) -> tuple[float | None, float | None]:
    """Return (ANI, query aligned fraction %) from a skani dist table, or Nones."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        return None, None  # skani emits no row below its ANI floor
    if len(rows) != 1:
        raise FinalizeError(f"Expected one skani row in {path}, found {len(rows)}")
    return float(rows[0]["ANI"]), float(rows[0]["Align_fraction_query"])


def lineage_from_mlst(sequence_type: str, status: str) -> str:
    st = sequence_type.strip()
    if st.isdigit() and int(st) > 0 and status.strip().upper() in LINEAGE_STATUSES:
        return f"ST{int(st)}"
    return ""


def check_mlst_scheme(mlst: dict[str, str], expected: dict[str, Any], where: Path) -> None:
    """Fail closed unless the report used the frozen scheme and loci."""
    loci = ALLELE_LOCUS.findall(mlst.get("alleles", ""))
    if mlst.get("scheme") != expected["scheme"] or sorted(loci) != sorted(expected["loci"]):
        raise FinalizeError(
            f"MLST report {where} used scheme {mlst.get('scheme')!r} with loci {loci}; "
            f"expected {expected['scheme']!r} with loci {expected['loci']}"
        )


def _profile(alleles: str) -> tuple[tuple[str, str], ...]:
    pairs = []
    for item in alleles.split(";"):
        locus, _, value = item.partition("(")
        pairs.append((locus.strip(), value.rstrip(")").strip()))
    return tuple(sorted(pairs))


def _st_number(st: str) -> int:
    return int(st.removeprefix("ST"))


def assign_slv_groups(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Replace per-row ST lineage with single-locus-variant groups, in place."""
    typed = [row for row in rows if row.get("lineage_group")]
    profiles: dict[str, tuple[tuple[str, str], ...]] = {}
    for row in typed:
        profile = _profile(str(row.get("mlst_alleles", "")))
        known = profiles.setdefault(row["lineage_group"], profile)
        if known != profile:
            raise FinalizeError(f"{row['lineage_group']} has inconsistent allele profiles")

    def slv(a: str, b: str) -> bool:
        pa, pb = profiles[a], profiles[b]
        if [locus for locus, _ in pa] != [locus for locus, _ in pb]:
            raise FinalizeError(f"{a} and {b} were typed on different loci")
        return sum(x != y for (_, x), (_, y) in zip(pa, pb)) == 1

    defining = [row for row in typed if row["cohort"] == "development" and row["eligible_for_modeling"] == "true"]
    frequency = Counter(row["lineage_group"] for row in defining)
    development_sts = sorted(frequency, key=_st_number)
    parent = {st: st for st in development_sts}

    def find(st: str) -> str:
        while parent[st] != st:
            parent[st] = parent[parent[st]]
            st = parent[st]
        return st

    for index, a in enumerate(development_sts):
        for b in development_sts[index + 1:]:
            if slv(a, b):
                parent[find(a)] = find(b)
    members: dict[str, list[str]] = {}
    for st in development_sts:
        members.setdefault(find(st), []).append(st)
    group_of: dict[str, str] = {}
    for sts in members.values():
        name = "SLV:" + min(sts, key=lambda st: (-frequency[st], _st_number(st)))
        for st in sts:
            group_of[st] = name

    bases: Counter[str] = Counter()
    for row in typed:
        st = row["lineage_group"]
        if st in group_of:
            row["lineage_group"], row["lineage_basis"] = group_of[st], "development_slv_graph"
        else:
            neighbours = {group_of[other] for other in development_sts if slv(st, other)}
            if len(neighbours) == 1:
                row["lineage_group"], row["lineage_basis"] = neighbours.pop(), "slv_of_development_group"
            else:
                row["lineage_group"], row["lineage_basis"] = "SLV:" + st, "outside_development_groups"
        bases[row["lineage_basis"]] += 1
    return {
        "development_sequence_types": len(development_sts),
        "development_lineage_groups": len(members),
        "assignment_basis_counts": dict(sorted(bases.items())),
    }


def finalize_genome(
    request: dict[str, str],
    analysis_root: Path,
    thresholds: dict[str, float],
    reference_accession: str,
    mlst_config: dict[str, Any],
    mlst_root: Path | None = None,
) -> dict[str, Any]:
    cohort, isolate = request["cohort"], request["isolate_id"]
    source_type = request["source_type"]
    directory = analysis_root / cohort / isolate
    exclusions: list[str] = []
    warnings: list[str] = []
    row: dict[str, Any] = {
        "isolate_id": isolate,
        "cohort": cohort,
        "sequence_id": request["sequence_id"],
        "biosample_accession": request.get("biosample_accession", ""),
        "genome_source": source_type,
        "genome_source_accession": request.get("sequence_accession", ""),
        # Raw-read genomes were assembled here and have no NCBI assembly accession.
        "assembly_accession": request.get("sequence_accession", "") if source_type == "registered_assembly" else "",
        "species_method": (
            f"skani ANI vs E. coli type strain {reference_accession}; "
            f"ANI>={SPECIES_MIN_ANI:g} and query AF>={SPECIES_MIN_QUERY_AF:g}%"
        ),
    }

    if not (directory / "COMPLETE").is_file():
        exclusions.append("ANALYSIS_INCOMPLETE")
    elif not _verify_artifacts(directory):
        exclusions.append("ANALYSIS_CHECKSUM_MISMATCH")
    else:
        provenance = dict(
            line.split("=", 1)
            for line in (directory / "provenance.txt").read_text(encoding="utf-8").splitlines()
            if "=" in line
        )
        row["assembly_sha256"] = provenance.get("assembly_sha256", "")

        # --- assembly QC (thresholds as in build_genome_manifest.finalize_analyses) ---
        quast_exit = _read_exit_code(directory / "quast" / "exit_code.txt")
        metrics: dict[str, float] = {}
        if quast_exit != 0:
            exclusions.append("QUAST_EXECUTION_FAILED")
        elif not (directory / "quast" / "report.tsv").is_file():
            exclusions.append("QUAST_REPORT_MISSING")
        else:
            metrics = parse_quast_report(directory / "quast" / "report.tsv")
        total = _quast_value(metrics, "Total length", "Total length (>= 0 bp)")
        contigs = _quast_value(metrics, "# contigs", "# contigs (>= 0 bp)")
        n50 = _quast_value(metrics, "N50")
        gc = _quast_value(metrics, "GC (%)")
        ns = _quast_value(metrics, "# N's per 100 kbp")
        row.update({
            "quast_total_length": "" if total is None else int(total),
            "quast_contigs": "" if contigs is None else int(contigs),
            "quast_n50": "" if n50 is None else int(n50),
            "quast_gc_percent": "" if gc is None else gc,
            "quast_ns_per_100kb": "" if ns is None else ns,
        })
        qc_fail = [reason for reason in exclusions if reason.startswith("QUAST")]
        if quast_exit == 0 and metrics:
            if None in (total, contigs, n50, ns):
                qc_fail.append("QUAST_REQUIRED_METRICS_MISSING")
            else:
                if total < thresholds["min_total_length"]:
                    qc_fail.append("ASSEMBLY_LENGTH_BELOW_MIN")
                if total > thresholds["max_total_length"]:
                    qc_fail.append("ASSEMBLY_LENGTH_ABOVE_MAX")
                if contigs > thresholds["max_contigs"]:
                    qc_fail.append("ASSEMBLY_CONTIG_COUNT_ABOVE_MAX")
                if n50 < thresholds["min_n50"]:
                    qc_fail.append("ASSEMBLY_N50_BELOW_MIN")
                if ns > thresholds["max_ns_per_100kb"]:
                    qc_fail.append("ASSEMBLY_NS_ABOVE_MAX")
        exclusions.extend(reason for reason in qc_fail if reason not in exclusions)
        row["genome_qc_status"] = "FAIL" if qc_fail else "PASS"

        # --- species ---
        species_exit = _read_exit_code(directory / "species" / "exit_code.txt")
        if species_exit != 0:
            exclusions.append("SPECIES_CHECK_FAILED")
            row["species_result"] = "not_confirmed"
        else:
            ani, af = parse_species(directory / "species" / "skani.tsv")
            row["species_ani"] = "" if ani is None else ani
            row["species_query_aligned_fraction"] = "" if af is None else af
            confirmed = ani is not None and ani >= SPECIES_MIN_ANI and af >= SPECIES_MIN_QUERY_AF
            row["species_result"] = SPECIES_NAME if confirmed else "not_confirmed"
            if not confirmed:
                exclusions.append("SPECIES_NOT_CONFIRMED")

        # --- lineage ---
        mlst_dir = directory / "mlst"
        if mlst_root is not None:
            mlst_dir = mlst_root / cohort / isolate
            if not (mlst_dir / "COMPLETE").is_file() or not _verify_artifacts(mlst_dir):
                raise FinalizeError(f"MLST re-typing output missing or failed checksum: {mlst_dir}")
        mlst_exit = _read_exit_code(mlst_dir / "exit_code.txt")
        if mlst_exit != 0 or not (mlst_dir / "mlst.tsv").is_file():
            warnings.append("MLST_EXECUTION_FAILED")
            row["lineage_group"] = ""
        else:
            try:
                mlst = parse_mlst_report(mlst_dir / "mlst.tsv")
            except ValueError:
                warnings.append("MLST_REPORT_INVALID")
                row["lineage_group"] = ""
            else:
                check_mlst_scheme(mlst, mlst_config, mlst_dir / "mlst.tsv")
                row["mlst_sequence_type"] = mlst["sequence_type"]
                row["mlst_status"] = mlst["status"]
                row["mlst_alleles"] = mlst["alleles"]
                row["lineage_group"] = lineage_from_mlst(mlst["sequence_type"], mlst["status"])
                if not row["lineage_group"]:
                    warnings.append(f"LINEAGE_UNRESOLVED:{mlst['status'] or 'UNKNOWN'}")

        # --- AMR determinants ---
        amr_exit = _read_exit_code(directory / "amrfinder" / "exit_code.txt")
        if amr_exit != 0:
            exclusions.append("AMRFINDER_EXECUTION_FAILED")
        elif not (directory / "amrfinder" / "amrfinder.tsv").is_file():
            exclusions.append("AMRFINDER_REPORT_MISSING")
        else:
            row["amrfinder_hit_count"] = count_amrfinder_hits(directory / "amrfinder" / "amrfinder.tsv")

    row.setdefault("genome_qc_status", "FAIL" if exclusions else "")
    row.setdefault("species_result", "not_confirmed")
    row.setdefault("lineage_group", "")
    row["eligible_for_modeling"] = "true" if not exclusions else "false"
    row["overall_status"] = "PASS" if not exclusions else "EXCLUDED"
    row["exclusion_reasons"] = ";".join(exclusions)
    row["warning_reasons"] = ";".join(warnings)
    return row


def finalize(
    requests_path: Path,
    analysis_root: Path,
    study_config_path: Path,
    output_path: Path,
    manifest_path: Path,
    *,
    reference_accession: str = "GCF_003697165.2",
    mlst_root: Path | None = None,
) -> dict[str, Any]:
    with requests_path.open("r", encoding="utf-8", newline="") as handle:
        requests = list(csv.DictReader(handle, delimiter="\t"))
    if not requests:
        raise FinalizeError("Request table is empty")
    config = json.loads(study_config_path.read_text(encoding="utf-8"))
    qc = config["genomics"]["assembly_qc"]
    thresholds = {key: float(qc[key]) for key in (
        "min_total_length", "max_total_length", "max_contigs", "min_n50", "max_ns_per_100kb")}

    mlst_config = config["genomics"].get("mlst")
    if not mlst_config or not mlst_config.get("scheme") or not mlst_config.get("loci"):
        raise FinalizeError("config/study.json genomics.mlst must name the scheme and its loci")
    rows = [finalize_genome(request, analysis_root, thresholds, reference_accession, mlst_config, mlst_root)
            for request in requests]
    rows.sort(key=lambda row: (row["cohort"], row["isolate_id"]))
    lineage_summary = assign_slv_groups(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    def reasons(field: str) -> dict[str, int]:
        counter: Counter[str] = Counter()
        for row in rows:
            for reason in filter(None, str(row[field]).split(";")):
                counter[reason.split(":")[0]] += 1
        return dict(sorted(counter.items()))

    by_cohort: dict[str, dict[str, int]] = {}
    for row in rows:
        entry = by_cohort.setdefault(row["cohort"], Counter())
        entry["genomes"] += 1
        entry["eligible"] += row["eligible_for_modeling"] == "true"
        entry["lineage_assigned"] += bool(row["lineage_group"])
        entry["species_confirmed"] += row["species_result"] == SPECIES_NAME
        entry["qc_pass"] += row["genome_qc_status"] == "PASS"
    payload = {
        "schema_version": "1.0.0",
        "operation": "finalize_genome_analyses",
        "generated_at_utc": utc_now(),
        "rules": {
            "species": f"ANI>={SPECIES_MIN_ANI} and query AF>={SPECIES_MIN_QUERY_AF}% vs {reference_accession}",
            "sequence_type": "ST<n> when MLST status in PERFECT/OKAY; otherwise blank",
            "lineage": "single-locus-variant groups over eligible development STs (amendment 004 revision)",
            "lineage_summary": lineage_summary,
            "mlst_scheme": mlst_config["scheme"],
            "mlst_loci": mlst_config["loci"],
            "mlst_source": str(mlst_root) if mlst_root is not None else "analysis directory",
            "assembly_qc": qc,
        },
        "inputs": {
            "requests": {"path": str(requests_path), "sha256": sha256_file(requests_path)},
            "study_config": {"path": str(study_config_path), "sha256": sha256_file(study_config_path)},
        },
        "counts": {cohort: dict(values) for cohort, values in sorted(by_cohort.items())},
        "lineages": len({row["lineage_group"] for row in rows if row["lineage_group"]}),
        "exclusion_reason_counts": reasons("exclusion_reasons"),
        "warning_reason_counts": reasons("warning_reasons"),
        "outputs": {"genome_manifest": {"path": str(output_path), "sha256": sha256_file(output_path)}},
        "scientific_boundary": {"phenotypes_read": False},
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--requests", type=Path, required=True)
    parser.add_argument("--analysis-root", type=Path, default=Path("results/genomics/analysis"))
    parser.add_argument("--study-config", type=Path, default=Path("config/study.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--reference-accession", default="GCF_003697165.2")
    parser.add_argument("--mlst-root", type=Path,
                        help="Read MLST from this re-typing tree (<root>/<cohort>/<isolate>/) instead.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        payload = finalize(args.requests, args.analysis_root, args.study_config, args.output,
                           args.manifest, reference_accession=args.reference_accession,
                           mlst_root=args.mlst_root)
    except FinalizeError as exc:
        raise SystemExit(f"Genome finalization failed: {exc}") from exc
    print(json.dumps({k: payload[k] for k in ("counts", "lineages", "exclusion_reason_counts",
                                              "warning_reason_counts")}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
