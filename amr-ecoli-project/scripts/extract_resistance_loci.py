"""Extract resistance-relevant sequence units per genome for Evo 2 embedding.

docs/FOUNDATION_MODEL_PROTOCOL.md revision 0.2 (adopted 2026-09-27). Rules,
fixed before any embedding exists:

- Chromosomal panel: PANEL genes of E. coli K-12 MG1655, each with 300 bp
  upstream (promoters/attenuators), located in each assembly by blastn. The
  best hit counts when identity >= 80% and it covers >= 90% of the reference
  unit; the unit is extended to the reference's full extent, oriented to the
  reference strand. Otherwise the unit is recorded missing, never imputed.
- Acquired determinants: every AMRFinderPlus Type AMR / Subtype AMR call that
  is complete (same rule as build_amr_features.py), with 500 bp flanks,
  oriented to the gene's strand.

Genotype only; no phenotype is read. Subcommands:
    reference  build the panel FASTA from the reference genome + GFF3
    extract    extract units for a batch of genomes (array task)
    collect    deduplicate sequences across genomes for embedding
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

try:
    from scripts.build_amr_features import is_present, parse_amrfinder
except ModuleNotFoundError:  # Direct execution
    from build_amr_features import is_present, parse_amrfinder  # type: ignore[no-redef]

PANEL = ("gyrA", "gyrB", "parC", "parE", "marR", "acrR", "soxR", "ampC", "ompF", "ompC", "ompR", "envZ")
UPSTREAM = 300
FLANK = 500
MIN_IDENTITY = 80.0
MIN_COVERAGE = 0.90
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class LocusError(ValueError):
    """Raised when inputs do not support the extraction contract."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = None
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in records:
                    raise LocusError(f"Duplicate FASTA record {current} in {path}")
                records[current] = []
            elif current is None:
                raise LocusError(f"Sequence before header in {path}")
            else:
                records[current].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


# ---------------------------------------------------------------- reference

def panel_from_reference(genome: dict[str, str], gff: Path) -> list[tuple[str, str]]:
    """Each PANEL gene with UPSTREAM bp upstream, on the gene's own strand."""
    found: dict[str, tuple[str, int, int, str]] = {}
    with gff.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            name = attributes.get("gene") or attributes.get("Name")
            if name in PANEL:
                if name in found:
                    raise LocusError(f"Reference annotates {name} more than once")
                found[name] = (fields[0], int(fields[3]), int(fields[4]), fields[6])
    missing = [name for name in PANEL if name not in found]
    if missing:
        raise LocusError(f"Reference lacks panel genes: {missing}")
    units = []
    for name in PANEL:
        contig, start, end, strand = found[name]
        sequence = genome[contig]
        if strand == "+":
            units.append((name, sequence[max(0, start - 1 - UPSTREAM):end]))
        else:
            units.append((name, reverse_complement(sequence[start - 1:min(len(sequence), end + UPSTREAM)])))
    return units


# ------------------------------------------------------------------ extract

BLAST_FIELDS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen"]


def run_blastn(blastn: str, query: Path, subject: Path) -> str:
    result = subprocess.run(
        [blastn, "-query", str(query), "-subject", str(subject), "-outfmt", "6 " + " ".join(BLAST_FIELDS),
         "-evalue", "1e-20", "-max_hsps", "1", "-max_target_seqs", "5"],
        check=True, capture_output=True, text=True)
    return result.stdout


def best_panel_hits(blast_tsv: str) -> dict[str, dict[str, Any]]:
    best: dict[str, dict[str, Any]] = {}
    for line in blast_tsv.splitlines():
        if not line.strip():
            continue
        row = dict(zip(BLAST_FIELDS, line.split("\t")))
        hit = {"contig": row["sseqid"], "pident": float(row["pident"]), "qstart": int(row["qstart"]),
               "qend": int(row["qend"]), "sstart": int(row["sstart"]), "send": int(row["send"]),
               "bitscore": float(row["bitscore"]), "qlen": int(row["qlen"])}
        current = best.get(row["qseqid"])
        if current is None or (hit["bitscore"], hit["contig"]) > (current["bitscore"], current["contig"]):
            best[row["qseqid"]] = hit
    return best


def panel_unit(hit: dict[str, Any], contigs: dict[str, str]) -> tuple[str | None, str]:
    """(sequence or None, status) for one panel locus hit."""
    coverage = (hit["qend"] - hit["qstart"] + 1) / hit["qlen"]
    if hit["pident"] < MIN_IDENTITY or coverage < MIN_COVERAGE:
        return None, f"below_threshold(identity={hit['pident']:.1f},coverage={coverage:.3f})"
    contig = contigs[hit["contig"]]
    head, tail = hit["qstart"] - 1, hit["qlen"] - hit["qend"]
    if hit["sstart"] <= hit["send"]:
        start, end = hit["sstart"] - head, hit["send"] + tail
        if start < 1 or end > len(contig):
            return None, "truncated_at_contig_end"
        return contig[start - 1:end], "found"
    start, end = hit["send"] - tail, hit["sstart"] + head
    if start < 1 or end > len(contig):
        return None, "truncated_at_contig_end"
    return reverse_complement(contig[start - 1:end]), "found"


def acquired_units(amrfinder: Path, contigs: dict[str, str]) -> list[tuple[str, str, str]]:
    """(unit name, sequence, status) for complete acquired AMR gene calls."""
    units = []
    seen: dict[str, int] = defaultdict(int)
    for hit in parse_amrfinder(amrfinder):
        if hit["Type"].strip() != "AMR" or hit["Subtype"].strip() != "AMR" or not is_present(hit["Method"]):
            continue
        symbol = hit["Element symbol"].strip()
        seen[symbol] += 1
        name = f"{symbol}#{seen[symbol]}" if seen[symbol] > 1 else symbol
        contig_id = hit["Contig id"].strip()
        if contig_id not in contigs:
            raise LocusError(f"AMRFinderPlus contig {contig_id} not in assembly")
        contig = contigs[contig_id]
        start, stop = int(hit["Start"]), int(hit["Stop"])
        sequence = contig[max(0, start - 1 - FLANK):min(len(contig), stop + FLANK)]
        if hit["Strand"].strip() == "-":
            sequence = reverse_complement(sequence)
        units.append((name, sequence, "found"))
    return units


def extract_genome(isolate: str, cohort: str, assembly: Path, amrfinder: Path, panel: Path,
                   blast: Callable[[Path, Path], str]) -> list[dict[str, str]]:
    contigs = read_fasta(assembly)
    hits = best_panel_hits(blast(panel, assembly))
    records = []
    for name in PANEL:
        if name not in hits:
            sequence, status = None, "not_found"
        else:
            sequence, status = panel_unit(hits[name], contigs)
        records.append({"isolate_id": isolate, "cohort": cohort, "unit": f"panel:{name}", "status": status,
                        "sequence": sequence or ""})
    for name, sequence, status in acquired_units(amrfinder, contigs):
        records.append({"isolate_id": isolate, "cohort": cohort, "unit": f"acquired:{name}", "status": status,
                        "sequence": sequence})
    for record in records:
        record["length"] = str(len(record["sequence"]))
        record["sequence_sha256"] = sha256_text(record["sequence"]) if record["sequence"] else ""
    return records


# ------------------------------------------------------------------ collect

def collect(unit_tables: list[Path]) -> tuple[dict[str, str], list[dict[str, str]]]:
    unique: dict[str, str] = {}
    mapping: list[dict[str, str]] = []
    for table in unit_tables:
        with table.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["sequence"]:
                    digest = sha256_text(row["sequence"])
                    if digest != row["sequence_sha256"]:
                        raise LocusError(f"Sequence hash mismatch in {table} for {row['isolate_id']} {row['unit']}")
                    unique.setdefault(digest, row["sequence"])
                mapping.append({key: row[key] for key in
                                ("isolate_id", "cohort", "unit", "status", "length", "sequence_sha256")})
    return unique, mapping


UNIT_FIELDS = ["isolate_id", "cohort", "unit", "status", "length", "sequence_sha256", "sequence"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    ref = sub.add_parser("reference")
    ref.add_argument("--genome", type=Path, required=True)
    ref.add_argument("--gff", type=Path, required=True)
    ref.add_argument("--output", type=Path, required=True)
    ext = sub.add_parser("extract")
    ext.add_argument("--genome-manifest", type=Path, required=True)
    ext.add_argument("--panel", type=Path, required=True)
    ext.add_argument("--blastn", required=True)
    ext.add_argument("--assemblies-root", type=Path, default=Path("results/genomics/disjointness/assemblies"))
    ext.add_argument("--analysis-root", type=Path, default=Path("results/genomics/analysis"))
    ext.add_argument("--start", type=int, required=True)
    ext.add_argument("--count", type=int, required=True)
    ext.add_argument("--output", type=Path, required=True)
    col = sub.add_parser("collect")
    col.add_argument("--units", type=Path, nargs="+", required=True)
    col.add_argument("--unique-fasta", type=Path, required=True)
    col.add_argument("--mapping", type=Path, required=True)
    col.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()

    try:
        if args.command == "reference":
            units = panel_from_reference(read_fasta(args.genome), args.gff)
            write_fasta(args.output, units)
            print(json.dumps({name: len(seq) for name, seq in units}))
        elif args.command == "extract":
            with args.genome_manifest.open("r", encoding="utf-8", newline="") as handle:
                genomes = [r for r in csv.DictReader(handle, delimiter="\t") if r["eligible_for_modeling"] == "true"]
            batch = genomes[args.start:args.start + args.count]
            rows = []
            for genome in batch:
                cohort, isolate = genome["cohort"], genome["isolate_id"]
                rows.extend(extract_genome(
                    isolate, cohort, args.assemblies_root / cohort / isolate / "assembly.fna",
                    args.analysis_root / cohort / isolate / "amrfinder" / "amrfinder.tsv", args.panel,
                    lambda query, subject: run_blastn(args.blastn, query, subject)))
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=UNIT_FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
        else:
            unique, mapping = collect(args.units)
            write_fasta(args.unique_fasta, sorted(unique.items()))
            with args.mapping.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=UNIT_FIELDS[:-1], delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(mapping)
            statuses: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
            for row in mapping:
                statuses[row["unit"].split(":")[0] if row["unit"].startswith("acquired") else row["unit"]][row["status"].split("(")[0]] += 1
            total = sum(len(s) for s in unique.values())
            payload = {"schema_version": "1.0.0", "operation": "collect_resistance_loci",
                       "generated_at_utc": utc_now(), "genomes": len({r["isolate_id"] for r in mapping}),
                       "unit_rows": len(mapping), "unique_sequences": len(unique), "unique_bases": total,
                       "status_counts": {k: dict(v) for k, v in sorted(statuses.items())},
                       "rules": {"panel": list(PANEL), "upstream": UPSTREAM, "flank": FLANK,
                                 "min_identity": MIN_IDENTITY, "min_coverage": MIN_COVERAGE}}
            args.manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            print(json.dumps({k: payload[k] for k in ("genomes", "unit_rows", "unique_sequences", "unique_bases")}))
    except LocusError as exc:
        raise SystemExit(f"Locus extraction failed: {exc}") from exc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
