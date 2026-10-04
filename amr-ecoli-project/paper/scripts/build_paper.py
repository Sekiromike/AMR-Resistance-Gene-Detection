"""Extract published-in-repository summaries and render the manuscript.

Reads the verified result exports created by verify_results.py.
No model fitting or external selection is performed by either script.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper"
DOCS = ROOT / "docs"
STATIC_SOURCES = [
    "docs/EXTERNAL_RESULTS.md", "docs/POSTHOC_FAMILY_ANALYSIS.md",
    "docs/LEAVE_VARIANT_OUT_PLAN.md", "docs/MIC_EVALUATION_SPEC.md",
    "docs/EXTERNAL_EVALUATION_PLAN.md", "docs/FOUNDATION_MODEL_PROTOCOL.md",
    "docs/ENDPOINT_AMENDMENT.md", "docs/EXTERNAL_THRESHOLD_AMENDMENT.md",
    "docs/DEDUPLICATION_POLICY.md", "docs/GROUPING_VARIABLE_AMENDMENT.md",
    "docs/EVALUATION_DESIGN_AMENDMENT.md", "docs/COUNTRY_AND_POPULATION_AMENDMENT.md",
    "docs/JARBS_EXTERNAL_COHORT_AUDIT.md", "docs/GENOMIC_DISJOINTNESS_PROTOCOL.md",
    "config/study.json", "config/external_prediction_lock.tsv",
    "config/external_posthoc_lock.tsv", "scripts/run_mic_models.py",
    "scripts/run_censored_regression.py", "scripts/evaluate_mic_predictions.py",
    "scripts/embedding_features.py", "scripts/embed_resistance_units.py",
    "scripts/extract_resistance_loci.py", "scripts/build_hierarchical_features.py",
    "scripts/leave_variant_out.py", "scripts/predict_external.py",
    "scripts/finalize_genome_analyses.py", "hpc/unity/source-robustness.sbatch",
    "hpc/unity/disjointness-raw-array.sbatch",
]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def clean(text: str) -> str:
    return text.replace("**", "").replace("−", "-").strip()


def numbers(text: str) -> list[float]:
    return [float(x) for x in re.findall(r"[+-]?\d+\.\d+", clean(text))]


def rows(name: str) -> list[tuple[int, list[str]]]:
    result = []
    for line_no, line in enumerate((DOCS / name).read_text().splitlines(), 1):
        if line.startswith("|"):
            result.append((line_no, [clean(c) for c in line.strip("|").split("|")]))
    return result


def evidence_snapshot() -> dict:
    path = PAPER / "provenance/evidence_snapshot.json"
    if path.exists():
        return json.loads(path.read_text())
    state_path = DOCS / "PROJECT_STATUS.json"
    state = json.loads(state_path.read_text())
    records = []
    for phase in state["phases"]:
        if phase["id"] in {"cohort", "frontier_models"}:
            for entry in phase.get("evidence", []):
                if "2026-09-27" <= entry["at_utc"] < "2026-09-30T02":
                    records.append({"phase": phase["id"], **entry})
    out = {"source": "docs/PROJECT_STATUS.json", "source_sha256": sha(state_path),
           "source_updated_at_utc": state["updated_at_utc"], "records": records}
    path.write_text(json.dumps(out, indent=2) + "\n")
    return out


def extract() -> dict[str, list[dict]]:
    from audited_tables import extract as extract_verified
    return extract_verified()


def validate(tables: dict) -> None:
    for name, data in tables.items():
        for row in data:
            for metric in ("share", "common_share", "ea", "undercall"):
                if metric in row:
                    assert 0 <= row[metric] <= 1, (name, row)
            if row.get("ci_low", "") != "":
                point = row.get("ea", row.get("difference"))
                assert row["ci_low"] <= point <= row["ci_high"], (name, row)
    primary = [r for r in tables["leave_variant_out"] if r["features"] == "family"
               and r["learner"] == "ridge" and r["target"].startswith("blaCTX")]
    assert sum(r["n_above"] for r in primary) == 568
    assert round(min(r["share"] for r in primary), 3) == 0.972
    assert all(r["ci_low"] > 0 for r in tables["leave_variant_out_hypotheses"])
    manuscript = (PAPER / "manuscript.md").read_text()
    bib = (PAPER / "references.bib").read_text()
    cites = set(re.findall(r"(?<![\w.])@([A-Za-z][\w-]*)", manuscript))
    entries = set(re.findall(r"@\w+\{([^,]+),", bib))
    assert cites == entries, ("citation mismatch", cites - entries, entries - cites)
    assert "target–isolate observations" in manuscript
    assert "post hoc" in manuscript and "not models retrained" in manuscript
    # Check that the main external table reproduces all six source rows.
    for r in tables["locked_external"]:
        assert f'{r["ea"]:.3f} [{r["ci_low"]:.3f}, {r["ci_high"]:.3f}]' in manuscript
    # Validate prediction-lock structure without pretending predictions are local.
    for file, count in (("external_prediction_lock.tsv", 16), ("external_posthoc_lock.tsv", 2)):
        with (ROOT / "config" / file).open() as handle:
            lock = list(csv.DictReader(handle, delimiter="\t"))
        assert len(lock) == count and len({r["file"] for r in lock}) == count
        assert all(re.fullmatch(r"[a-f0-9]{64}", r["sha256"]) for r in lock)


def write_tables(tables: dict) -> None:
    for name, data in tables.items():
        with (PAPER / "tables" / f"{name}.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    (PAPER / "tables/reported_results.json").write_text(json.dumps(tables, indent=2) + "\n")


def provenance() -> None:
    manifest = {
        "basis": "Original completed-run artifacts; scores and recorded bootstrap intervals reproduced from frozen predictions",
        "repository_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "sources": {p: sha(ROOT / p) for p in STATIC_SOURCES},
        "evidence_snapshot_sha256": sha(PAPER / "provenance/evidence_snapshot.json"),
        "audit_sha256": sha(PAPER / "verification/audit.json"),
    }
    (PAPER / "provenance/source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def figures(tables: dict) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/amr-paper-matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.ticker import PercentFormatter

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 220})
    colors = {"ridge": "#0072B2", "xgboost": "#D55E00"}

    def save(fig, name):
        for ext in ("png", "pdf", "svg"):
            fig.savefig(PAPER / "figures" / f"{name}.{ext}", bbox_inches="tight", facecolor="white")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 3.5), layout="constrained")
    stages = ["Development", "Development excluding source", "Locked external"]
    for learner, offset, marker in (("ridge", -.11, "o"), ("xgboost", .11, "s")):
        data = [next(r for r in tables["gentamicin_gain"] if r["stage"] == s and r["learner"] == learner) for s in stages]
        point = np.array([r["difference"] for r in data]) * 100
        low = np.array([r["ci_low"] for r in data]) * 100
        high = np.array([r["ci_high"] for r in data]) * 100
        ax.errorbar(point, np.arange(3) + offset, xerr=[point-low, high-point], fmt=marker,
                    color=colors[learner], capsize=3, label="Ridge" if learner == "ridge" else "XGBoost")
    ax.set_yticks(range(3), ["Development", "Development, excluding\nPRJNA1297298 evaluation rows", "Locked external"])
    ax.invert_yaxis()
    ax.axvline(0, color="#555555", linestyle="--", linewidth=1)
    ax.set_xlabel("Change in essential agreement (percentage points)")
    ax.set_xlim(-6.5, 5)
    ax.set_ylim(2.55, -.55)
    ax.grid(axis="x", alpha=.2)
    ax.legend(loc="lower right", frameon=False)
    save(fig, "figure1_gentamicin_gain")

    fig, ax = plt.subplots(figsize=(7.2, 3.4), layout="constrained")
    features = ["allele", "amr+evo2", "evo2", "family"]
    for learner, offset, marker in (("ridge", -.10, "o"), ("xgboost", .10, "s")):
        data = [next(r for r in tables["posthoc_external_variants"] if r["features"] == f and r["learner"] == learner) for f in features]
        x = [r["share"] for r in data]
        ax.scatter(x, np.arange(4)+offset, marker=marker, s=55, color=colors[learner], label="Ridge" if learner == "ridge" else "XGBoost")
    ax.set_yticks(range(4), ["Allele features", "Allele + Evo 2", "Evo 2", "Allele + family + drug class"])
    ax.invert_yaxis()
    ax.set_xlim(0, 1.04)
    ax.set_ylim(3.6, -.7)
    ax.xaxis.set_major_formatter(PercentFormatter(1))
    ax.set_xlabel("Share of 102 isolates with a prediction >2 mg/L")
    ax.grid(axis="x", alpha=.2)
    ax.legend(loc="upper center", bbox_to_anchor=(.5, 1.19), ncol=2, frameon=False)
    save(fig, "figure2_external_variants")

    targets = ["blaCTX-M-15", "blaCTX-M-27", "blaCTX-M-14", "blaCTX-M-55", "blaCMY-2"]
    features = ["allele", "amr+evo2", "evo2", "family"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True, layout="constrained")
    for ax, learner in zip(axes, ("ridge", "xgboost")):
        matrix = np.array([[next(r["share"] for r in tables["leave_variant_out"]
                                if r["target"] == t and r["features"] == f and r["learner"] == learner)
                            for f in features] for t in targets])
        im = ax.imshow(matrix, vmin=0, vmax=1, cmap="cividis", aspect="auto")
        for i in range(5):
            for j in range(4):
                ax.text(j, i, f"{matrix[i,j]:.3f}", ha="center", va="center",
                        color="white" if matrix[i,j] < .50 else "#111111", fontsize=10)
        ax.set_xticks(range(4), ["Allele", "Allele\n+ Evo 2", "Evo 2", "Allele + family\n+ drug class"], fontsize=9)
        ax.set_yticks(range(5), ["CTX-M-15 (361)", "CTX-M-27 (134)", "CTX-M-14 (37)", "CTX-M-55 (36)", "CMY-2 control (22)"])
        ax.set_title("Ridge" if learner == "ridge" else "XGBoost", fontweight="bold")
        ax.axhline(3.5, color="white", linewidth=2)
        ax.tick_params(length=0)
    fig.colorbar(im, ax=axes, shrink=.82, label="Share with a prediction >2 mg/L", pad=.025)
    save(fig, "figure3_leave_variant_out")


def render(tectonic: str | None = None) -> None:
    import pypandoc
    from weasyprint import HTML
    os.chdir(PAPER)
    common = ["--standalone", "--citeproc", "--bibliography=references.bib"]
    pypandoc.convert_file("manuscript.md", "html5", outputfile="manuscript.html",
                         extra_args=common + ["--css=style.css", "--mathml"])
    pypandoc.convert_file("manuscript.md", "latex", outputfile="manuscript.tex", extra_args=common)
    tex_path = PAPER / "manuscript.tex"
    tex = tex_path.read_text()
    for name in ("figure1_gentamicin_gain", "figure2_external_variants", "figure3_leave_variant_out"):
        tex = tex.replace(name + ".png", name + ".pdf")
    tex = tex.replace("≥", r"\ensuremath{\geq}").replace("≤", r"\ensuremath{\leq}").replace("−", r"\ensuremath{-}")
    tex = re.sub(r"(\\caption\{)Figure [123]\. ", r"\1", tex)
    tex_path.write_text(tex)
    document = HTML(filename=str(PAPER / "manuscript.html")).render()
    if not tectonic:
        document.write_pdf(PAPER / "manuscript.pdf")
    else:
        (PAPER / "tex-build").mkdir(exist_ok=True)
        command = [tectonic, "--keep-logs", "--outdir", "tex-build", "manuscript.tex"]
        completed = subprocess.run(command, check=True, text=True, capture_output=True)
        (PAPER / "provenance/latex_build.txt").write_text(completed.stdout + completed.stderr)
        import shutil
        shutil.copyfile(PAPER / "tex-build/manuscript.pdf", PAPER / "manuscript.pdf")
    # A PNG contact sheet of page-size layout previews can be rendered separately
    # when a PDF rasterizer is present; the scientific plots are always standalone.
    import matplotlib, weasyprint
    report = {"html_pdf_pages": len(document.pages), "pandoc": str(pypandoc.get_pandoc_version()),
              "matplotlib": matplotlib.__version__, "weasyprint": weasyprint.__version__,
              "latex_compiled": bool(tectonic), "raw_predictions_recomputed": True,
              "tables_checked": 9, "figures": 3, "references": 8}
    (PAPER / "provenance/build_report.json").write_text(json.dumps(report, indent=2) + "\n")
    pypandoc.convert_file("supplement.md", "latex", outputfile="supplement.tex", extra_args=["--standalone"])
    if tectonic:
        completed = subprocess.run([tectonic, "--keep-logs", "--outdir", "tex-build", "supplement.tex"], check=True, text=True, capture_output=True)
        (PAPER / "provenance/supplement_latex_build.txt").write_text(completed.stdout + completed.stderr)
        shutil.copyfile(PAPER / "tex-build/supplement.pdf", PAPER / "supplement.pdf")
        for name in ("manuscript", "supplement"):
            info = subprocess.check_output(["pdfinfo", str(PAPER / (name + ".pdf"))], text=True)
            report[name + "_pdf_pages"] = int(re.search(r"Pages:\s+(\d+)", info).group(1))
        report['tectonic'] = subprocess.check_output([tectonic, '--version'], text=True).strip()
        (PAPER / "provenance/build_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


def artifact_hashes() -> None:
    lines = []
    destination = PAPER / "provenance/artifact.sha256"
    for path in sorted(PAPER.rglob("*")):
        if path.is_file() and path != destination and not ({"__pycache__", "artifacts", "release", "tex-build"} & set(path.relative_to(PAPER).parts)):
            lines.append(f"{sha(path)}  {path.relative_to(PAPER)}")
    destination.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true", help="Read and verify existing assets; no dependencies or output writes")
    parser.add_argument("--tectonic", help="Path to Tectonic; compile the actual LaTeX PDF")
    args = parser.parse_args()
    if not args.verify_only:
        for folder in ("tables", "figures", "provenance"):
            (PAPER / folder).mkdir(parents=True, exist_ok=True)
    elif not (PAPER / "provenance/evidence_snapshot.json").exists():
        raise SystemExit("Build the paper first to capture the ledger evidence snapshot.")
    data = extract()
    validate(data)
    if args.verify_only:
        expected = json.loads((PAPER / "tables/reported_results.json").read_text())
        assert data == expected, "Extracted values differ from rendered data"
        manifest = json.loads((PAPER / "provenance/source_manifest.json").read_text())
        assert all(sha(ROOT / p) == h for p, h in manifest["sources"].items()), "A source has changed"
        assert manifest["evidence_snapshot_sha256"] == sha(PAPER / "provenance/evidence_snapshot.json")
        for name, records in data.items():
            with (PAPER / "tables" / f"{name}.csv").open() as handle:
                actual = list(csv.DictReader(handle))
            assert actual == [{k: str(v) for k, v in row.items()} for row in records], name
        for name in re.findall(r"!\[.*?\]\(([^)]+)\)", (PAPER / "manuscript.md").read_text()):
            assert (PAPER / name).is_file(), name
        for line in (PAPER / "provenance/artifact.sha256").read_text().splitlines():
            digest, name = line.split("  ", 1)
            assert sha(PAPER / name) == digest, ("Draft artifact changed", name)
        print("PASS: 9 verified result tables, cohort arithmetic, paired intervals, citation keys, source hashes, prediction-lock structure, and figure links")
        return
    write_tables(data)
    from build_supplement import build
    build(data)
    provenance()
    figures(data)
    render(args.tectonic)
    artifact_hashes()


if __name__ == "__main__":
    main()
