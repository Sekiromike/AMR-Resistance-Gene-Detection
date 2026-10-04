# Reproducing the completed MIC comparison

Release: [v1.0.0](https://github.com/Sekiromike/AMR-Resistance-Gene-Detection/releases/tag/v1.0.0).
Commands below assume the `amr-ecoli-project` directory within the checkout.

## Download and verify

Download these release assets together with `SHA256SUMS`:

| Asset | Contents |
|---|---|
| `amr-results-v1.0.0.tar.gz` | Cohort, exclusions, metadata, splits, allele/family features, sequence units and mappings, frozen predictions, evaluation and model-selection manifests, environment exports, final-run logs, selected source snapshots and normalized inputs |
| `amr-embeddings-v1.0.0.tar.gz` | 64 frozen Evo 2 embedding shards |
| `paper-and-tables-v1.0.0.tar.gz` | Manuscript, supplement, figures, CSV tables, source and audit reports |
| `arxiv-source-v1.0.0.tar.gz` | Minimal LaTeX upload files and vector PDF figures |
| `research-code-v1.0.0.tar.gz` | Project source at the release commit |

```bash
sha256sum --check SHA256SUMS
mkdir -p paper/artifacts
tar -xzf amr-results-v1.0.0.tar.gz -C paper/artifacts --strip-components=1
tar -xzf amr-embeddings-v1.0.0.tar.gz -C paper/artifacts --strip-components=1
```

Both research assets share the `amr-research-artifacts/` top directory. Large
raw read collections, model weights, and the AMRFinderPlus database are accessed
through upstream manifests rather than bundled. Selected original source and
normalized inputs are under `source_inputs/data/` inside the results asset.

## Verify numerical results

```bash
uv run --no-project --with-requirements paper/verify-requirements.txt \
  python paper/scripts/verify_results.py
```

This checks recorded artifact hashes and the 18 prediction locks, reproduces
all 54 archived evaluation files and 438 paired comparisons (including their
bootstrap intervals), verifies 192 source-exclusion summaries, independently
checks EA with scalar arithmetic, and repeats the external-variant and
withholding read-outs. It writes the full-precision exports in
`paper/verification/`. It does not train or select any model.

`release_manifest.json` records original and distributed hashes for each
artifact. Only private runtime paths and cloud-query project identifiers are
redacted in public copies. The verifier checks distributed bytes first and
uses the retained original hash when a manifest refers to a redacted file.
Every numerical prediction CSV is distributed byte-for-byte unchanged.
Public-copy hashes differ for redacted manifests; this is explicit, not a
replacement of historical evidence.

## Locate the final analyses

| Analysis | Artifact directory |
|---|---|
| Final cohort | `results/cohort/64948707/` |
| Features, splits, constant and nearest-relative baselines | `results/models/64960056/` |
| Final bounded ridge development | `results/models/65006995/` |
| Final XGBoost development | `results/models/64983779/` |
| Panel/acquired-unit ablation | `results/models/65006996/` |
| Source sensitivity | `results/comparisons/65016028/` |
| Original locked external predictions | `results/external/predictions/65019273/` |
| Original external scores | `results/external/evaluation/65020832/` |
| Post-hoc family predictions | `results/posthoc/family/65022960/` |
| Post-hoc family scores | `results/posthoc/family/65022960/scoring-65023163/` |
| Variant withholding | `results/leave-variant-out/65051616/` |
| Selected units | `results/genomics/loci/64948707/` |
| Frozen embeddings | `results/genomics/embeddings/loci/64948707/evo2_20b/` |
| Pretraining-overlap audit | `results/genomics/pretraining-overlap/64994989/` |

Earlier development runs are retained for provenance and labelled by run ID.
Use `tables/final_development.csv` for the final model matrix;
`verification/all_evaluations.csv` also includes superseded development fits.
Family results are exploratory despite their separate prediction lock.

## Rebuild the manuscript

```bash
uv run --no-project --with-requirements paper/render-requirements.txt \
  python paper/scripts/build_paper.py --tectonic /path/to/tectonic
python3 paper/scripts/build_paper.py --verify-only
```

The recorded build uses Tectonic 0.17.0, Pandoc 3.9, and pinned Python rendering
packages. Tectonic downloads its TeX support bundle on first use. The generated
LaTeX contains resolved references and uses vector figures. Rendering can vary
slightly across TeX environments; inspect the submission system's final PDF.

## Model or upstream regeneration

This release makes exact downstream scoring possible without GPUs or private
cluster access. To refit development models, use `scripts/run_mic_models.py`,
`run_censored_regression.py`, `predict_external.py`, and `leave_variant_out.py`
with the archived feature, split, cohort, and embedding inputs. The exact
original commands and resource requests are in `hpc/unity/` and the archived
run manifests/logs; replace private path placeholders with local locations.
Original environments are recorded in `provenance/*explicit.txt` and
`workflow/envs/`. Keep any rerun separate from the frozen published outputs.

The external labels are now public. Repeating this release is a reproduction,
not a new untouched external validation. A revised model needs a new external
cohort for a confirmatory claim. The family model was designed after inspection
of the original external result.

Production embedding manifests omitted a checkpoint-content hash. The cached
revision recovered during release is recorded in
`provenance/checkpoint_recovery.json`; original shard hashes identify the actual
released embeddings. No cross-environment embedding-equivalence test is claimed.

## Software checks

The existing repository suite passed all 406 tests. The command and result are
recorded in `provenance/test_report.json`. These tests establish software checks;
the scientific numerical verification is separately recorded in
`verification/audit.json`.

The current study is on the [`main` branch](https://github.com/Sekiromike/AMR-Resistance-Gene-Detection/tree/main/amr-ecoli-project). The repository history was replaced at the owner's request on 4 October 2026. [study_chronology.json](provenance/study_chronology.json) preserves the current study's original plan and prediction-lock contents, hashes, and local Git dates. These records are not independent timestamps or public preregistration. Release checksums identify the current downloadable packages.
