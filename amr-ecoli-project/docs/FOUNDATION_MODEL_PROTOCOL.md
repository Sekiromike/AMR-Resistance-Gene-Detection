# Genome Foundation Model Protocol

Status: revision 0.2 ADOPTED 2026-09-27 by the decision owner; not yet executed
Version: 0.1 (2026-09-20), 0.2 (2026-09-27)
Ladder position: rung 7 of `docs/STUDY_PROTOCOL.md`

This protocol governs sequence foundation-model representations. It is written
before any embedding exists so that the evaluation cannot be shaped by results.

## Why this runs before the cohort is rebuilt

Embedding extraction consumes **assemblies only**. It reads no MIC, no
category, no split, and no external membership. It therefore passes no gate and
skips no gate: it is phenotype-free feature computation, and it is started
early because it carries the most unresolved engineering risk in the programme
(context length over multi-megabase genomes, aggregation strategy, GPU memory,
throughput over 13,493 assemblies).

**Nothing produced here may be evaluated against phenotype until the locked
mechanistic and interpretable baselines in rungs 1-4 are complete.** Computing
a representation is not crediting a model.

## Model selection

Revised 2026-09-21. An earlier draft named Evo 1 because the existing Evo 2
work on Unity had not been found. That draft is superseded.

**Evo 2 is the model, using the environments and pipeline already built and
validated in `evo2_project` on Unity.** Both `evo2_7b` and `evo2_20b` have been
run there against an independent benchmark.

| Asset | State |
|---|---|
| `~/.conda/envs/evo2` | evo2 0.5.5, flash-attn 2.8.0.post2, torch 2.7.1 |
| `~/.conda/envs/evo2_1b` | evo2 0.5.5 + vortex, torch 2.11.0+cu130 |
| `evo2_project/scripts/01_embed_rpob_with_evo2.py` | working extraction: `Evo2(model)` -> named layer -> mean-pool over tokens |
| `evo2_project/results/` | 7B vs 20B comparison and a five-point layer sweep |

DNABERT-2 is retained only as a short-context comparator, not a second headline.

### Two environment defects, both recorded

1. `import evo2` fails with `ldconfig -p | grep 'libnvrtc'` returning
   non-zero. `transformer_engine._load_nvrtc` globs `CUDA_HOME/**/libnvrtc.so*`
   first and only then falls back to `ldconfig`, where `check_output` raises
   when grep matches nothing, so the `LD_LIBRARY_PATH` fallback is never
   reached. Unity supplies CUDA through modules that set `LD_LIBRARY_PATH` but
   do not populate the `ldconfig` cache, so the check fails even with
   `cuda/12.6` loaded. **Fix, verified in job 64663130:** export
   `CUDA_HOME=<env>/lib/python3.12/site-packages/nvidia`, whose bundled
   `libnvrtc.so.12` satisfies the glob. Import then succeeds with CUDA
   available.
2. transformer_engine reports `Supported flash-attn versions are >= 2.1.1,
   <= 2.7.4.post1. Found flash-attn 2.8.0.post2.` This is a warning at import.
   It must be revisited before trusting throughput or attention numerics, and
   the observed version is recorded with every run.

### Layer choice, and the limits of the RpoB evidence

An earlier draft of this protocol named `blocks.6.mlp.l3` primary **because it
scored highest on the RpoB phylogeny sweep**. That reasoning is withdrawn. The
sweep is useful, but not for choosing a layer on the strength of its
phylogenetic correlation.

**Why phylogenetic fidelity is the wrong selection criterion here.** The RpoB
experiment asks whether Evo 2 embeddings preserve evolutionary distance, i.e.
whether they encode lineage. For AMR prediction, lineage is the **confounder**,
not the target. Biased sampling driven by population structure is the
documented reason genotype-to-phenotype models collapse out of distribution,
and this study's entire contribution is reliability under lineage shift. A
layer chosen to maximise phylogenetic signal may be exactly the layer that best
encodes population structure, which is the failure mode we are trying to
measure and control.

Three further reasons the transfer is weak:

- **Scale.** RpoB is a single ~4 kb conserved gene. Here the unit is a ~5 Mb
  draft assembly cut into windows; composition and context are different.
- **Evolutionary dynamics.** RpoB is a core gene under purifying selection.
  AMR determinants are frequently horizontally acquired, plasmid-borne or
  transposon-associated accessory sequence. Representing core-gene phylogeny
  well implies little about detecting a recently acquired cassette.
- **Effect size.** The best correlation is r ~ 0.43, modest for a task the
  model should find easy.

**What the sweep does support.** The collapse from r = 0.433 at
`blocks.6.mlp.l3` to 0.045 at `blocks.23.mlp.l3` is most plausibly a property
of the representation rather than of the task: late layers of an autoregressive
genome model specialise toward next-token prediction and lose general-purpose
features, which matches published reports of a stability boundary in deeper
layers. That is adequate evidence to **exclude** deep layers, and nothing more.

| Layer | RpoB archaea cosine r | Use here |
|---|---:|---|
| `blocks.6.mlp.l3` | 0.433 | candidate |
| `blocks.12.mlp.l3` | 0.397 | candidate |
| `blocks.17.mlp.l3` | 0.421 | candidate |
| `blocks.21.mlp.l3` | 0.108 | excluded, degenerate |
| `blocks.23.mlp.l3` | 0.045 | excluded, degenerate |

**Revised rule.** The candidate set is `blocks.6`, `blocks.12` and
`blocks.17` (`.mlp.l3`). Selection among the three happens **inside training
folds on development data only**, as nested selection, never on the external
set and never on the phylogeny score. Extraction runs for all three candidates
because re-extraction is the expensive step; selection remains a
training-fold decision.

**Diagnostic, not objective.** Phylogenetic structure in the chosen embedding
is still worth measuring, but as a *warning sign*: if nearest-neighbour
structure in embedding space largely recapitulates lineage, that predicts poor
out-of-distribution behaviour and must be reported alongside performance rather
than treated as a quality score.

## Revision 0.2 (ADOPTED 2026-09-27): embed resistance loci, not whole genomes

**Why the whole-genome design is withdrawn.**

1. *Cost, measured.* The existing Evo 2 20B runs on Unity (jobs 56562342 and
   56562343, H100) embedded 0.42 Mb in 130 s and 1.41 Mb in 266 s, a marginal
   rate of about 7,300 bases per second after a ~70 s model load. Whole
   assemblies on both strands are about 138 Gb here: roughly 19 million GPU
   seconds, **about 5,000 GPU-hours**. That is not a feasible or proportionate
   use of the allocation.
2. *Science.* A resistance determinant is ~1 kb of a ~5 Mb genome, about
   0.02%. A mean over the whole genome is dominated by the core genome, i.e. by
   lineage, which amendments 004 and 007 exist to keep out of the estimate.
   Such a representation would mostly measure population structure.

**Proposed design.** For every eligible genome, embed a fixed set of
resistance-relevant sequences, each as its own unit:

| Unit | Content | Located by |
|---|---|---|
| Chromosomal panel | gyrA, gyrB, parC, parE (fluoroquinolone targets); marR, acrR, soxR (efflux regulators); ampC with its promoter/attenuator; ompF, ompC, ompR, envZ (porins and their regulators) | blastn of E. coli K-12 MG1655 reference loci (each with 300 bp upstream) against the assembly; best hit with ≥ 80% identity over ≥ 90% of the reference, else recorded missing |
| Acquired determinants | every AMRFinderPlus Type AMR gene call, with 500 bp flanks | AMRFinderPlus coordinates already produced |

About 60 kb per genome, 1.6 Gb on both strands, **about 60 GPU-hours** at the
measured rate. The existing rules otherwise stand: candidate layers
blocks.6/12/17 `.mlp.l3`, both strands averaged, mean pooling within a unit,
and every choice among them made inside training folds only.

**Question tested.** Per antibiotic, on identical rows and nested CV:
(a) AMRFinderPlus features alone; (b) AMRFinderPlus features plus Evo 2 locus
embeddings, reduced by PCA fitted inside training folds; (c) Evo 2 locus
embeddings alone. The foundation model is credited only if (b) beats (a) on
the paired lineage bootstrap. (c) shows how much of the known mechanism the
model recovers without the curated catalogue.

**Unchanged gates.** The pretraining-overlap audit still runs before any
claim. Evo 2 was trained on public genomes that may include these isolates;
overlap cannot leak phenotype, since Evo 2 saw no phenotypes, but it must be
reported. The external set stays sealed.

## Credit rule (adopted 2026-09-28, before all-sources Evo 2 results or any external MIC were read)

Evo 2 is credited for an antibiotic only if **all** of the following hold:

1. In the primary population (human clinical), development nested CV,
   AMRFinderPlus + Evo 2 beats AMRFinderPlus alone for **both** learners
   (ridge_censored and xgboost_aft), each paired lineage-bootstrap 95%
   interval of the essential-agreement difference excluding zero.
2. The same learner comparison, fixed before unsealing, is positive on the
   locked external evaluation with its interval excluding zero.
3. The pretraining-overlap audit is reported with the result.

All-sources results are reported as a sensitivity analysis and cannot create
credit on their own. Results failing any condition are reported as "no
demonstrated gain", not as a trend.

At adoption, human-clinical development results met condition 1 for
gentamicin only (+0.027 [+0.010, +0.038] ridge; +0.024 [+0.001, +0.043]
XGBoost). Ceftriaxone (+0.009 ridge, −0.015 XGBoost) and ciprofloxacin
(−0.003, −0.005) did not.

## Windowing and aggregation

Prespecified before extraction so aggregation cannot be tuned on phenotype:

- **Windowing:** non-overlapping windows of the model's native context across
  each assembly, in the given contig order, both strands considered.
- **Reverse complement:** the representation must be RC-consistent. Where the
  architecture is not RC-invariant, both strands are embedded and averaged, and
  the choice is recorded.
- **Aggregation:** mean pooling is the prespecified default. Any alternative
  (max, attention pooling, per-window sequence summarisation) is declared
  before it is fitted and is selected **inside training folds only**.
- **Layer choice:** declared before extraction and fixed. Layer selection by
  downstream performance is a tuning decision and is confined to training folds.

## Compute and provenance

Per `docs/UNITY_HPC.md`: start at `--gpus=1`, measure peak memory and wall time
on a bounded canary before any production array, use `gpu-preempt` with
checkpointing for short work and `gpu` for long non-preemptible runs, and
choose the lowest `vram<N>` constraint that fits the measured peak.

Every run records GPU model, CUDA/driver/framework versions, model checkpoint
revision, peak memory, wall time, window and aggregation parameters, and a
SHA-256 for every output shard. Matched compute cost against the baselines is
part of the comparison, not an afterthought.

## Execution order

1. GPU preflight (`hpc/unity/gpu-smoke.sbatch`) — device provenance.
2. Environment build and checkpoint download with a recorded revision hash.
3. **Bounded canary**: a handful of assemblies. Measure peak VRAM, throughput,
   and RC consistency. No production array before this passes.
4. Pretraining-overlap audit.
5. Production extraction across the 13,493 assemblies, checkpointed and
   resumable, writing sharded embeddings with per-shard checksums.
6. Hold. Embeddings wait, unevaluated, until rungs 1-4 are locked.

## What this protocol does not authorise

- No evaluation against any phenotype before the locked baselines exist.
- No use of the external-locked set for layer, window, or aggregation choice.
- No claim that a foundation model improves AMR prediction. That claim requires
  beating matched baselines on identical isolates with the overlap audit
  attached.
