# Unity HPC Execution

Usernames, allocation names, SSH aliases, account strings, and key locations are
private deployment configuration and must not be committed. Create the ignored
local configuration on Unity:

```bash
cp config/unity.env.example config/unity.local.env
chmod 600 config/unity.local.env
# Replace the placeholders, then load it in each login shell.
source config/unity.local.env
```

Use `/work` for active jobs and versioned outputs. High-I/O, reproducible
intermediates may be staged in an allocated `/scratch/workspace` directory, but
scratch is expiring and has no snapshots. The only copy of a cohort manifest,
split manifest, breakpoint artifact, or final result must not live in scratch.

Unity uses Slurm. Compute is submitted to `cpu` by default; the smoke test uses
`cpu-preempt` plus the short QoS because it is bounded to twenty minutes. The
Snakemake controller runs on the login node only as a lightweight scheduler and
submits every compute rule through the Slurm executor plugin.

## Login Node Policy

Unity states that work on login nodes "can cause Unity to become sluggish for
the entire user base", that CPU and memory there are deliberately limited, and
that "disruptive processes on the login nodes may be terminated without
notice". Use `salloc` for interactive work and `sbatch` for everything else.

Permitted on a login node: `squeue`, `sacct`, `sinfo`, the `unity-slurm-*`
helper scripts, editing files, submitting jobs, and small reads of small files.

Never on a login node:

- building or solving a Conda environment;
- scanning, sorting, or hashing large result files;
- `du` over a large tree;
- any `nohup ... &` background process.

A Conda environment that needs CUDA cannot be solved on a login node in any
case: conda's `__cuda` virtual package is detected only where a GPU is present,
so a CUDA build fails to resolve there. Build GPU environments with
`hpc/unity/embeddings-env.sbatch`, which requests a GPU for exactly this reason.

Analysis of a large result file belongs in a job. Summarize inside the job and
copy back only the summary.

## Conda Package Cache

Conda defaults its package cache to `~/.conda/pkgs`, inside the 100 GB home
quota. Every environment job must export `CONDA_PKGS_DIRS` into the `/work`
allocation, as `hpc/unity/disjointness-env.sbatch` and
`hpc/unity/embeddings-env.sbatch` do. The cache is needed only at build time and
can be deleted to reclaim home space.

## Partitions, Time, and QoS

Verified against the current Unity documentation:

- `cpu`, `gpu`, `cpu-preempt`, and `gpu-preempt` all allow up to 14 days; the
  default is 1 hour, so always set `--time` explicitly. `--qos=long` is required
  to reach the 14-day maximum.
- A job on a preempt partition can be killed after two hours. Use a preempt
  partition only for work that is bounded under two hours or checkpointed.
- GPU constraints are `vram8|11|12|16|23|32|40|48|80|102`,
  `sm_52|61|70|75|80|86|89|90`, and model names such as `a100`, `a100-80g`,
  `l40s`, `a16`, `gh200`. Combine with `&` and offer alternatives with `|`,
  for example `--constraint=sm_80&vram40`.
- Prefer `--constraint` over naming one GPU model, so the scheduler can satisfy
  the request from a larger pool.

## Bootstrap

From the deployed project root:

```bash
bash hpc/unity/bootstrap.sh
sbatch --export=ALL hpc/unity/smoke.sbatch
squeue --me
```

The bootstrap places Conda environments and package caches under the provided
`/work` root, not the 100 GB home directory. It records an explicit resolved
controller environment in `provenance/controller-conda-explicit.txt`.

## Production Run

First acquire the two NCBI source exports, freeze the independently generated
isolate metadata, and place the official EUCAST workbook plus independently
verified machine-readable rules at the paths in `config/study.json`. Replace
the remaining `REQUIRED_BEFORE_COHORT_FREEZE` value with the archived
AMRFinderPlus database release. Then run:

```bash
bash hpc/unity/run.sh --dry-run
bash hpc/unity/run.sh
```

The cohort validator intentionally stops the DAG when the breakpoint artifact,
site, quantitative AST provenance, external membership, lineage, or genomic
near-neighbor cluster is absent or inconsistent.

## GPU Allocation

Do not allocate GPUs to AST acquisition/normalization, breakpoint
interpretation, cohort construction, QUAST, MLST, AMRFinderPlus, pangenome or
unitig construction, GWAS, or classical models. These are CPU workloads.

Use a GPU only for measured GPU-capable workloads such as sequence foundation
model embedding or fine-tuning. Before the first such run, submit:

```bash
sbatch --export=ALL hpc/unity/gpu-smoke.sbatch
```

The preflight requests one GPU on `gpu-preempt` for 20 minutes with the portable
minimum constraints `sm_70&vram16`, then records allocation and device metadata
under `provenance/`. For checkpointed jobs below two hours, use `gpu-preempt`.
For jobs that cannot be preempted and run up to 48 hours, use `gpu`. Request
`--gpus=1` unless scaling evidence justifies more, and choose the lowest
`vram<N>`/compute-capability constraint that fits a measured peak. Foundation
model rules must record GPU model, CUDA/driver/framework versions, peak memory,
checkpoint hash, and wall time.

`AMR_WORK_ROOT` is mandatory in every submitted environment. Do not add a
private `#SBATCH --account`, username, allocation path, or credential file to a
tracked script. Supply any site-required account at submission time or in a
private Slurm configuration outside the repository.

## Cluster References

- Unity job submission: https://docs.unity.rc.umass.edu/documentation/jobs/
- Unity batch jobs: https://docs.unity.rc.umass.edu/documentation/jobs/sbatch/
- Unity storage: https://docs.unity.rc.umass.edu/documentation/cluster_specs/storage/
- Unity scratch workspaces: https://docs.unity.rc.umass.edu/documentation/managing-files/hpc-workspace/
- Unity Conda: https://docs.unity.rc.umass.edu/documentation/software/conda/
- Unity GPUs: https://docs.unity.rc.umass.edu/documentation/tools/gpus/
- Snakemake Slurm executor: https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html
