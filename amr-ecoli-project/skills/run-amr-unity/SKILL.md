---
name: run-amr-unity
description: Deploy, execute, monitor, and verify this AMR Snakemake workflow on UMass Unity HPC. Use for ssh unity work, Slurm profiles and jobs, Conda or Apptainer environments, project synchronization, storage placement, resource requests, smoke tests, provenance capture, and recovery from failed cluster runs.
---

# Run AMR On Unity

Read [unity-layout.md](references/unity-layout.md) before submitting work.

## Procedure

1. Inspect `docs/PROJECT_STATUS.json` and local changes before synchronization.
2. Require `AMR_WORK_ROOT` to name the private durable `/work` allocation; never
   commit the allocation, username, account, SSH configuration, or key path.
3. Keep source, frozen manifests, breakpoint artifacts, and final outputs in
   `/work`; use allocated scratch only for reproducible temporary I/O.
4. Probe tools before installing. Use `module load conda/latest`; never modify
   system Python.
5. Run `bash hpc/unity/bootstrap.sh`, then submit `hpc/unity/smoke.sbatch`
   with `--export=ALL` so the private environment reaches the job.
6. Inspect `sacct`, the complete Slurm log, generated hashes, and environment
   provenance before declaring success.
7. Dry-run production with `bash hpc/unity/run.sh --dry-run`; then submit through
   the Snakemake Slurm executor. Do not run compute on a login node.
8. Record job IDs, tree state, environment export, database releases, inputs,
   and outputs in the project ledger.

Never delete or release scratch workspaces automatically. Never cancel jobs or
overwrite remote results without explicit scope confirmation.
