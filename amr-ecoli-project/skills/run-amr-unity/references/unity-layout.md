# Unity Layout And Guardrails

- Login alias: private local SSH configuration, referenced as `<ssh-alias>`
- User and allocation: never committed
- Durable root: `${AMR_WORK_ROOT}`
- Project: `${AMR_PROJECT_ROOT:-${AMR_WORK_ROOT}/amr-ecoli-project}`
- Controller environment: `${AMR_CONTROLLER_ENV:-${AMR_WORK_ROOT}/conda/envs/amr-ecoli-controller}`
- Rule environments: `${AMR_WORK_ROOT}/conda/snakemake`
- General CPU partition: `cpu`
- Short/preemptible CPU partition: `cpu-preempt`
- General GPU partition: `gpu`
- Short/preemptible GPU partition: `gpu-preempt`

Use `sbatch` or the Snakemake Slurm executor for compute. Request the minimum
credible CPU, memory, and wall time. Use `--qos=short` only for one bounded job
under four hours. Capture `sacct -j JOBID` after completion.

GPUs are rule-scoped, never profile defaults. Keep curation, QC, annotation,
pangenome/unitig generation, GWAS, and classical modeling on CPU. Request one
GPU for verified GPU-capable foundation-model embedding or fine-tuning. Use
`gpu-preempt` only for checkpointed work under two hours; otherwise use `gpu`
for nonpreemptible jobs up to 48 hours. Start with the least powerful device that
meets measured compute-capability and VRAM needs. Run
`sbatch hpc/unity/gpu-smoke.sbatch` before the first model allocation and retain
its device record.

Official documentation:

- https://docs.unity.rc.umass.edu/documentation/jobs/
- https://docs.unity.rc.umass.edu/documentation/jobs/sbatch/
- https://docs.unity.rc.umass.edu/documentation/cluster_specs/storage/
- https://docs.unity.rc.umass.edu/documentation/managing-files/hpc-workspace/
- https://docs.unity.rc.umass.edu/documentation/software/conda/
- https://docs.unity.rc.umass.edu/documentation/tools/gpus/
