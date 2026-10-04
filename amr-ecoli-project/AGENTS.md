# Repository Agent Instructions

These instructions apply to the entire repository.

## Required Project Skills

- For substantial tasks, read and follow `skills/manage-amr-research/SKILL.md`.
- For cohort/AST work, also read `skills/curate-amr-cohort/SKILL.md`.
- For Unity/Slurm work, also read `skills/run-amr-unity/SKILL.md`.
- For reviews, metrics, models, or claims, also read
  `skills/audit-amr-benchmarks/SKILL.md`.

## Durable State

Before planning or editing, run:

```bash
python skills/manage-amr-research/scripts/status.py check
python skills/manage-amr-research/scripts/status.py show-next
```

Before finishing substantial work, attach verified evidence to the active phase,
render `docs/PROJECT_STATUS.md`, and rerun the status check. Do not mark phases
complete from prose or intention alone.

## Scientific Integrity

- Never use legacy CARD-window or 40-isolate smoke metrics as phenotype evidence.
- Never tune on the external cohort.
- Never fill absent AST method, site, breakpoint, lineage, or cluster with proxies.
- Preserve source snapshots, exclusions, hashes, tool/database versions, splits,
  and model/evaluation manifests.
- Keep software tests, development analyses, locked external results,
  replication, and laboratory validation as distinct evidence tiers.

## Parallel Work

Give agents non-overlapping file ownership. Only the coordinating agent updates
the canonical plan, external-cohort definition, `workflow/Snakefile`,
`config/study.json`, and project status unless ownership is explicitly delegated.
