---
name: manage-amr-research
description: Govern this E. coli genotype-to-phenotype AMR research program, preserve project state across sessions, coordinate specialist agents, enforce phase exit criteria, and select the next valid task. Use for any substantial work in this repository, especially planning, status reporting, milestone changes, multi-agent work, or deciding whether a scientific result is ready to trust.
---

# Manage AMR Research

Treat `docs/PROJECT_STATUS.json` as canonical state and
`docs/PROJECT_STATUS.md` as its rendered human view.

## Start Every Substantial Turn

1. Run `python skills/manage-amr-research/scripts/status.py check`.
2. Run `python skills/manage-amr-research/scripts/status.py show-next`.
3. Inspect `git status --short`; preserve user changes.
4. Read the active phase, exit criteria, blockers, and latest evidence.
5. Update the live task plan. Work only on actions that advance an exit criterion
   or resolve a recorded blocker.

Read [governance.md](references/governance.md) before changing a scientific
claim, phase status, primary endpoint, external cohort, or evaluation design.

## Execute And Record

- Give each subagent a bounded artifact and non-overlapping file ownership.
- Keep the external cohort sealed from feature, model, threshold, and epoch
  selection.
- Record commands, hashes, counts, versions, and failure reports as evidence.
- Mark a phase complete only when every exit criterion has concrete evidence.
- Keep absent data absent. Never replace site, AST method, breakpoint version,
  lineage, or genomic clusters with convenient proxies.

After a verified milestone, update and render the state:

```bash
python skills/manage-amr-research/scripts/status.py set-phase \
  --id governance --status complete \
  --evidence "Unit tests and the tiny end-to-end workflow passed"
python skills/manage-amr-research/scripts/status.py render
```

Use `--next-action` to replace the immediate action and `--blocker` to append a
concrete blocker. Do not record aspirations as evidence.

## Finish Every Substantial Turn

1. Run focused tests, then the broadest affordable checks.
2. Update the active phase with new evidence, blockers, and next action.
3. Render the Markdown status and rerun `status.py check`.
4. Report what changed, what was verified, what remains blocked, and the single
   next scientific milestone.
