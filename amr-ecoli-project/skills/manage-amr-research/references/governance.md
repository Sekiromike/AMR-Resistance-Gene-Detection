# Research Governance

## Evidence Tiers

- `SMOKE_TEST_ONLY`: exercises software; no scientific inference.
- `DEVELOPMENT_ONLY`: development estimates; no external-generalization claim.
- `LOCKED_EXTERNAL`: model and threshold frozen before one external evaluation.
- `INDEPENDENT_REPLICATION`: a separate source and analysis confirms the result.
- `LAB_VALIDATED`: targeted experimental validation supports a mechanism.

Never promote a claim past its evidence tier.

## Phase Rules

- A phase is `complete` only when all exit criteria have evidence paths or hashes.
- A phase is `blocked` only for a concrete external dependency; record it and
  continue independent work.
- Reopen a phase when an upstream cohort, split, breakpoint, code, or database
  hash changes.
- External results are append-only. A model changed after external evaluation is
  a new analysis and needs a new untouched cohort for a confirmatory claim.

## Agent Coordination

- Assign exclusive file ownership or review-only work.
- Require commands, tests, sources, limitations, and integration hooks.
- Integrate centrally; parallel agents do not redefine the question or split.
