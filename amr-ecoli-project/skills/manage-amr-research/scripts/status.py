"""Validate, display, update, and render the canonical project phase ledger."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STATE_PATH = REPO_ROOT / "docs" / "PROJECT_STATUS.json"
MARKDOWN_PATH = REPO_ROOT / "docs" / "PROJECT_STATUS.md"
ALLOWED_STATUSES = {"pending", "in_progress", "complete", "blocked"}


def load_state() -> dict:
    return json.loads(STATE_PATH.read_text(encoding="utf-8"))


def validate_state(state: dict) -> list[str]:
    errors: list[str] = []
    for key in ("schema_version", "objective", "updated_at_utc", "current_phase", "phases"):
        if key not in state:
            errors.append(f"Missing top-level key: {key}")
    phases = state.get("phases", [])
    ids = [phase.get("id") for phase in phases]
    if len(ids) != len(set(ids)):
        errors.append("Phase IDs are not unique.")
    in_progress = []
    for phase in phases:
        phase_id = phase.get("id", "<missing>")
        status = phase.get("status")
        if status not in ALLOWED_STATUSES:
            errors.append(f"Phase {phase_id} has invalid status: {status}")
        if status == "in_progress":
            in_progress.append(phase_id)
        if not phase.get("exit_criteria"):
            errors.append(f"Phase {phase_id} has no exit criteria.")
        if status == "complete" and not phase.get("evidence"):
            errors.append(f"Complete phase {phase_id} has no evidence.")
    if len(in_progress) > 1:
        errors.append(f"More than one phase is in progress: {in_progress}")
    current = state.get("current_phase")
    if current not in ids:
        errors.append(f"current_phase does not name a phase: {current}")
    elif in_progress and current != in_progress[0]:
        errors.append(f"current_phase={current} but in-progress phase={in_progress[0]}")
    return errors


def active_phase(state: dict) -> dict:
    current = state["current_phase"]
    return next(phase for phase in state["phases"] if phase["id"] == current)


def render_markdown(state: dict) -> str:
    phase = active_phase(state)
    lines = [
        "# Project Status",
        "",
        f"Updated: `{state['updated_at_utc']}`  ",
        f"Current phase: **{phase['title']}** (`{phase['status']}`)",
        "",
        "## Objective",
        "",
        state["objective"],
        "",
        "## Immediate Action",
        "",
        phase.get("next_action") or "No action recorded.",
        "",
        "## Phase Ledger",
        "",
        "| Phase | Status | Evidence |",
        "| --- | --- | --- |",
    ]
    for item in state["phases"]:
        evidence_count = len(item.get("evidence", []))
        lines.append(f"| {item['title']} | `{item['status']}` | {evidence_count} item(s) |")
    blockers = [
        (item["title"], blocker)
        for item in state["phases"]
        for blocker in item.get("blockers", [])
    ]
    lines.extend(["", "## Blockers", ""])
    if blockers:
        lines.extend(f"- **{title}:** {blocker}" for title, blocker in blockers)
    else:
        lines.append("- None recorded.")
    lines.extend(["", "## Current Exit Criteria", ""])
    lines.extend(f"- {criterion}" for criterion in phase["exit_criteria"])
    lines.extend(["", "## Latest Evidence", ""])
    evidence = [
        (item["title"], record)
        for item in state["phases"]
        for record in item.get("evidence", [])
    ][-12:]
    if evidence:
        for title, record in evidence:
            lines.append(f"- `{record['at_utc']}` **{title}:** {record['detail']}")
    else:
        lines.append("- No evidence recorded.")
    lines.extend(
        [
            "",
            "This file is generated from `docs/PROJECT_STATUS.json`. Update it with",
            "`python skills/manage-amr-research/scripts/status.py`.",
            "",
        ]
    )
    return "\n".join(lines)


def save_state(state: dict) -> None:
    errors = validate_state(state)
    if errors:
        raise ValueError("Invalid state:\n- " + "\n- ".join(errors))
    STATE_PATH.write_text(json.dumps(state, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    MARKDOWN_PATH.write_text(render_markdown(state), encoding="utf-8")


def choose_current_phase(state: dict) -> str:
    for phase in state["phases"]:
        if phase["status"] == "in_progress":
            return phase["id"]
    for phase in state["phases"]:
        if phase["status"] == "pending":
            return phase["id"]
    return state["phases"][-1]["id"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("check")
    subparsers.add_parser("show-next")
    subparsers.add_parser("render")
    update = subparsers.add_parser("set-phase")
    update.add_argument("--id", required=True)
    update.add_argument("--status", choices=sorted(ALLOWED_STATUSES))
    update.add_argument("--evidence")
    update.add_argument("--next-action")
    update.add_argument("--blocker")
    update.add_argument(
        "--clear-blockers",
        action="store_true",
        help="Remove resolved blockers from the selected phase before adding a new blocker",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    state = load_state()
    if args.command == "check":
        errors = validate_state(state)
        if errors:
            for error in errors:
                print(f"FAIL: {error}")
            raise SystemExit(1)
        print(
            f"Project state PASS | current={state['current_phase']} | "
            f"updated={state['updated_at_utc']}"
        )
        return
    if args.command == "show-next":
        phase = active_phase(state)
        print(f"{phase['id']}: {phase['title']} [{phase['status']}]")
        print(phase.get("next_action") or "No action recorded.")
        if phase.get("blockers"):
            print("Blockers:")
            for blocker in phase["blockers"]:
                print(f"- {blocker}")
        return
    if args.command == "render":
        MARKDOWN_PATH.write_text(render_markdown(state), encoding="utf-8")
        print(f"Rendered {MARKDOWN_PATH.relative_to(REPO_ROOT)}")
        return

    phase = next((item for item in state["phases"] if item["id"] == args.id), None)
    if phase is None:
        raise SystemExit(f"Unknown phase: {args.id}")
    now = datetime.now(timezone.utc).isoformat()
    if args.status:
        phase["status"] = args.status
    if args.evidence:
        phase.setdefault("evidence", []).append({"at_utc": now, "detail": args.evidence})
    if args.next_action:
        phase["next_action"] = args.next_action
    if args.clear_blockers:
        phase["blockers"] = []
    if args.blocker and args.blocker not in phase.setdefault("blockers", []):
        phase["blockers"].append(args.blocker)
    state["updated_at_utc"] = now
    state["current_phase"] = choose_current_phase(state)
    save_state(state)
    print(f"Updated phase {args.id}; current phase is {state['current_phase']}")


if __name__ == "__main__":
    main()
