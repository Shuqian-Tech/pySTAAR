# Release Handoff

Status date: 2026-02-07

## Scope

This handoff covers full functional migration of the planned STAAR workflows (`STAAR-1` through `STAAR-44`) with parity validation on the reference backend.

## Current Validation Status

- Parity suite: `pytest tests/parity -q` -> `46 passed`
- Full suite: `pytest -q` -> `111 passed`
- Related baseline parity scenarios run on pure-Python paths (`use_precomputed_artifacts: false` in baseline related specs).

## Deviation Summary

- Active deviations: none
- Historical deviation record: `DEV-001` (`reports/deviations.md`)
- Scientific owner approvals recorded:
  - 2026-02-06 (`xiaozhouwang`)
  - 2026-02-07 (`xiaozhouwang`)
- Final related tolerance policy from `STAAR-42`: up to `rtol=3.5e-4` (with existing per-sentinel `atol` values retained).
- Unaffected deterministic sentinels (for example `num_variant`, `cMAC`, stable mapping fields) remain strict.

## Runtime Behavior for Users

- Default related workflows execute fully in Python.
- Baseline parity does not require precomputed R-derived covariance artifacts.
- Optional compatibility mode remains available via `use_precomputed_artifacts=True` in workflow APIs but is not part of baseline parity gating.

## Release Notes Guidance

Include the following points in release notes:

- Full STAAR migration scope implemented (core, conditional, AI, individual score, binary SPA).
- Reference-backend parity is validated on pure-Python related paths.
- `DEV-001` is closed and retained as historical context; related parity runs are fully passing with tightened `rtol=3.5e-4` in affected related sentinels.

## Suggested Post-Release Follow-up

- Continue profiling and optimizing high-cost related workflows while retaining current parity gates.
- Expand parity/performance validation beyond `example` with additional representative datasets/scenarios.
