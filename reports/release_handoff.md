# Release Handoff

Status date: 2026-02-07

## Scope

This handoff covers full functional migration of the planned STAAR workflows (`STAAR-1` through `STAAR-46`) with parity validation on the reference backend.

## Current Validation Status

- Parity suite: `pytest tests/parity -q` -> `47 passed`
- Full suite: `pytest -q` -> `112 passed`
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
- `STAAR-46` optimization reduced `staar_related_sparse_binary_spa` pure-path median runtime from `2.343567s` to `0.909528s` (`~2.58x`) versus the post-`STAAR-42` probe.
- `STAAR-45` adds non-`example` parity/performance coverage via runtime directory-dataset scenario `staar_unrelated_glm_nonexample_dir`.

## Suggested Post-Release Follow-up

- `STAAR-47`: add a distinct (non-clone) non-`example` R baseline extraction and matching parity scenario family for stronger release gating beyond runtime directory-copy coverage.
