# Release Handoff

Status date: 2026-02-07

## Scope

This handoff covers full functional migration of the planned STAAR workflows (`STAAR-1` through `STAAR-49`) with parity validation on the reference backend.

## Current Validation Status

- Parity suite: `pytest tests/parity -q` -> `50 passed`
- Full suite: `pytest -q` -> `115 passed`
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
- `STAAR-47` adds distinct non-clone non-`example` R-baseline parity via `dataset=nonexample601` and scenario `staar_unrelated_glm_nonexample601`.
- `STAAR-48` extends distinct non-`example` parity into related workflows via `staar_related_sparse_glmmkin_nonexample601`.
- `STAAR-49` extends distinct non-`example` parity into related binary SPA via `staar_related_sparse_binary_spa_nonexample601`.

## Suggested Post-Release Follow-up

- `STAAR-50`: add distinct non-`example` R-baseline parity for a related conditional workflow family (for example `staar_related_sparse_glmmkin_cond`).
