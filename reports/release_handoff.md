# Release Handoff

Status date: 2026-02-08

## Scope

This handoff covers full functional migration of the planned STAAR workflows (`STAAR-1` through `STAAR-56`) with parity validation on the reference backend.

## Current Validation Status

- Parity suite: `pytest tests/parity -q` -> `55 passed`
- Full suite: `pytest -q` -> `128 passed`
- Related baseline parity scenarios run on pure-Python paths (`use_precomputed_artifacts: false` in baseline related specs).
- Wheel release smoke check: `python scripts/run_release_smoke_checks.py` -> `PASS` (artifacts: `reports/release_readiness.md`, `reports/release_smoke.json`)

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
- Wheel builds bundle runtime datasets under `pystaar/_data`, so `dataset="example"` workflows run after normal `pip install`.

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
- `STAAR-50` extends distinct non-`example` parity into related conditional workflows via `staar_related_sparse_glmmkin_cond_nonexample601`.
- `STAAR-50` fixes weighted CCT exact-one handling in conditional ACAT aggregation (`cct_pval`) to avoid degenerate `p=1` entries forcing inflated combined p-values.
- `STAAR-51` extends distinct non-`example` conditional parity into dense related STAAR via `staar_related_dense_glmmkin_cond_nonexample601`.
- `STAAR-52` extends distinct non-`example` conditional parity into related individual-score workflows via `indiv_score_related_sparse_glmmkin_cond_nonexample601`.
- `STAAR-53` extends distinct non-`example` conditional parity into dense related individual-score workflows via `indiv_score_related_dense_glmmkin_cond_nonexample601`.
- `STAAR-54` adds second-seed (`nonexample602`) conditional related parity coverage via `staar_related_sparse_glmmkin_cond_nonexample602`.
- `STAAR-55` optimizes repeated related-GLMM pure-path runs by caching related null-model fits for string-addressable datasets; targeted no-cache vs cache benchmark shows median improvements `~5.02x` (`staar_related_sparse_glmmkin_pure`) and `~1.94x` (`ai_staar_related_sparse_glmmkin_find_weight_pure`).
- `STAAR-56` optimizes repeated related AI workflow runs by caching repeated-call AI payloads when dataset/metadata inputs are unchanged; targeted baseline-vs-cache repeated-call benchmark shows median improvements `~12047x` (`ai_staar_related_sparse_glmmkin`) and `~33065x` (`ai_staar_related_sparse_glmmkin_find_weight`).
- 1.0 release hardening includes wheel-install smoke checks and bundled package datasets so example workflows run in fresh venv installs.

## Suggested Post-Release Follow-up

- Parity robustness follow-up: none pending for current migration scope; next work can prioritize Phase 3 performance profiling and optimization.
