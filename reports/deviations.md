# Deviation Log

## DEV-001: Pure-Path Related-Scenario Tolerance Relaxation

- Status: Closed (historical record; retired from active release gating on 2026-02-07)
- Date recorded: 2026-02-06
- Date closed: 2026-02-07
- Affected scenario families:
  - Related GLMM: `staar_related_*_glmmkin*`
  - Related conditional GLMM: `staar_related_*_glmmkin_cond*`
  - Related individual score: `indiv_score_related_*`
  - Related binary SPA: `staar_related_*_binary_spa*`
  - Related AI-STAAR: `ai_staar_related_*`
- Affected code:
  - `src/pystaar/workflows.py`
  - `src/pystaar/models.py`
  - `src/pystaar/staar_core.py`
  - related scenario specs under `specs/`

### What Changed

Related GLMM, conditional, individual-score, AI, and binary-SPA workflows now run parity on fully computed Python paths (`use_precomputed_artifacts: false` in related baseline specs). Compatibility artifacts remain available only for optional runtime fallback and are not required by parity tests.

Initial pure-path closure used relaxed tolerances up to `rtol=5e-4` and tracked 18 pure-shadow scenarios as `xfail`. In `STAAR-42`, REML estimation in `estimate_tau_reml` was aligned to a Powell-first optimization path (with Nelder-Mead fallback), related-scenario tolerances were tightened to `rtol=3.5e-4`, and all pure-shadow `xfail` blocks were removed.

### Why

Pure-Python GLMM(kinship) estimation and sparse linear algebra produce reproducible numeric drift versus the R/GMMAT extraction path on this backend. The deviation was introduced to keep parity coverage active while removing precomputed-artifact dependencies from related workflows.

### Evidence

- Previous parity state (2026-02-07, before `STAAR-42` closure): `pytest tests/parity -q` -> `28 passed, 18 xfailed`.
- Current parity state (2026-02-07, after `STAAR-49`): `pytest tests/parity -q` -> `50 passed`.
- Current full-suite state (2026-02-07): `pytest -q` -> `115 passed`.
- Related pure-shadow scenarios now run as required parity checks (no `xfail`) under the tightened `rtol=3.5e-4` envelope.

### Scientific Impact

- No change to scientific intent or user-facing workflow definitions.
- Numeric drift remains bounded by explicit, scenario-level tolerances.
- Unaffected deterministic sentinels (for example `num_variant`, `cMAC`, stable mapping keys/shapes) remain strict.

### Acceptability Criteria and Closure

- Closure criteria for `STAAR-42` were met on 2026-02-07:
  - tightened related pure-path tolerances (`5e-4` -> `3.5e-4`)
  - removed pure-shadow `xfail` gating
  - achieved fully passing parity/full test suites on the reference backend
- This entry remains as historical traceability for tolerance-policy evolution.

### Approval Record

- Scientific owner approval: **Approved**.
- Scientific owner: `xiaozhouwang`.
- Approval reference (initial deviation): migration planning/approval thread on 2026-02-06.
- Additional scientific-owner approval: 2026-02-07 (`xiaozhouwang`) for related-scenario tolerance relaxations and pure-path parity operation without parity-spec precomputed artifacts.

### Handoff Record

- Release handoff artifact: `reports/release_handoff.md`.
