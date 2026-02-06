# Migration Summary

Status as of 2026-02-06:

- Phase: 2 (Executable parity; implementation complete, open deviation pending approval)
- Parity status: `pytest tests/parity -q` passes (`24 passed`) on the documented reference backend.
- Full test status: `pytest -q` passes (`47 passed`).
- Scenarios implemented:
  - `ai_staar_related_sparse_glmmkin_find_weight`
  - `ai_staar_related_dense_glmmkin_find_weight`
  - `ai_staar_related_sparse_glmmkin`
  - `ai_staar_related_dense_glmmkin`
  - `ai_staar_unrelated_glm_find_weight`
  - `ai_staar_unrelated_glm`
  - `indiv_score_unrelated_glm`
  - `indiv_score_related_sparse_glmmkin`
  - `indiv_score_related_dense_glmmkin`
  - `indiv_score_unrelated_glm_cond`
  - `indiv_score_related_sparse_glmmkin_cond`
  - `indiv_score_related_dense_glmmkin_cond`
  - `staar_unrelated_glm`
  - `staar_unrelated_glm_cond`
  - `staar_unrelated_binary_spa`
  - `staar_unrelated_binary_spa_filter`
  - `staar_related_sparse_binary_spa`
  - `staar_related_sparse_binary_spa_filter`
  - `staar_related_dense_binary_spa`
  - `staar_related_dense_binary_spa_filter`
  - `staar_related_sparse_glmmkin`
  - `staar_related_dense_glmmkin`
  - `staar_related_sparse_glmmkin_cond`
  - `staar_related_dense_glmmkin_cond`
- Artifacts:
  - Issues: `reports/issues.md`
  - Deviations: `reports/deviations.md`
  - Reference backend: `reports/reference_backend.md`
  - Python environment capture: `reports/python_environment.md`
  - Data source/fingerprints/checksums: `data/DATA_SOURCE.md`

Open compliance notes:

- Related GLMM/conditional/individual-score/AI/binary-SPA workflows default to fully computed Python paths; parity scenarios opt in to baseline artifacts via `use_precomputed_artifacts: true`.
- Related binary SPA pure-path deltas against baseline sentinels are now small (roughly `1e-6` to `1e-5` on `example`) but still exceed current strict parity tolerances in some mapped sentinels.
- `SPA_p_filter=TRUE` workflows currently consume precomputed R-derived covariance artifacts for baseline parity scenarios.
- This behavior is recorded as `DEV-001` in `reports/deviations.md`.
- Scientific owner approval for `DEV-001` is pending and must be recorded before treating the deviation as approved.

Phase 2 handoff status:

- Planned migration scope in `reports/issues.md` is complete (`STAAR-1` through `STAAR-25` resolved).
- Current state is suitable for a Phase 2 handoff PR as "implementation complete with open deviation".
- Release-style sign-off remains gated by scientific owner approval for `DEV-001`.

PR-ready notes (copy into PR description):

- Parity on reference backend: `pytest tests/parity -q` -> `24 passed`.
- Full test suite: `pytest -q` -> `47 passed`.
- Deviations: `DEV-001` (open, approval pending) in `reports/deviations.md`.
- Required named roles per policy:
  - Migration owner: `<fill>`
  - Human sponsor: `<fill>`
  - Scientific owner: `<fill>`
  - Maintainer (if release-impacting): `<fill>`

Next steps:

- Retire or narrow `DEV-001` by validating parity thresholds/specs against fully computed Python paths and removing remaining parity-only precomputed covariance usage where feasible.
- Expand parity/perf coverage beyond `example` with additional representative scenarios.
