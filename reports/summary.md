# Migration Summary

Status as of 2026-02-06:

- Phase: 3 (Performance and backends; baseline benchmark collection started)
- Parity status: `pytest tests/parity -q` passes (`26 passed`) on the documented reference backend.
- Full test status: `pytest -q` passes (`52 passed`).
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
  - API contract map: `reports/api_contract.md`
  - Reference backend: `reports/reference_backend.md`
  - Python environment capture: `reports/python_environment.md`
  - Data source/fingerprints/checksums: `data/DATA_SOURCE.md`
  - Phase 3 performance report: `reports/performance.md` (Python vs R baseline)
  - Phase 3 Python benchmark logs: `benchmarks/phase3_baseline_raw.csv`, `benchmarks/phase3_baseline_summary.csv`, `benchmarks/phase3_baseline_meta.json`
  - Phase 3 R benchmark logs: `benchmarks/phase3_baseline_r_raw.csv`, `benchmarks/phase3_baseline_r_summary.csv`, `benchmarks/phase3_baseline_r_meta.json`
  - Phase 3 cross-language comparison: `benchmarks/phase3_cross_language_comparison.csv`

Open compliance notes:

- Related GLMM/conditional/individual-score/AI/binary-SPA workflows default to fully computed Python paths; parity scenarios opt in to baseline artifacts via `use_precomputed_artifacts: true`.
- Related binary SPA pure-path deltas against baseline sentinels are now small (roughly `1e-6` to `1e-5` on `example`) but still exceed current strict parity tolerances in some mapped sentinels.
- `SPA_p_filter=TRUE` workflows currently consume precomputed R-derived covariance artifacts for baseline parity scenarios.
- This behavior is recorded as `DEV-001` in `reports/deviations.md`.
- Scientific owner approval for `DEV-001` is recorded on 2026-02-06 (`xiaozhouwang`); deviation remains temporary and tracked.
- Cross-language baseline benchmark is complete on the reference backend; geometric-mean Python speedup vs R across measured scenarios is approximately `1.64x` (see `reports/performance.md`).
- Strict parity coverage has been expanded with additional R-backed parameter scenarios:
  - `staar_unrelated_glm_rare_maf_0_01`
  - `staar_unrelated_binary_spa_case_q90`

Phase 2 handoff status:

- Planned migration scope in `reports/issues.md` is complete (`STAAR-1` through `STAAR-25` resolved).
- Current state has completed Phase 2 implementation and moved into Phase 3 baseline measurement.
- Release-style sign-off should continue to call out approved temporary `DEV-001` until retired/narrowed.

PR-ready notes (copy into PR description):

- Parity on reference backend: `pytest tests/parity -q` -> `26 passed`.
- Full test suite: `pytest -q` -> `52 passed`.
- Deviations: `DEV-001` (approved temporary) in `reports/deviations.md`.
- Required named roles per policy:
  - Migration owner: `<fill>`
  - Human sponsor: `<fill>`
  - Scientific owner: `<fill>`
  - Maintainer (if release-impacting): `<fill>`

Next steps:

- Add optimized variants for selected workflows and measure speedups against `benchmarks/phase3_baseline_summary.csv` and `benchmarks/phase3_cross_language_comparison.csv`.
- Prioritize scenarios where Python is currently slower than R (notably `staar_unrelated_binary_spa`, `ai_staar_unrelated_glm`).
- Re-run parity after each optimization change and retain `DEV-001` tracking until retired/narrowed.
- Expand performance and parity coverage beyond `example` with additional representative scenarios.
