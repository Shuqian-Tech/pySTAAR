# Migration Summary

Status as of 2026-02-07:

- Phase: 3 (Performance and backends; baseline benchmark collection started)
- Parity status: `pytest tests/parity -q` passes (`46 passed`) on the documented reference backend.
- Full test status: `pytest -q` passes (`111 passed`).
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
  - `staar_unrelated_glm_rare_maf_0_01`
  - `staar_unrelated_glm_cond`
  - `staar_unrelated_binary_spa`
  - `staar_unrelated_binary_spa_case_q90`
  - `staar_unrelated_binary_spa_filter`
  - `staar_related_sparse_binary_spa`
  - `staar_related_sparse_binary_spa_filter`
  - `staar_related_dense_binary_spa`
  - `staar_related_dense_binary_spa_filter`
  - `staar_related_sparse_glmmkin`
  - `staar_related_sparse_glmmkin_rare_maf_0_01`
  - `staar_related_dense_glmmkin`
  - `staar_related_dense_glmmkin_rare_maf_0_01`
  - `staar_related_sparse_glmmkin_cond`
  - `staar_related_dense_glmmkin_cond`
- Artifacts:
  - Issues: `reports/issues.md`
  - Deviations: `reports/deviations.md`
  - Release handoff: `reports/release_handoff.md`
  - API contract map: `reports/api_contract.md`
  - Reference backend: `reports/reference_backend.md`
  - Python environment capture: `reports/python_environment.md`
  - Data source/fingerprints/checksums: `data/DATA_SOURCE.md`
  - Phase 3 performance report: `reports/performance.md` (Python vs R baseline)
  - Phase 3 Python benchmark logs: `benchmarks/phase3_baseline_raw.csv`, `benchmarks/phase3_baseline_summary.csv`, `benchmarks/phase3_baseline_meta.json`
  - Phase 3 R benchmark logs: `benchmarks/phase3_baseline_r_raw.csv`, `benchmarks/phase3_baseline_r_summary.csv`, `benchmarks/phase3_baseline_r_meta.json`
  - Phase 3 cross-language comparison: `benchmarks/phase3_cross_language_comparison.csv`
  - Phase 3 optimization artifacts (`STAAR-43`/`STAAR-44`): `benchmarks/phase3_opt_43_44_raw.csv`, `benchmarks/phase3_opt_43_44_summary.csv`, `benchmarks/phase3_opt_43_44_meta.json`, `benchmarks/phase3_opt_43_44_comparison.csv`, `reports/performance_opt_43_44.md`
  - Post-`STAAR-42` optimization targeting probe: `benchmarks/phase3_post42_probe_raw.csv`, `benchmarks/phase3_post42_probe_summary.csv`, `benchmarks/phase3_post42_probe_meta.json`, `reports/performance_post42_probe.md`

Open compliance notes:

- Related GLMM/conditional/individual-score/AI/binary-SPA workflows run baseline parity on fully computed Python paths (`use_precomputed_artifacts: false` in related baseline specs).
- Related binary SPA pure-path deltas against baseline sentinels are small (roughly `1e-7` to `1e-6` on `example`) but can exceed strict `1e-6` checks in mapped sentinels.
- Unrelated `SPA_p_filter=TRUE` now runs parity on a fully computed Python covariance path (no precomputed covariance artifact).
- Related `SPA_p_filter=TRUE` now computes covariance in Python from reconstructed fitted values + kinship (no precomputed `*_cov_filter.csv` artifacts).
- Related binary optional precomputed path reconstructs fitted values from shared precomputed scaled residuals and computes `XW`/`XXWX_inv` in Python (no precomputed `*_fitted.csv`, `*_XW.csv`, or `*_XXWX_inv.csv` artifacts).
- Related conditional core/individual-score parity paths now compute conditional covariance fully in Python (no `example_glmmkin_cov_cond_*.csv` artifacts).
- Related AI dense/sparse parity paths now compute covariance fully in Python (no `example_ai_cov_*.csv` artifacts).
- Related GLMM core STAAR parity now runs without loading GLMM covariance artifacts for both baseline and lower rare-MAF cutoffs.
- Affected related baseline and pure-shadow scenario sentinels are now tightened to `rtol=3.5e-4` (from `rtol=5e-4`) after `STAAR-42` numeric-alignment updates.
- `STAAR-42` removed all pure-shadow xfails; parity now passes for all 46 scenarios on the reference backend.
- `DEV-001` is retained as a historical record in `reports/deviations.md` and is no longer active release gating.
- Cross-language baseline benchmark is complete on the reference backend; geometric-mean Python speedup vs R across measured scenarios is approximately `1.64x` (see `reports/performance.md`).
- Strict parity coverage has been expanded with additional R-backed parameter scenarios:
  - `staar_unrelated_glm_rare_maf_0_01`
  - `staar_unrelated_binary_spa_case_q90`
  - `staar_related_sparse_glmmkin_rare_maf_0_01`
  - `staar_related_dense_glmmkin_rare_maf_0_01`

Phase 2 handoff status:

- Planned migration scope in `reports/issues.md` is complete (`STAAR-1` through `STAAR-29` resolved).
- Full migration checklist in `reports/issues.md` is complete (`STAAR-1` through `STAAR-44` resolved).
- Current state has completed Phase 2 implementation and moved into Phase 3 baseline measurement.
- Release sign-off can reference `DEV-001` as closed historical context only.

PR-ready notes (copy into PR description):

- Parity on reference backend: `pytest tests/parity -q` -> `46 passed`.
- Full test suite: `pytest -q` -> `111 passed`.
- Deviations: no active deviations; `DEV-001` is closed and retained historically in `reports/deviations.md`.
- Required named roles per policy:
  - Migration owner: `<fill>`
  - Human sponsor: `<fill>`
  - Scientific owner: `<fill>`
  - Maintainer (if release-impacting): `<fill>`

Next steps:

- Execute `STAAR-45`: expand performance/parity coverage beyond `example` with additional representative scenarios and exported R baselines.
- Execute `STAAR-46`: optimize `staar_related_sparse_binary_spa` pure path (post-`STAAR-42` probe median: `2.343567s` in `benchmarks/phase3_post42_probe_summary.csv`).
- Keep historical deviation context in sync with any future tolerance/backend adjustments.
