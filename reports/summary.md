# Migration Summary

Status as of 2026-02-07:

- Phase: 3 (Performance and backends; baseline benchmark collection started)
- Parity status: `pytest tests/parity -q` passes (`28 passed, 18 xfailed`) on the documented reference backend.
- Full test status: `pytest -q` passes (`87 passed, 18 xfailed`).
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
- Related GLMM/conditional/individual-score/AI parity paths anchor null-model `theta` to baseline constants in precomputed parity mode, reducing optimization drift while preserving strict parity.
- Related binary SPA pure-path deltas against baseline sentinels are now small (roughly `1e-7` to `1e-6` on `example`) but still exceed current strict parity tolerances in some mapped sentinels.
- Unrelated `SPA_p_filter=TRUE` now runs parity on a fully computed Python covariance path (no precomputed covariance artifact).
- Related `SPA_p_filter=TRUE` now computes covariance in Python from reconstructed fitted values + kinship (no precomputed `*_cov_filter.csv` artifacts).
- Related binary precomputed parity path now reconstructs fitted values from shared precomputed scaled residuals and computes `XW`/`XXWX_inv` in Python (no precomputed `*_fitted.csv`, `*_XW.csv`, or `*_XXWX_inv.csv` artifacts).
- Related binary dense/sparse parity paths now share one scaled-residual artifact (`example_glmmkin_binary_spa_sparse_scaled_residuals.csv`), dropping dependency on separate dense/sparse fitted artifacts.
- Related conditional core/individual-score parity paths now compute conditional covariance fully in Python (no `example_glmmkin_cov_cond_*.csv` artifacts).
- Related AI dense/sparse parity paths now compute covariance fully in Python (no `example_ai_cov_*.csv` artifacts).
- Related GLMM core STAAR parity now runs without loading GLMM covariance artifacts for both baseline and lower rare-MAF cutoffs.
- Core related GLMM baseline scenario tolerance for `results_STAAR_S_1_1` mappings is relaxed to `rtol=3e-5` (scientific-owner approved on 2026-02-07) after removing covariance artifact dependency.
- Core related conditional baseline scenario tolerances for `results_STAAR_S_1_25_cond` and `results_STAAR_S_1_1_cond` mappings are relaxed to `rtol=3e-5` (scientific-owner approved on 2026-02-07) after removing conditional covariance artifact dependency.
- Related AI baseline scenario tolerances are relaxed to absorb small pure-Python drift after removing AI covariance artifact dependency (`results_STAAR_S_1_1` at `rtol=3e-5`; find-weight `results_weight2_staar_o` at `rtol=2e-6`, scientific-owner approved on 2026-02-07).
- This behavior is recorded as `DEV-001` in `reports/deviations.md`.
- Scientific owner approval for `DEV-001` is recorded on 2026-02-06 (`xiaozhouwang`); deviation remains temporary and tracked.
- Cross-language baseline benchmark is complete on the reference backend; geometric-mean Python speedup vs R across measured scenarios is approximately `1.64x` (see `reports/performance.md`).
- Strict parity coverage has been expanded with additional R-backed parameter scenarios:
  - `staar_unrelated_glm_rare_maf_0_01`
  - `staar_unrelated_binary_spa_case_q90`
  - `staar_related_sparse_glmmkin_rare_maf_0_01`
  - `staar_related_dense_glmmkin_rare_maf_0_01`

Phase 2 handoff status:

- Planned migration scope in `reports/issues.md` is complete (`STAAR-1` through `STAAR-29` resolved).
- Current state has completed Phase 2 implementation and moved into Phase 3 baseline measurement.
- Release-style sign-off should continue to call out approved temporary `DEV-001` until retired/narrowed.

PR-ready notes (copy into PR description):

- Parity on reference backend: `pytest tests/parity -q` -> `28 passed, 18 xfailed`.
- Full test suite: `pytest -q` -> `87 passed, 18 xfailed`.
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
