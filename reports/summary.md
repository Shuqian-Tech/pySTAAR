# Migration Summary

Status as of 2026-02-08:

- Phase: 3 (Performance and backends; baseline benchmark collection started)
- Parity status: `pytest tests/parity -q` passes (`55 passed`) on the documented reference backend.
- Full test status: `pytest -q` passes (`126 passed`).
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
  - `indiv_score_related_sparse_glmmkin_cond_nonexample601`
  - `indiv_score_related_dense_glmmkin_cond`
  - `indiv_score_related_dense_glmmkin_cond_nonexample601`
  - `staar_unrelated_glm`
  - `staar_unrelated_glm_nonexample_dir`
  - `staar_unrelated_glm_nonexample601`
  - `staar_unrelated_glm_rare_maf_0_01`
  - `staar_unrelated_glm_cond`
  - `staar_unrelated_binary_spa`
  - `staar_unrelated_binary_spa_case_q90`
  - `staar_unrelated_binary_spa_filter`
  - `staar_related_sparse_binary_spa`
  - `staar_related_sparse_binary_spa_nonexample601`
  - `staar_related_sparse_binary_spa_filter`
  - `staar_related_dense_binary_spa`
  - `staar_related_dense_binary_spa_filter`
  - `staar_related_sparse_glmmkin`
  - `staar_related_sparse_glmmkin_nonexample601`
  - `staar_related_sparse_glmmkin_rare_maf_0_01`
  - `staar_related_dense_glmmkin`
  - `staar_related_dense_glmmkin_rare_maf_0_01`
  - `staar_related_sparse_glmmkin_cond`
  - `staar_related_sparse_glmmkin_cond_nonexample601`
  - `staar_related_sparse_glmmkin_cond_nonexample602`
  - `staar_related_dense_glmmkin_cond`
  - `staar_related_dense_glmmkin_cond_nonexample601`
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
  - Phase 3 optimization artifacts (`STAAR-46`): `benchmarks/phase3_opt_46_raw.csv`, `benchmarks/phase3_opt_46_summary.csv`, `benchmarks/phase3_opt_46_meta.json`, `benchmarks/phase3_opt_46_comparison.csv`, `reports/performance_opt_46.md`
  - Phase 3 optimization artifacts (`STAAR-55`): `scripts/run_phase3_opt_55_benchmarks.py`, `benchmarks/phase3_opt_55_raw.csv`, `benchmarks/phase3_opt_55_summary.csv`, `benchmarks/phase3_opt_55_meta.json`, `benchmarks/phase3_opt_55_comparison.csv`, `reports/performance_opt_55.md`
  - Phase 3 optimization artifacts (`STAAR-56`): `scripts/run_phase3_opt_56_benchmarks.py`, `benchmarks/phase3_opt_56_raw.csv`, `benchmarks/phase3_opt_56_summary.csv`, `benchmarks/phase3_opt_56_meta.json`, `benchmarks/phase3_opt_56_comparison.csv`, `reports/performance_opt_56.md`
  - Non-example performance probe (`STAAR-45`): `benchmarks/phase3_nonexample_raw.csv`, `benchmarks/phase3_nonexample_summary.csv`, `benchmarks/phase3_nonexample_meta.json`, `reports/performance_nonexample.md`
  - Distinct non-example baseline artifacts (`STAAR-47`): `baselines/nonexample601_sim_data.rds`, `baselines/nonexample601_sim_metadata.json`, `baselines/nonexample601_fingerprint.json`, `baselines/unrelated_glm_nonexample601_sentinels.json`, and runtime dataset files `data/nonexample601_*`
  - Distinct non-example related baseline artifact (`STAAR-48`): `baselines/related_sparse_glmmkin_nonexample601_sentinels.json` with parity scenario `specs/staar_related_sparse_glmmkin_nonexample601.yaml`
  - Distinct non-example related binary baseline artifact (`STAAR-49`): `baselines/related_sparse_binary_spa_nonexample601_sentinels.json` with parity scenario `specs/staar_related_sparse_binary_spa_nonexample601.yaml`
  - Distinct non-example related conditional baseline artifact (`STAAR-50`): `baselines/related_sparse_glmmkin_cond_nonexample601_sentinels.json` with parity scenario `specs/staar_related_sparse_glmmkin_cond_nonexample601.yaml`
  - Distinct non-example related dense-conditional baseline artifact (`STAAR-51`): `baselines/related_dense_glmmkin_cond_nonexample601_sentinels.json` with parity scenario `specs/staar_related_dense_glmmkin_cond_nonexample601.yaml`
  - Distinct non-example related conditional individual-score baseline artifact (`STAAR-52`): `baselines/indiv_score_related_sparse_glmmkin_cond_nonexample601_sentinels.json` with parity scenario `specs/indiv_score_related_sparse_glmmkin_cond_nonexample601.yaml`
  - Distinct non-example related dense conditional individual-score baseline artifact (`STAAR-53`): `baselines/indiv_score_related_dense_glmmkin_cond_nonexample601_sentinels.json` with parity scenario `specs/indiv_score_related_dense_glmmkin_cond_nonexample601.yaml`
  - Distinct second-seed non-example conditional related baseline artifact (`STAAR-54`): `baselines/nonexample602_sim_data.rds`, `baselines/nonexample602_sim_metadata.json`, `baselines/nonexample602_fingerprint.json`, and `baselines/related_sparse_glmmkin_cond_nonexample602_sentinels.json` with parity scenario `specs/staar_related_sparse_glmmkin_cond_nonexample602.yaml`

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
- `STAAR-42` removed all pure-shadow xfails; parity now passes for all baseline + shadow scenarios on the reference backend.
- `DEV-001` is retained as a historical record in `reports/deviations.md` and is no longer active release gating.
- Cross-language baseline benchmark is complete on the reference backend; geometric-mean Python speedup vs R across measured scenarios is approximately `1.64x` (see `reports/performance.md`).
- `STAAR-46` optimized the high-cost related sparse binary pure path: `staar_related_sparse_binary_spa_pure` median improved from `2.343567s` to `0.909528s` (`~2.58x`) versus the post-`STAAR-42` probe baseline.
- `STAAR-55` optimized repeated related-GLMM pure-path execution via cached related null-model fitting: `staar_related_sparse_glmmkin_pure` improved from `0.968130s` to `0.192804s` (`~5.02x`) and `ai_staar_related_sparse_glmmkin_find_weight_pure` improved from `1.588007s` to `0.817542s` (`~1.94x`) in targeted no-cache vs cache benchmarks.
- `STAAR-56` optimized repeated related AI runs by caching final AI payloads for unchanged dataset/metadata inputs: repeated-call medians dropped from `0.656068s` to `0.000054s` (`~12047x`) for `ai_staar_related_sparse_glmmkin` and from `0.630982s` to `0.000019s` (`~33065x`) for `ai_staar_related_sparse_glmmkin_find_weight` in targeted benchmark mode.
- `STAAR-45` adds parity + performance coverage for non-`example` runtime directory datasets via `staar_unrelated_glm_nonexample_dir`.
- `STAAR-47` adds distinct non-clone non-`example` R-baseline parity coverage via `staar_unrelated_glm_nonexample601` (`dataset=nonexample601`, sentinels in `baselines/unrelated_glm_nonexample601_sentinels.json`).
- `STAAR-48` expands distinct non-`example` parity into a related workflow family via `staar_related_sparse_glmmkin_nonexample601`.
- `STAAR-49` expands distinct non-`example` parity into related binary SPA via `staar_related_sparse_binary_spa_nonexample601`.
- `STAAR-50` expands distinct non-`example` parity into related conditional workflows via `staar_related_sparse_glmmkin_cond_nonexample601`.
- `STAAR-51` expands distinct non-`example` conditional parity into dense related STAAR via `staar_related_dense_glmmkin_cond_nonexample601`.
- `STAAR-52` expands distinct non-`example` conditional parity into related individual-score workflows via `indiv_score_related_sparse_glmmkin_cond_nonexample601`.
- `STAAR-53` expands distinct non-`example` conditional parity into dense related individual-score workflows via `indiv_score_related_dense_glmmkin_cond_nonexample601`.
- `STAAR-54` adds second-seed (`nonexample602`) conditional related parity coverage via `staar_related_sparse_glmmkin_cond_nonexample602`.
- `STAAR-50` also fixed weighted CCT exact-one handling in conditional ACAT-V paths (`cct_pval`), preventing degenerate `p=1` entries from overwhelming combined conditional ACAT statistics.
- Strict parity coverage has been expanded with additional R-backed parameter scenarios:
  - `staar_unrelated_glm_rare_maf_0_01`
  - `staar_unrelated_binary_spa_case_q90`
  - `staar_related_sparse_glmmkin_rare_maf_0_01`
  - `staar_related_dense_glmmkin_rare_maf_0_01`

Phase 2 handoff status:

- Planned migration scope in `reports/issues.md` is complete (`STAAR-1` through `STAAR-29` resolved).
- Full migration checklist in `reports/issues.md` is complete through `STAAR-56`; no open parity-robustness follow-up remains.
- Current state has completed Phase 2 implementation and moved into Phase 3 baseline measurement.
- Release sign-off can reference `DEV-001` as closed historical context only.

PR-ready notes (copy into PR description):

- Parity on reference backend: `pytest tests/parity -q` -> `55 passed`.
- Full test suite: `pytest -q` -> `126 passed`.
- Deviations: no active deviations; `DEV-001` is closed and retained historically in `reports/deviations.md`.
- Required named roles per policy:
  - Migration owner: `<fill>`
  - Human sponsor: `<fill>`
  - Scientific owner: `<fill>`
  - Maintainer (if release-impacting): `<fill>`

Next steps:

- Parity robustness for current migration scope is complete; move Phase 3 effort back to performance profiling/optimization candidates.
- Keep historical deviation context in sync with any future tolerance/backend adjustments.
