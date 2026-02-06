# Deviation Log

## DEV-001: Parity Uses Precomputed R-Derived Intermediates

- Status: Approved (temporary)
- Date recorded: 2026-02-06
- Affected scenarios:
  - `indiv_score_related_sparse_glmmkin`
  - `indiv_score_related_dense_glmmkin`
  - `indiv_score_related_sparse_glmmkin_cond`
  - `indiv_score_related_dense_glmmkin_cond`
  - `staar_related_sparse_glmmkin`
  - `staar_related_sparse_glmmkin_rare_maf_0_01`
  - `staar_related_dense_glmmkin`
  - `staar_related_dense_glmmkin_rare_maf_0_01`
  - `staar_related_sparse_glmmkin_cond`
  - `staar_related_dense_glmmkin_cond`
  - `staar_related_sparse_binary_spa`
  - `staar_related_sparse_binary_spa_filter`
  - `staar_related_dense_binary_spa`
  - `staar_related_dense_binary_spa_filter`
  - `ai_staar_related_sparse_glmmkin`
  - `ai_staar_related_dense_glmmkin`
  - `ai_staar_related_sparse_glmmkin_find_weight`
  - `ai_staar_related_dense_glmmkin_find_weight`
- Affected code:
  - `src/pystaar/workflows.py`
  - `src/pystaar/models.py`
  - `src/pystaar/staar_core.py`

### What Changed

Related GLMM, conditional, individual-score, AI, and binary-SPA workflows now default to fully computed Python paths and only consume baseline precomputed artifacts when `use_precomputed_artifacts=TRUE` is explicitly requested (used by parity specs for `example`).

The following precomputed artifacts are still used for baseline parity when `use_precomputed_artifacts=TRUE` is enabled in specs:

- `data/example_glmmkin_cov.csv`
- `data/example_glmmkin_scaled_residuals.csv`
- `data/example_glmmkin_cov_cond_sparse.csv`
- `data/example_glmmkin_cov_cond_dense.csv`
- `data/example_glmmkin_binary_spa_sparse_fitted.csv`
- `data/example_ai_cov_sparse_s1_b1.csv`
- `data/example_ai_cov_sparse_s1_b2.csv`
- `data/example_ai_cov_sparse_s2_b1.csv`
- `data/example_ai_cov_sparse_s2_b2.csv`
- `data/example_ai_cov_dense_s1_b1.csv`
- `data/example_ai_cov_dense_s1_b2.csv`
- `data/example_ai_cov_dense_s2_b1.csv`
- `data/example_ai_cov_dense_s2_b2.csv`

These are passed into the Python null-model/STAAR pipeline to reduce backend-specific numeric drift and match baseline sentinels for the `example` scenario.

For related binary-SPA parity paths, Python now reconstructs `scaled_residuals`, `XW`, and `XXWX_inv` from precomputed fitted values and computes `SPA_p_filter` covariance from fitted values + kinship; precomputed `*_cov_filter.csv`, `*_scaled_residuals.csv`, `*_XW.csv`, and `*_XXWX_inv.csv` artifacts are no longer loaded. Dense and sparse related-binary parity now share a single fitted artifact (`example_glmmkin_binary_spa_sparse_fitted.csv`).
For related GLMM parity paths with `rare_maf_cutoff` below the baseline (`0.05`), Python now derives the covariance submatrix directly from `example_glmmkin_cov.csv`; cutoff-specific covariance artifacts are no longer loaded.

### Why

Pure-Python GLMM(kinship) estimation and sparse linear algebra produce small, reproducible numeric differences against the R/GMMAT path on this backend. Those differences are large enough to fail some strict sentinel tolerances in `specs/`.

### Evidence

Observed on 2026-02-06 (reference backend):

- Baseline (`related_sparse_glmmkin_sentinels.json`) `results_STAAR_O`: `0.1396628042746395`
- Pure-Python related path `results_STAAR_O`: `0.13966320353776843`
  - Delta: `+3.9926312891958027e-07`
- Hybrid path using precomputed artifacts `results_STAAR_O`: `0.13966280427468974`
  - Delta: `+5.0237591864288333e-14`
- Largest observed pure-Python delta in one mapped sentinel:
  - `results_STAAR_S_1_25["SKAT(1,25)-Z2"]`: `+5.2818928803877174e-05`
- Largest observed pure-Python delta in related conditional sentinel:
  - `results_STAAR_S_1_25_cond["SKAT(1,25)-Z2"]`: `-2.6103357830986607e-05`
- Related binary SPA now has a fully computed Python PQL-style fallback path and uses that path by default.
- Unrelated binary SPA prefilter now computes covariance directly from the Python null model (no precomputed covariance artifact needed) while preserving strict parity.
- Related binary SPA prefilter now computes covariance in Python from precomputed fitted values + kinship (no precomputed covariance artifact needed) while preserving strict parity.
- Related binary SPA precomputed parity path now reconstructs `scaled_residuals`, `XW`, and `XXWX_inv` from fitted values (no precomputed component artifacts needed) while preserving strict parity.
- Related GLMM parity at `rare_maf_cutoff=0.01` now derives covariance from baseline `example_glmmkin_cov.csv` (no precomputed `example_glmmkin_cov_rare_maf_0_01.csv` dependency) while preserving strict parity.
- Current related binary pure-path deltas against baseline sentinels (`example`):
  - Sparse `results_STAAR_B`: `0.23360206101344283` vs baseline `0.2336049736705653` (delta `-2.92187712246116e-06`)
  - Dense `results_STAAR_B`: `0.2336021825223668` vs baseline `0.233605092772099` (delta `-2.9102497322207713e-06`)
  - Sparse filter `results_STAAR_B`: `0.5771799157416181` vs baseline `0.5771882813570353` (delta `-8.365615417165462e-06`)
  - Dense filter `results_STAAR_B`: `0.577180478848294` vs baseline `0.5771888444992558` (delta `-8.36565096178465e-06`)
  - Largest observed related-binary mapped sentinel delta:
    - `results_STAAR_B_1_25["STAAR-B(1,25)"]`: `~9.984e-06` in filter scenarios.
- Related AI-STAAR pure-Python path showed mapped sentinel drift above baseline tolerance before using precomputed AI covariance artifacts.

Parity test status with current hybrid path:

- `pytest tests/parity -q` -> `28 passed`

### Scientific Impact

- No algorithmic intent change to the test statistic definitions.
- For the `example` parity scenarios, reported outputs are anchored to R-derived intermediates for related workflows.
- Generalization risk: workflows on datasets without corresponding precomputed artifacts may show small numeric differences relative to R baselines.

### Acceptability Criteria

- Temporary acceptance only for Phase 2 parity closure on the `example` scenarios.
- Related binary SPA default path is already fully computed in Python; remaining work is to reduce/remove remaining parity-only precomputed artifact usage (notably related binary fitted-value artifacts and related GLMM/conditional/AI covariance artifacts), or explicitly re-baseline/approve.

### Approval Record

- Scientific owner approval: **Approved**.
- Scientific owner: `xiaozhouwang`.
- Approval reference: migration planning/approval thread on 2026-02-06 (instruction to record DEV-001 approval and proceed to Phase 3).
- Approval date: 2026-02-06.
- Approved scope: Phase 2 parity acceptance for current `example` scenarios using parity-spec opt-in precomputed artifacts; continue Phase 3 with explicit tracking and follow-up to retire/narrow the deviation.
- Merge/release gating note: approved as a temporary deviation; keep explicit mention in release notes until retired or narrowed.

### PR Handoff Checklist

- PR status label for this deviation: `Approved temporary deviation`.
- PR description should explicitly state:
  - `DEV-001` is approved as temporary and remains tracked.
  - Parity pass status is based on parity specs that opt in to precomputed artifacts.
  - Follow-up work remains planned to retire or narrow `DEV-001`.
