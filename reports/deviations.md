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

- `data/example_glmmkin_scaled_residuals.csv`
- `data/example_glmmkin_binary_spa_sparse_scaled_residuals.csv`

These are passed into the Python null-model/STAAR pipeline to reduce backend-specific numeric drift and match baseline sentinels for the `example` scenario.
Related precomputed parity paths also anchor GLMM null-model `theta` to baseline constants (sparse/dense) to reduce residual backend optimization drift.

For related binary-SPA parity paths, Python now reconstructs fitted values (`fitted = Y - scaled_residuals`) from a shared precomputed scaled-residual artifact and computes `XW`, `XXWX_inv`, and `SPA_p_filter` covariance in Python; precomputed `*_cov_filter.csv`, `*_fitted.csv`, `*_XW.csv`, and `*_XXWX_inv.csv` artifacts are no longer loaded. Dense and sparse related-binary parity now share a single scaled-residual artifact (`example_glmmkin_binary_spa_sparse_scaled_residuals.csv`).
Core related GLMM STAAR parity now computes covariance fully in Python for all cutoffs; baseline parity no longer loads `example_glmmkin_cov.csv`.
To keep baseline-cutoff strict parity (`rare_maf_cutoff=0.05`) after removing the covariance artifact dependency, the `results_STAAR_S_1_1` mapping tolerance in related sparse/dense core STAAR specs is relaxed from `rtol=1e-6` to `rtol=3e-5` with scientific-owner approval.
Core related `staar_*_glmmkin_cond` and related individual conditional score-test parity now compute conditional covariance fully in Python (no `example_glmmkin_cov_cond_*.csv` artifact usage).
To keep baseline-cutoff strict parity (`rare_maf_cutoff=0.05`) after removing conditional covariance artifacts, related sparse/dense conditional SKAT mapping tolerances (`results_STAAR_S_1_25_cond`, `results_STAAR_S_1_1_cond`) are relaxed from `rtol=1e-6` to `rtol=3e-5` with scientific-owner approval.
Related AI sparse/dense parity now compute covariance fully in Python (no `example_ai_cov_*.csv` artifact usage).
To keep baseline-cutoff strict parity (`rare_maf_cutoff=0.05`) after removing AI covariance artifacts, related AI `results_STAAR_S_1_1` mapping tolerances are relaxed from `rtol=1e-6` to `rtol=3e-5`, and related AI find-weight `results_weight2_staar_o` mapping tolerances are relaxed from `rtol=1e-6` to `rtol=2e-6`, with scientific-owner approval.

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
- Related binary SPA prefilter now computes covariance in Python from reconstructed fitted values + kinship (no precomputed covariance artifact needed) while preserving strict parity.
- Related binary SPA precomputed parity path now reconstructs fitted values from shared precomputed scaled residuals and derives `XW`/`XXWX_inv` in Python (no precomputed fitted/XW/XXWX_inv artifacts needed) while preserving strict parity.
- Related GLMM parity at `rare_maf_cutoff=0.01` now runs without loading any GLMM covariance artifact while preserving strict parity.
- After removing `example_glmmkin_cov.csv` from baseline-cutoff core related GLMM parity (`rare_maf_cutoff=0.05`), the largest remaining mapped drift is `results_STAAR_S_1_1["SKAT(1,1)-Z8"]` at approximately `1.12e-5` versus baseline; approved tolerance relaxation for this mapping (`rtol=3e-5`) restores parity compliance.
- After removing `example_glmmkin_cov_cond_sparse.csv` from baseline-cutoff core related conditional parity (`rare_maf_cutoff=0.05`), the largest remaining mapped drifts are `results_STAAR_S_1_25_cond["SKAT(1,25)-Z2"]` at approximately `1.72e-6` (relative) and `results_STAAR_S_1_1_cond["SKAT(1,1)-Z8"]` at approximately `1.58e-5` (relative); approved tolerance relaxation for these mappings (`rtol=3e-5`) restores parity compliance.
- Related GLMM/AI/individual-score precomputed parity paths now use baseline `theta` constants (sparse/dense) to reduce null-model fit drift while preserving strict parity.
- After removing `example_ai_cov_*.csv` artifacts from baseline-cutoff related AI parity (`rare_maf_cutoff=0.05`), the largest remaining mapped drifts are:
  - `results_STAAR_S_1_1["SKAT(1,1)-Z2"]`: approximately `1.83e-5` relative.
  - `results_weight2_staar_o["B2"]`: approximately `1.01e-6` relative.
  Approved tolerance relaxations for these mappings (`rtol=3e-5` and `rtol=2e-6`, respectively) restore parity compliance.
- Current related binary pure-path deltas against baseline sentinels (`example`):
  - Sparse `results_STAAR_B`: `0.23360463525923016` vs baseline `0.2336049736705653` (delta `-3.3861133513779507e-07`)
  - Dense `results_STAAR_B`: `0.23360475727667856` vs baseline `0.233605092772099` (delta `-3.354954204448646e-07`)
  - Sparse filter `results_STAAR_B`: `0.5771872755012124` vs baseline `0.5771882813570353` (delta `-1.0058558228553949e-06`)
  - Dense filter `results_STAAR_B`: `0.5771878385971415` vs baseline `0.5771888444992558` (delta `-1.0059021143815627e-06`)
- Largest observed related-binary mapped sentinel delta:
    - `results_STAAR_B_1_25["STAAR-B(1,25)"]`: `~1.199e-06` in filter scenarios.
- Related AI-STAAR pure-Python path showed mapped sentinel drift above baseline tolerance before using precomputed AI covariance artifacts.

Parity test status with current hybrid path:

- `pytest tests/parity -q` -> `28 passed, 18 xfailed`

### Scientific Impact

- No algorithmic intent change to the test statistic definitions.
- For the `example` parity scenarios, reported outputs are anchored to R-derived intermediates for related workflows.
- Generalization risk: workflows on datasets without corresponding precomputed artifacts may show small numeric differences relative to R baselines.

### Acceptability Criteria

- Temporary acceptance only for Phase 2 parity closure on the `example` scenarios.
- Related binary SPA default path is already fully computed in Python; remaining work is to reduce/remove remaining parity-only precomputed artifact usage (notably related binary scaled-residual artifacts and baseline-cutoff core related-GLMM scaled/theta anchoring), or explicitly re-baseline/approve.

### Approval Record

- Scientific owner approval: **Approved**.
- Scientific owner: `xiaozhouwang`.
- Approval reference: migration planning/approval thread on 2026-02-06 (instruction to record DEV-001 approval and proceed to Phase 3).
- Approval date: 2026-02-06.
- Approved scope: Phase 2 parity acceptance for current `example` scenarios using parity-spec opt-in precomputed artifacts; continue Phase 3 with explicit tracking and follow-up to retire/narrow the deviation.
- Merge/release gating note: approved as a temporary deviation; keep explicit mention in release notes until retired or narrowed.
- Additional scientific-owner approval: 2026-02-07 (`xiaozhouwang`) to accept baseline-cutoff core related GLMM tolerance relaxation for `results_STAAR_S_1_1` mappings (`rtol=3e-5`) after removing covariance artifact dependency.
- Additional scientific-owner approval: 2026-02-07 (`xiaozhouwang`) to accept baseline-cutoff core related conditional SKAT mapping tolerance relaxation (`results_STAAR_S_1_25_cond`, `results_STAAR_S_1_1_cond`, `rtol=3e-5`) after removing conditional covariance artifact dependency.
- Additional scientific-owner approval: 2026-02-07 (`xiaozhouwang`) to accept baseline-cutoff related AI mapping tolerance relaxations (`results_STAAR_S_1_1`: `rtol=3e-5`; find-weight `results_weight2_staar_o`: `rtol=2e-6`) after removing AI covariance artifact dependency.

### PR Handoff Checklist

- PR status label for this deviation: `Approved temporary deviation`.
- PR description should explicitly state:
  - `DEV-001` is approved as temporary and remains tracked.
  - Parity pass status is based on parity specs that opt in to precomputed artifacts.
  - Follow-up work remains planned to retire or narrow `DEV-001`.
