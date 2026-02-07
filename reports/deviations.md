# Deviation Log

## DEV-001: Pure-Path Parity Requires Relaxed Tolerances

- Status: Approved (release-documented, active)
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

Related GLMM, conditional, individual-score, AI, and binary-SPA workflows now default to fully computed Python paths. Baseline parity specs now run with `use_precomputed_artifacts: false` for related workflows and validate the pure-Python path directly.

The following precomputed artifacts remain available only for optional compatibility mode when `use_precomputed_artifacts=TRUE` is explicitly requested at runtime (not used by baseline parity specs):

- `data/example_glmmkin_scaled_residuals.csv`
- `data/example_glmmkin_binary_spa_sparse_scaled_residuals.csv`

These optional artifacts are not needed for CI parity and are retained only as a temporary fallback path.

For related binary-SPA parity paths, Python now reconstructs fitted values (`fitted = Y - scaled_residuals`) from a shared precomputed scaled-residual artifact and computes `XW`, `XXWX_inv`, and `SPA_p_filter` covariance in Python; precomputed `*_cov_filter.csv`, `*_fitted.csv`, `*_XW.csv`, and `*_XXWX_inv.csv` artifacts are no longer loaded. Dense and sparse related-binary parity now share a single scaled-residual artifact (`example_glmmkin_binary_spa_sparse_scaled_residuals.csv`).
Core related GLMM STAAR parity now computes covariance fully in Python for all cutoffs; baseline parity no longer loads `example_glmmkin_cov.csv`.
To keep baseline-cutoff pure-path parity (`rare_maf_cutoff=0.05`) after removing artifact dependencies, affected related scenario sentinels are relaxed up to `rtol=5e-4` with scientific-owner approval.
Core related `staar_*_glmmkin_cond` and related individual conditional score-test parity now compute conditional covariance fully in Python (no `example_glmmkin_cov_cond_*.csv` artifact usage).
This includes related core GLMM, conditional GLMM, individual score-test, related binary SPA, and related AI scenarios now exercised without precomputed parity artifacts.
Related AI sparse/dense parity now compute covariance fully in Python (no `example_ai_cov_*.csv` artifact usage).
Tolerances for unaffected sentinels (for example, `num_variant`, `cMAC`, and deterministic weight-matrix mappings) remain strict.

### Why

Pure-Python GLMM(kinship) estimation and sparse linear algebra produce reproducible numeric differences against the R/GMMAT path on this backend. Those differences are large enough to fail strict `rtol=1e-6` checks across multiple related scenarios, so parity now uses an approved relaxed tolerance envelope (up to `rtol=5e-4`) while staying on pure-Python execution paths.

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
- Related baseline parity specs have been switched to pure-path (`use_precomputed_artifacts: false`) and pass in CI.
- Largest observed pure-path relative drift in related scenarios is on the order of `~4.7e-4` (for example, related AI `results_STAAR_S_1_1` mappings) and is covered by the approved `rtol=5e-4` ceiling.
- Current related binary pure-path deltas against baseline sentinels (`example`):
  - Sparse `results_STAAR_B`: `0.23360463525923016` vs baseline `0.2336049736705653` (delta `-3.3861133513779507e-07`)
  - Dense `results_STAAR_B`: `0.23360475727667856` vs baseline `0.233605092772099` (delta `-3.354954204448646e-07`)
  - Sparse filter `results_STAAR_B`: `0.5771872755012124` vs baseline `0.5771882813570353` (delta `-1.0058558228553949e-06`)
  - Dense filter `results_STAAR_B`: `0.5771878385971415` vs baseline `0.5771888444992558` (delta `-1.0059021143815627e-06`)
- Largest observed related-binary mapped sentinel delta:
    - `results_STAAR_B_1_25["STAAR-B(1,25)"]`: `~1.199e-06` in filter scenarios.
- Related AI-STAAR pure-Python path showed mapped sentinel drift above stricter tolerances, now addressed by approved relaxed tolerances.

Parity test status with current pure-path baseline specs:

- `pytest tests/parity -q` -> `28 passed, 18 xfailed`

### Scientific Impact

- No algorithmic intent change to the test statistic definitions.
- For the `example` parity scenarios, outputs are now produced by pure-Python related workflows; no parity spec depends on precomputed intermediate artifacts.
- Generalization risk remains numeric: results can show backend-dependent drift relative to R baselines, bounded in CI by approved tolerances.

### Acceptability Criteria

- Temporary acceptance for pure-path parity closure on the `example` scenarios with approved tolerance relaxations up to `rtol=5e-4` in affected related sentinels.
- Remaining work is to tighten tolerances (or re-baseline with scientific-owner approval) if stricter equivalence is required for release gating.

### Approval Record

- Scientific owner approval: **Approved**.
- Scientific owner: `xiaozhouwang`.
- Approval reference: migration planning/approval thread on 2026-02-06 (instruction to record DEV-001 approval and proceed to Phase 3).
- Approval date: 2026-02-06.
- Approved scope: Phase 2/3 parity acceptance for current `example` scenarios on pure-Python related paths with approved tolerance relaxations.
- Merge/release gating note: approved as a temporary deviation; keep explicit mention in release notes until retired or narrowed.
- Additional scientific-owner approval: 2026-02-07 (`xiaozhouwang`) to accept pure-path related-scenario tolerance relaxations up to `rtol=5e-4` for affected parity sentinels and proceed without parity-spec precomputed artifacts.

### Handoff Record

- Final strict-parity documentation + release handoff completed: 2026-02-07.
- Handoff artifact: `reports/release_handoff.md`.

### PR Handoff Checklist

- PR status label for this deviation: `Approved temporary deviation`.
- PR description should explicitly state:
  - `DEV-001` is approved as temporary and remains tracked.
  - Parity pass status is based on pure-path related workflows with approved relaxed tolerances.
  - Follow-up work remains planned to tighten tolerances or re-baseline as needed.
