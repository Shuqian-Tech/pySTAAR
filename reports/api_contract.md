# API Contract (R to Python)

Status as of 2026-02-07.

This file tracks the strict API/UX contract alignment against exported functions in `../STAAR/NAMESPACE`.

## R Export Mapping

| R export | Python public symbol | Status |
|---|---|---|
| `CCT` | `pystaar.CCT` | Mapped alias to `pystaar.staar_stats.cct` |
| `fit_null_glm` | `pystaar.fit_null_glm` | Mapped |
| `fit_null_glmmkin` | `pystaar.fit_null_glmmkin` | Mapped |
| `fit_null_glm_Binary_SPA` | `pystaar.fit_null_glm_Binary_SPA` | Mapped alias to `fit_null_glm_binary_spa` |
| `fit_null_glmmkin_Binary_SPA` | `pystaar.fit_null_glmmkin_Binary_SPA` | Mapped alias to `fit_null_glmmkin_binary_spa` |
| `STAAR` | `pystaar.STAAR` | Mapped alias to `pystaar.staar_core.staar` |
| `STAAR_cond` | `pystaar.STAAR_cond` | Mapped alias to `pystaar.staar_core.staar_cond` |
| `Indiv_Score_Test_Region` | `pystaar.Indiv_Score_Test_Region` | Mapped alias to `pystaar.staar_core.indiv_score_test_region` |
| `Indiv_Score_Test_Region_cond` | `pystaar.Indiv_Score_Test_Region_cond` | Mapped alias to `pystaar.staar_core.indiv_score_test_region_cond` |
| `STAAR_Binary_SPA` | `pystaar.STAAR_Binary_SPA` | Mapped alias to `pystaar.staar_core.staar_binary_spa` |
| `AI_STAAR` | `pystaar.AI_STAAR` | Mapped alias to `pystaar.staar_core.ai_staar` |

## Current Caveats for Strict Equivalence

- Full API name mapping is in place, but strict behavioral equivalence still depends on scenario/test coverage breadth.
- Current parity remains strongest on `example` baselines, with additional non-`example` coverage via `staar_unrelated_glm_nonexample_dir` (runtime directory mode), distinct R-baselines `staar_unrelated_glm_nonexample601` and `staar_related_sparse_glmmkin_cond_nonexample602`, and distinct related baselines `staar_related_sparse_glmmkin_nonexample601`, `staar_related_sparse_binary_spa_nonexample601`, `staar_related_sparse_glmmkin_cond_nonexample601`, `staar_related_dense_glmmkin_cond_nonexample601`, `indiv_score_related_sparse_glmmkin_cond_nonexample601`, and `indiv_score_related_dense_glmmkin_cond_nonexample601` in `specs/`.
- No active parity-blocking deviations remain; `DEV-001` is closed and retained in `reports/deviations.md` as historical context for related-scenario tolerance policy.
