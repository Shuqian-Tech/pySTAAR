# Phase 3 Performance Snapshot (After STAAR-46)

- Generated: 2026-02-07T07:57:20Z
- Dataset fingerprint: `baselines/example_fingerprint.json`
- Warm-up policy: 1 warm-up run(s) discarded
- Measured runs: 5 per scenario
- Reported statistic: median seconds
- Reference backend details: `reports/reference_backend.md`
- Python environment details: `reports/python_environment.md`

## Snapshot Results

| Scenario | Median (s) | Min (s) | Max (s) | Speedup vs baseline |
|---|---:|---:|---:|---:|
| staar_unrelated_glm | 0.160165 | 0.141690 | 0.163169 | 1.00x |
| staar_related_sparse_glmmkin_pure | 1.189287 | 1.175614 | 1.206127 | 1.00x |
| staar_unrelated_binary_spa | 0.167775 | 0.151349 | 0.173336 | 1.00x |
| staar_related_sparse_binary_spa_pure | 0.909528 | 0.888602 | 0.970494 | 1.00x |
| staar_unrelated_glm_cond | 0.204292 | 0.182871 | 0.323650 | 1.00x |
| indiv_score_unrelated_glm | 0.054707 | 0.052240 | 0.069593 | 1.00x |
| ai_staar_unrelated_glm | 0.575461 | 0.542046 | 0.595572 | 1.00x |
| ai_staar_related_sparse_glmmkin_find_weight_pure | 1.814874 | 1.776217 | 1.898069 | 1.00x |

## Notes

- This report captures the post-optimization runtime snapshot after `STAAR-46`.
- Target scenario comparison is recorded in `benchmarks/phase3_opt_46_comparison.csv`.
- `staar_related_sparse_binary_spa_pure` median improved from `2.343567s` (post-`STAAR-42` probe) to `0.909528s` (`~2.58x`).
