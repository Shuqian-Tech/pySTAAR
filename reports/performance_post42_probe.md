# Phase 3 Performance Baseline

- Generated: 2026-02-07T07:46:15Z
- Dataset fingerprint: `baselines/example_fingerprint.json`
- Warm-up policy: 1 warm-up run(s) discarded
- Measured runs: 5 per scenario
- Reported statistic: median seconds
- Reference backend details: `reports/reference_backend.md`
- Python environment details: `reports/python_environment.md`

## Baseline Results

| Scenario | Median (s) | Min (s) | Max (s) | Speedup vs baseline |
|---|---:|---:|---:|---:|
| staar_unrelated_glm | 0.166589 | 0.086462 | 0.207595 | 1.00x |
| staar_related_sparse_glmmkin_pure | 1.103220 | 1.074899 | 1.133861 | 1.00x |
| staar_unrelated_binary_spa | 0.150487 | 0.143573 | 0.233007 | 1.00x |
| staar_related_sparse_binary_spa_pure | 2.343567 | 2.323226 | 2.401468 | 1.00x |
| staar_unrelated_glm_cond | 0.253488 | 0.216589 | 0.287521 | 1.00x |
| indiv_score_unrelated_glm | 0.059947 | 0.057400 | 0.066295 | 1.00x |
| ai_staar_unrelated_glm | 0.587930 | 0.565201 | 0.619542 | 1.00x |
| ai_staar_related_sparse_glmmkin_find_weight_pure | 1.768918 | 1.735303 | 1.867076 | 1.00x |

## Notes

- This report captures baseline Phase 3 performance before optimization.
- Relative speedups are 1.00x by definition until optimized variants are added.
