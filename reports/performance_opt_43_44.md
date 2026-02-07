# Phase 3 Performance Snapshot (After STAAR-43/44)

- Generated: 2026-02-07T07:00:33Z
- Dataset fingerprint: `baselines/example_fingerprint.json`
- Warm-up policy: 1 warm-up run(s) discarded
- Measured runs: 5 per scenario
- Reported statistic: median seconds
- Reference backend details: `reports/reference_backend.md`
- Python environment details: `reports/python_environment.md`

## Snapshot Results

| Scenario | Median (s) | Min (s) | Max (s) | Speedup vs baseline |
|---|---:|---:|---:|---:|
| staar_unrelated_glm | 0.165556 | 0.150045 | 0.214349 | 1.00x |
| staar_related_sparse_glmmkin_pure | 0.916582 | 0.902721 | 0.983989 | 1.00x |
| staar_unrelated_binary_spa | 0.161430 | 0.145308 | 0.176931 | 1.00x |
| staar_related_sparse_binary_spa_pure | 2.531439 | 2.511083 | 2.643385 | 1.00x |
| staar_unrelated_glm_cond | 0.210886 | 0.199111 | 0.279937 | 1.00x |
| indiv_score_unrelated_glm | 0.092078 | 0.059392 | 0.103564 | 1.00x |
| ai_staar_unrelated_glm | 0.856711 | 0.637385 | 1.157548 | 1.00x |
| ai_staar_related_sparse_glmmkin_find_weight_pure | 1.512776 | 1.452528 | 1.736953 | 1.00x |

## Notes

- This report captures the post-optimization runtime snapshot after `STAAR-43` and `STAAR-44`.
- Comparison to the original Phase 3 baseline is recorded in `benchmarks/phase3_opt_43_44_comparison.csv`.
