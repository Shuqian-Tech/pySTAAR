# Phase 3 Non-Example Probe

- Generated: 2026-02-07T08:33:44Z
- Scenario: `staar_unrelated_glm` on runtime non-example directory dataset copy
- Dataset source: runtime copy of `data/example_*` files into directory layout (`geno.mtx`, etc.)
- Warm-up policy: 1 warm-up run(s) discarded
- Measured runs: 5
- Reported statistic: median seconds
- Reference backend details: `reports/reference_backend.md`
- Python environment details: `reports/python_environment.md`

## Result

| Scenario | Median (s) | Min (s) | Max (s) |
|---|---:|---:|---:|
| staar_unrelated_glm_nonexample_dir | 0.165067 | 0.153121 | 0.169546 |

## Notes

- This probe is used for STAAR-45 non-example performance coverage tracking.
