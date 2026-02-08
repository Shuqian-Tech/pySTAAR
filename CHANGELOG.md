# Changelog

All notable changes to this project are documented in this file.

## [Unreleased]

## [1.0.0] - 2026-02-08

### Added

- Full workflow migration coverage through `STAAR-56` with parity scenarios and reports under `reports/`.
- Runtime cache management APIs for users:
  - `pystaar.get_runtime_cache_info()`
  - `pystaar.clear_runtime_caches(include_dataset_cache=True)`
- API stability/deprecation policy document: `docs/api/stability.md`.
- CI workflow matrix for Python `3.10`, `3.11`, `3.12`: `.github/workflows/ci.yml`.
- Release smoke-check script and artifacts:
  - `scripts/run_release_smoke_checks.py`
  - `reports/release_readiness.md`
  - `reports/release_smoke.json`
- Cold-process vs warm-inprocess latency benchmarking:
  - `scripts/run_phase3_cold_warm_benchmarks.py`
  - `reports/performance_cold_warm.md`
  - `benchmarks/phase3_cold_warm_*`

### Changed

- Package version bumped to `1.0.0` in `pyproject.toml`.
- Wheel packaging now bundles runtime datasets under `pystaar/_data`, and dataset resolution now supports both source-tree and installed-wheel layouts.
- Related workflow performance improvements from optimization passes:
  - `STAAR-55`: cached related null-model fits for repeated pure-path runs.
  - `STAAR-56`: cached repeated related AI workflow payloads for unchanged inputs.
- Documentation and release summaries updated for current status and test counts.

### Validation Snapshot

- Parity suite: `pytest tests/parity -q` -> `55 passed`
- Full suite: `pytest -q` -> `128 passed`
- Wheel smoke check: `python scripts/run_release_smoke_checks.py` -> `PASS`

### Notes

- Functional migration scope is complete for planned `R -> Python` parity targets.
- Strict near-bitwise parity with R is not claimed; historical deviation context is retained in `reports/deviations.md`.
