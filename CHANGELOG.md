# Changelog

All notable changes to this project are documented in this file.

## [Unreleased]

### Added

- Chinese-first documentation set under `docs/`, including quickstart, installation, tutorials, API references, and R migration guide.
- Output field reference table: `docs/api/output_fields.md`.
- Release handoff document: `reports/release_handoff.md`.
- Regression tests for CCT input validation, GLM residual-dof guardrails, AI covariance fast-path activation, and `matrix_flip` all-missing behavior.
- Runtime cache management APIs: `get_runtime_cache_info` and `clear_runtime_caches`.

### Changed

- Related baseline parity specs now run pure-path (`use_precomputed_artifacts: false`).
- Approved tolerance policy for affected related parity sentinels relaxed up to `rtol=5e-4` (tracked in `DEV-001`).
- Migration checklist completion documented through `STAAR-41`.
- Migration checklist completion documented through `STAAR-56`.
- `CCT`/`CCT_pval` now validate empty inputs and zero-sum/invalid weights explicitly.
- `fit_null_glm` now raises when residual degrees of freedom are non-positive (`n_samples <= n_parameters`).
- Phase 3 optimization pass for `STAAR-43`/`STAAR-44`:
  - Faster unrelated AI covariance construction via population-group sufficient statistics.
  - Cached dataset loading for repeated workflow benchmarking.
  - Cached file-backed AI metadata loading for example workflows.
  - Vectorized `matrix_flip` and lower-overhead tail-probability helpers.

### Notes

- Functional migration scope is complete.
- Strict near-bitwise R parity is not claimed; see `reports/deviations.md` for approved deviation details.
