# API Stability and Deprecation Policy (1.x)

This policy defines what `pySTAAR` considers stable in the `1.x` series.

## Stability Scope

The following are treated as stable public API in `1.x`:

- Top-level symbols exported by `pystaar.__all__`.
- Workflow function signatures in `src/pystaar/workflows.py`.
- Documented output sentinel keys used in parity specs under `specs/`.

Internal helpers prefixed with `_` are not part of the public contract and may change without notice.

## SemVer Expectations

- Patch releases (`1.0.x`): bug fixes and performance improvements; no intentional breaking API changes.
- Minor releases (`1.x.0`): backward-compatible additions only.
- Major releases (`2.0+`): breaking API changes are allowed with migration documentation.

## Deprecation Process

If a public symbol or behavior must change:

1. Introduce a backward-compatible transition path where possible.
2. Document deprecation in `CHANGELOG.md` with first deprecation version.
3. Keep deprecated behavior for at least one minor release before removal.
4. Provide explicit migration guidance in docs (`docs/migration_from_r.md` and/or API docs).

## Runtime Cache Behavior

`pySTAAR` uses runtime caches for repeated workflow calls. This behavior is considered operationally stable:

- Cache inspection: `pystaar.get_runtime_cache_info()`
- Cache clearing: `pystaar.clear_runtime_caches(include_dataset_cache=True)`

Users can disable cache accumulation between jobs by calling `clear_runtime_caches()`.
