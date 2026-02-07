# Changelog

All notable changes to this project are documented in this file.

## [Unreleased]

### Added

- Chinese-first documentation set under `docs/`, including quickstart, installation, tutorials, API references, and R migration guide.
- Output field reference table: `docs/api/output_fields.md`.
- Release handoff document: `reports/release_handoff.md`.

### Changed

- Related baseline parity specs now run pure-path (`use_precomputed_artifacts: false`).
- Approved tolerance policy for affected related parity sentinels relaxed up to `rtol=5e-4` (tracked in `DEV-001`).
- Migration checklist completion documented through `STAAR-41`.

### Notes

- Functional migration scope is complete.
- Strict near-bitwise R parity is not claimed; see `reports/deviations.md` for approved deviation details.
