# Reference Backend

- Document date: 2026-02-07
- Purpose: Reference backend for Phase 2 parity validation in this repository.

## Backend Identity

- OS: macOS 15.6.1 (build 24G90)
- Kernel: Darwin 24.6.0 (arm64)
- CPU: Apple M4
- Logical CPU cores: 10
- Physical CPU cores: 10
- RAM: 17179869184 bytes (16 GiB)
- GPU: Not used for parity workflows

## R Baseline Backend (Extraction)

From baseline artifacts:

- Hardware capture: `baselines/hardware.txt`
- R environment capture: `baselines/r_environment.json`
- Baseline source metadata: `baselines/SOURCE.md`

Key baseline identifiers:

- R version: 4.5.0
- BLAS (R): OpenBLAS (from `baselines/r_environment.json`)
- Baseline seed: 600

## Python Parity Backend

- Python: 3.13.5 (Anaconda build)
- NumPy: 2.1.3
- SciPy: 1.15.3
- pandas: 2.2.3
- PyYAML: 6.0.2
- pytest: 8.3.4
- BLAS/LAPACK (NumPy): OpenBLAS 0.3.21
- Threading environment variables:
  - `OMP_NUM_THREADS` unset
  - `OPENBLAS_NUM_THREADS` unset
  - `MKL_NUM_THREADS` unset
  - `VECLIB_MAXIMUM_THREADS` unset

Detailed Python capture is recorded in `reports/python_environment.md`.

## Determinism Settings

- Scenario seed: 600 (from `specs/*.yaml`)
- Workflow entry points set `np.random.seed(seed)` before computation.

## Validation Commands

- `pytest tests/parity -q`
- `pytest -q`

Latest run (2026-02-07):

- `pytest tests/parity -q` -> `28 passed, 18 xfailed`
- `pytest -q` -> `91 passed, 18 xfailed`
