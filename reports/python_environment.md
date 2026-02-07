# Python Parity Environment Capture

- Capture timestamp (UTC): 2026-02-06 02:00:10 UTC
- Capture timestamp (local): 2026-02-06 10:00:10 CST

## System

- `platform.platform()`: `macOS-15.6.1-arm64-arm-64bit-Mach-O`
- `sw_vers`:
  - ProductName: macOS
  - ProductVersion: 15.6.1
  - BuildVersion: 24G90
- `uname -a`: `Darwin ... 24.6.0 ... arm64`
- CPU: Apple M4
- Logical CPUs: 10
- Physical CPUs: 10
- Memory bytes: 17179869184

## Python Runtime

- Python: `3.13.5 | packaged by Anaconda, Inc.`
- Executable family: CPython

## Package Versions Used for Parity

- numpy==2.1.3
- scipy==1.15.3
- pandas==2.2.3
- PyYAML==6.0.2
- pytest==8.3.4
- patsy==1.0.1
- statsmodels==0.14.4
- rpy2==3.6.4

## Numerical Backend

From `numpy.__config__.show()`:

- BLAS: OpenBLAS 0.3.21
- LAPACK: OpenBLAS 0.3.21
- OpenBLAS config: `USE_OPENMP=0`, `MAX_THREADS=128`

## Threading/Backend Environment Variables

- `OMP_NUM_THREADS=None`
- `OPENBLAS_NUM_THREADS=None`
- `MKL_NUM_THREADS=None`
- `VECLIB_MAXIMUM_THREADS=None`
- `NUMEXPR_NUM_THREADS=None`

## Test Execution Snapshot

- `pytest tests/parity -q` -> `54 passed`
- `pytest -q` -> `121 passed`
