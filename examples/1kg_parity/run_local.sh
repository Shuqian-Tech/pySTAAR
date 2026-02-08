#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RAW_DIR="${SCRIPT_DIR}/data/raw"
DERIVED_DIR="${SCRIPT_DIR}/data/derived"
REPORT_DIR="${SCRIPT_DIR}/reports"

RDA_SOURCE="${RDA_SOURCE:-$HOME/Downloads/kgp3202_chr1_example_like.rda}"
RDA_LOCAL="${RDA_LOCAL:-${RAW_DIR}/kgp3202_chr1_example_like.rda}"
RUNS="${RUNS:-5}"
WARMUP="${WARMUP:-1}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
REPORT_TAG="${REPORT_TAG:-}"

if [[ -n "${REPORT_TAG}" ]]; then
  REPORT_SUFFIX="_${REPORT_TAG}"
else
  REPORT_SUFFIX=""
fi

R_BENCHMARK_OUT="${REPORT_DIR}/r_benchmark${REPORT_SUFFIX}.json"
PARITY_JSON_OUT="${REPORT_DIR}/parity_and_perf${REPORT_SUFFIX}.json"
PARITY_MD_OUT="${REPORT_DIR}/parity_and_perf${REPORT_SUFFIX}.md"

mkdir -p "${RAW_DIR}" "${DERIVED_DIR}" "${REPORT_DIR}"

if [[ ! -f "${RDA_LOCAL}" ]]; then
  if [[ ! -f "${RDA_SOURCE}" ]]; then
    echo "ERROR: RDA file not found."
    echo "Checked:"
    echo "  RDA_LOCAL=${RDA_LOCAL}"
    echo "  RDA_SOURCE=${RDA_SOURCE}"
    exit 1
  fi
  cp "${RDA_SOURCE}" "${RDA_LOCAL}"
  echo "Copied RDA to local example dir: ${RDA_LOCAL}"
else
  echo "Using local RDA: ${RDA_LOCAL}"
fi

echo "[1/4] Exporting from RDA..."
Rscript "${SCRIPT_DIR}/r/export_from_rda.R" \
  --input-rda "${RDA_LOCAL}" \
  --outdir "${DERIVED_DIR}"

MATRIX_PATH="${DERIVED_DIR}/genotype.mtx.gz"
if [[ ! -f "${MATRIX_PATH}" ]]; then
  MATRIX_PATH="${DERIVED_DIR}/genotype.mtx"
fi
if [[ ! -f "${MATRIX_PATH}" ]]; then
  echo "ERROR: exported Matrix Market file not found in ${DERIVED_DIR}"
  exit 1
fi

echo "[2/4] Running R benchmark..."
Rscript "${SCRIPT_DIR}/r/benchmark_matrix_ops.R" \
  --matrix "${MATRIX_PATH}" \
  --maf "${DERIVED_DIR}/maf.csv.gz" \
  --snploc "${DERIVED_DIR}/snploc.csv.gz" \
  --runs "${RUNS}" \
  --warmup "${WARMUP}" \
  --out "${R_BENCHMARK_OUT}"

echo "[3/4] Running Python compare + benchmark..."
echo "Using Python interpreter: ${PYTHON_BIN}"
"${PYTHON_BIN}" -c 'import sys; print("Python executable:", sys.executable); print("Python version:", sys.version.split()[0])'
PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}" \
"${PYTHON_BIN}" "${SCRIPT_DIR}/python/compare_and_benchmark.py" \
  --r-summary "${DERIVED_DIR}/r_summary.json" \
  --matrix "${MATRIX_PATH}" \
  --maf "${DERIVED_DIR}/maf.csv.gz" \
  --snploc "${DERIVED_DIR}/snploc.csv.gz" \
  --r-benchmark "${R_BENCHMARK_OUT}" \
  --runs "${RUNS}" \
  --warmup "${WARMUP}" \
  --out-json "${PARITY_JSON_OUT}" \
  --out-md "${PARITY_MD_OUT}"

echo "[4/4] Done."
echo "Main report: ${PARITY_MD_OUT}"
