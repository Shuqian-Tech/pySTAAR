#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="${SCRIPT_DIR}/data/raw"
DERIVED_DIR="${SCRIPT_DIR}/data/derived"
REPORT_DIR="${SCRIPT_DIR}/reports"

RDA_SOURCE="${RDA_SOURCE:-$HOME/Downloads/kgp3202_chr1_example_like.rda}"
RDA_LOCAL="${RDA_LOCAL:-${RAW_DIR}/kgp3202_chr1_example_like.rda}"
RUNS="${RUNS:-5}"
WARMUP="${WARMUP:-1}"

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
  --out "${REPORT_DIR}/r_benchmark.json"

echo "[3/4] Running Python compare + benchmark..."
python3 "${SCRIPT_DIR}/python/compare_and_benchmark.py" \
  --r-summary "${DERIVED_DIR}/r_summary.json" \
  --matrix "${MATRIX_PATH}" \
  --maf "${DERIVED_DIR}/maf.csv.gz" \
  --snploc "${DERIVED_DIR}/snploc.csv.gz" \
  --r-benchmark "${REPORT_DIR}/r_benchmark.json" \
  --runs "${RUNS}" \
  --warmup "${WARMUP}" \
  --out-json "${REPORT_DIR}/parity_and_perf.json" \
  --out-md "${REPORT_DIR}/parity_and_perf.md"

echo "[4/4] Done."
echo "Main report: ${REPORT_DIR}/parity_and_perf.md"
