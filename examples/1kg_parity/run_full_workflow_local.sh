#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

RAW_DIR="${SCRIPT_DIR}/data/raw"
WORKFLOW_DATASET_DIR="${SCRIPT_DIR}/data/workflow_dataset"
REPORT_DIR="${SCRIPT_DIR}/reports"

RDA_SOURCE="${RDA_SOURCE:-$HOME/Downloads/kgp3202_chr1_example_like.rda}"
RDA_LOCAL="${RDA_LOCAL:-${RAW_DIR}/kgp3202_chr1_example_like.rda}"

TARGET_NSAMPLE="${TARGET_NSAMPLE:-1200}"
TARGET_NVAR="${TARGET_NVAR:-1500}"
SAMPLE_MAF_MAX="${SAMPLE_MAF_MAX:-0.2}"
NUM_ANNOTATIONS="${NUM_ANNOTATIONS:-10}"
KIN_SPARSE_THRESHOLD="${KIN_SPARSE_THRESHOLD:-0.02}"

SEED="${SEED:-600}"
RARE_MAF_CUTOFF="${RARE_MAF_CUTOFF:-0.05}"
ADJ_VARIANTS="${ADJ_VARIANTS:-1}"
RUNS="${RUNS:-3}"
WARMUP="${WARMUP:-1}"
ATOL="${ATOL:-1e-6}"
RTOL="${RTOL:-1e-3}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
REPORT_TAG="${REPORT_TAG:-}"

if [[ -n "${REPORT_TAG}" ]]; then
  REPORT_SUFFIX="_${REPORT_TAG}"
else
  REPORT_SUFFIX=""
fi

R_RESULTS_OUT="${REPORT_DIR}/r_full_workflow_results${REPORT_SUFFIX}.json"
R_BENCHMARK_OUT="${REPORT_DIR}/r_full_workflow_benchmark${REPORT_SUFFIX}.json"
FULL_JSON_OUT="${REPORT_DIR}/full_workflow_parity_and_perf${REPORT_SUFFIX}.json"
FULL_MD_OUT="${REPORT_DIR}/full_workflow_parity_and_perf${REPORT_SUFFIX}.md"

mkdir -p "${RAW_DIR}" "${WORKFLOW_DATASET_DIR}" "${REPORT_DIR}"

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

echo "[1/4] Building reproducible simulated workflow dataset..."
Rscript "${SCRIPT_DIR}/r/build_simulated_workflow_dataset.R" \
  --input-rda "${RDA_LOCAL}" \
  --outdir "${WORKFLOW_DATASET_DIR}" \
  --seed "${SEED}" \
  --target-nsample "${TARGET_NSAMPLE}" \
  --target-nvar "${TARGET_NVAR}" \
  --sample-maf-max "${SAMPLE_MAF_MAX}" \
  --num-annotations "${NUM_ANNOTATIONS}" \
  --kin-sparse-threshold "${KIN_SPARSE_THRESHOLD}"

echo "[2/4] Running R full-workflow benchmark suite..."
Rscript "${SCRIPT_DIR}/r/run_full_workflow_benchmarks.R" \
  --sim-rds "${WORKFLOW_DATASET_DIR}/sim_workflow_data.rds" \
  --out-results "${R_RESULTS_OUT}" \
  --out-benchmark "${R_BENCHMARK_OUT}" \
  --runs "${RUNS}" \
  --warmup "${WARMUP}" \
  --seed "${SEED}" \
  --rare-maf-cutoff "${RARE_MAF_CUTOFF}" \
  --adj-variants "${ADJ_VARIANTS}"

echo "[3/4] Running Python suite + R/Python comparison..."
echo "Using Python interpreter: ${PYTHON_BIN}"
"${PYTHON_BIN}" -c 'import sys; print("Python executable:", sys.executable); print("Python version:", sys.version.split()[0])'
PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}" \
"${PYTHON_BIN}" "${SCRIPT_DIR}/python/run_full_workflow_compare.py" \
  --dataset-dir "${WORKFLOW_DATASET_DIR}" \
  --r-results "${R_RESULTS_OUT}" \
  --r-benchmark "${R_BENCHMARK_OUT}" \
  --runs "${RUNS}" \
  --warmup "${WARMUP}" \
  --seed "${SEED}" \
  --rare-maf-cutoff "${RARE_MAF_CUTOFF}" \
  --adj-variants "${ADJ_VARIANTS}" \
  --atol "${ATOL}" \
  --rtol "${RTOL}" \
  --out-json "${FULL_JSON_OUT}" \
  --out-md "${FULL_MD_OUT}"

echo "[4/4] Done."
echo "Main report: ${FULL_MD_OUT}"
