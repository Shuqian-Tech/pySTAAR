#!/usr/bin/env python3
"""Run Phase 3 non-example dataset benchmark probe."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import shutil
import statistics
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
DATA_DIR = ROOT / "data"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pystaar import workflows  # noqa: E402


EXAMPLE_DIRECTORY_FILE_MAP = {
    "example_geno.mtx": "geno.mtx",
    "example_phred.csv": "phred.csv",
    "example_pheno_unrelated.csv": "pheno_unrelated.csv",
    "example_pheno_related.csv": "pheno_related.csv",
    "example_kins_sparse.mtx": "kins_sparse.mtx",
    "example_kins_dense.mtx": "kins_dense.mtx",
}

SCENARIO_ID = "staar_unrelated_glm_nonexample_dir"


def _materialize_example_directory_copy(target_dir: Path) -> Path:
    target_dir.mkdir(parents=True, exist_ok=True)
    for src_name, dst_name in EXAMPLE_DIRECTORY_FILE_MAP.items():
        shutil.copy2(DATA_DIR / src_name, target_dir / dst_name)
    return target_dir


def _run_once(dataset_dir: Path) -> tuple[float, float]:
    start = time.perf_counter()
    result = workflows.staar_unrelated_glm(
        dataset=str(dataset_dir),
        seed=600,
        rare_maf_cutoff=0.05,
    )
    elapsed = time.perf_counter() - start
    sentinel = float(result["results_STAAR_O"])
    return elapsed, sentinel


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _write_report(path: Path, summary: dict[str, Any], metadata: dict[str, Any]) -> None:
    lines = [
        "# Phase 3 Non-Example Probe",
        "",
        f"- Generated: {metadata['generated_utc']}",
        "- Scenario: `staar_unrelated_glm` on runtime non-example directory dataset copy",
        "- Dataset source: runtime copy of `data/example_*` files into directory layout (`geno.mtx`, etc.)",
        f"- Warm-up policy: {metadata['warmup_runs']} warm-up run(s) discarded",
        f"- Measured runs: {metadata['measured_runs']}",
        "- Reported statistic: median seconds",
        "- Reference backend details: `reports/reference_backend.md`",
        "- Python environment details: `reports/python_environment.md`",
        "",
        "## Result",
        "",
        "| Scenario | Median (s) | Min (s) | Max (s) |",
        "|---|---:|---:|---:|",
        (
            f"| {summary['scenario_id']} | {summary['median_seconds']:.6f} "
            f"| {summary['min_seconds']:.6f} | {summary['max_seconds']:.6f} |"
        ),
        "",
        "## Notes",
        "",
        "- This probe is used for STAAR-45 non-example performance coverage tracking.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--warmup-runs", type=int, default=1)
    parser.add_argument("--measured-runs", type=int, default=5)
    parser.add_argument(
        "--raw-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_nonexample_raw.csv",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_nonexample_summary.csv",
    )
    parser.add_argument(
        "--meta-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_nonexample_meta.json",
    )
    parser.add_argument(
        "--report-output",
        type=Path,
        default=ROOT / "reports" / "performance_nonexample.md",
    )
    args = parser.parse_args()

    if args.warmup_runs < 0:
        raise ValueError("warmup-runs must be >= 0")
    if args.measured_runs < 1:
        raise ValueError("measured-runs must be >= 1")

    raw_rows: list[dict[str, Any]] = []
    measured: list[float] = []
    sentinel_checks: dict[str, float] = {}

    with tempfile.TemporaryDirectory(prefix="pystaar_nonexample_") as tmp:
        dataset_dir = _materialize_example_directory_copy(Path(tmp) / "dataset")

        for warmup_idx in range(args.warmup_runs):
            elapsed, sentinel = _run_once(dataset_dir)
            raw_rows.append(
                {
                    "scenario_id": SCENARIO_ID,
                    "run_type": "warmup",
                    "run_index": warmup_idx + 1,
                    "seconds": elapsed,
                    "sentinel_value": sentinel,
                }
            )
            sentinel_checks[f"warmup{warmup_idx + 1}"] = sentinel

        for run_idx in range(args.measured_runs):
            elapsed, sentinel = _run_once(dataset_dir)
            measured.append(elapsed)
            raw_rows.append(
                {
                    "scenario_id": SCENARIO_ID,
                    "run_type": "measured",
                    "run_index": run_idx + 1,
                    "seconds": elapsed,
                    "sentinel_value": sentinel,
                }
            )
            sentinel_checks[f"run{run_idx + 1}"] = sentinel

    summary_row = {
        "scenario_id": SCENARIO_ID,
        "warmup_runs": args.warmup_runs,
        "measured_runs": args.measured_runs,
        "median_seconds": statistics.median(measured),
        "mean_seconds": statistics.fmean(measured),
        "min_seconds": min(measured),
        "max_seconds": max(measured),
        "std_seconds": statistics.stdev(measured) if len(measured) > 1 else 0.0,
    }

    generated_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    metadata = {
        "generated_utc": generated_utc,
        "scenario_id": SCENARIO_ID,
        "dataset_mode": "runtime_directory_copy_of_example",
        "dataset_fingerprint_file": "baselines/example_fingerprint.json",
        "warmup_runs": args.warmup_runs,
        "measured_runs": args.measured_runs,
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "python_version": platform.python_version(),
        },
        "sentinel_checks": sentinel_checks,
    }

    _write_csv(
        args.raw_output,
        raw_rows,
        fieldnames=["scenario_id", "run_type", "run_index", "seconds", "sentinel_value"],
    )
    _write_csv(
        args.summary_output,
        [summary_row],
        fieldnames=[
            "scenario_id",
            "warmup_runs",
            "measured_runs",
            "median_seconds",
            "mean_seconds",
            "min_seconds",
            "max_seconds",
            "std_seconds",
        ],
    )
    _write_json(args.meta_output, metadata)
    _write_report(args.report_output, summary_row, metadata)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
