#!/usr/bin/env python3
"""Run Phase 3 baseline benchmarks for representative pySTAAR workflows."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import statistics
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pystaar import workflows  # noqa: E402


@dataclass(frozen=True)
class Scenario:
    scenario_id: str
    function_name: str
    kwargs: dict[str, Any]
    sentinel_key: str


SCENARIOS: tuple[Scenario, ...] = (
    Scenario(
        scenario_id="staar_unrelated_glm",
        function_name="staar_unrelated_glm",
        kwargs={"dataset": "example", "seed": 600, "rare_maf_cutoff": 0.05},
        sentinel_key="results_STAAR_O",
    ),
    Scenario(
        scenario_id="staar_related_sparse_glmmkin_pure",
        function_name="staar_related_sparse_glmmkin",
        kwargs={
            "dataset": "example",
            "seed": 600,
            "rare_maf_cutoff": 0.05,
            "use_precomputed_artifacts": False,
        },
        sentinel_key="results_STAAR_O",
    ),
    Scenario(
        scenario_id="staar_unrelated_binary_spa",
        function_name="staar_unrelated_binary_spa",
        kwargs={"dataset": "example", "seed": 600, "rare_maf_cutoff": 0.05},
        sentinel_key="results_STAAR_B",
    ),
    Scenario(
        scenario_id="staar_related_sparse_binary_spa_pure",
        function_name="staar_related_sparse_binary_spa",
        kwargs={
            "dataset": "example",
            "seed": 600,
            "rare_maf_cutoff": 0.05,
            "use_precomputed_artifacts": False,
        },
        sentinel_key="results_STAAR_B",
    ),
    Scenario(
        scenario_id="staar_unrelated_glm_cond",
        function_name="staar_unrelated_glm_cond",
        kwargs={"dataset": "example", "seed": 600, "rare_maf_cutoff": 0.05},
        sentinel_key="results_STAAR_O_cond",
    ),
    Scenario(
        scenario_id="indiv_score_unrelated_glm",
        function_name="indiv_score_unrelated_glm",
        kwargs={"dataset": "example", "seed": 600, "rare_maf_cutoff": 0.05},
        sentinel_key="pvalue_min",
    ),
    Scenario(
        scenario_id="ai_staar_unrelated_glm",
        function_name="ai_staar_unrelated_glm",
        kwargs={"dataset": "example", "seed": 600, "rare_maf_cutoff": 0.05},
        sentinel_key="results_STAAR_O",
    ),
    Scenario(
        scenario_id="ai_staar_related_sparse_glmmkin_find_weight_pure",
        function_name="ai_staar_related_sparse_glmmkin_find_weight",
        kwargs={
            "dataset": "example",
            "seed": 600,
            "rare_maf_cutoff": 0.05,
            "use_precomputed_artifacts": False,
        },
        sentinel_key="results_STAAR_O",
    ),
)


def _run_once(scenario: Scenario) -> tuple[float, float]:
    fn = getattr(workflows, scenario.function_name)
    start = time.perf_counter()
    result = fn(**scenario.kwargs)
    elapsed = time.perf_counter() - start
    if scenario.sentinel_key not in result:
        raise KeyError(
            f"{scenario.scenario_id} is missing expected key '{scenario.sentinel_key}'."
        )
    sentinel_value = float(result[scenario.sentinel_key])
    return elapsed, sentinel_value


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _write_report(path: Path, summary_rows: list[dict[str, Any]], metadata: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Phase 3 Performance Baseline")
    lines.append("")
    lines.append(f"- Generated: {metadata['generated_utc']}")
    lines.append(f"- Dataset fingerprint: `{metadata['dataset_fingerprint_file']}`")
    lines.append(f"- Warm-up policy: {metadata['warmup_runs']} warm-up run(s) discarded")
    lines.append(f"- Measured runs: {metadata['measured_runs']} per scenario")
    lines.append("- Reported statistic: median seconds")
    lines.append("- Reference backend details: `reports/reference_backend.md`")
    lines.append("- Python environment details: `reports/python_environment.md`")
    lines.append("")
    lines.append("## Baseline Results")
    lines.append("")
    lines.append("| Scenario | Median (s) | Min (s) | Max (s) | Speedup vs baseline |")
    lines.append("|---|---:|---:|---:|---:|")
    for row in summary_rows:
        lines.append(
            "| "
            + row["scenario_id"]
            + f" | {row['median_seconds']:.6f}"
            + f" | {row['min_seconds']:.6f}"
            + f" | {row['max_seconds']:.6f}"
            + f" | {row['relative_speedup_vs_baseline']:.2f}x |"
        )
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("- This report captures baseline Phase 3 performance before optimization.")
    lines.append("- Relative speedups are 1.00x by definition until optimized variants are added.")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--warmup-runs", type=int, default=1)
    parser.add_argument("--measured-runs", type=int, default=5)
    parser.add_argument(
        "--raw-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_baseline_raw.csv",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_baseline_summary.csv",
    )
    parser.add_argument(
        "--meta-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_baseline_meta.json",
    )
    parser.add_argument(
        "--report-output",
        type=Path,
        default=ROOT / "reports" / "performance.md",
    )
    args = parser.parse_args()

    if args.warmup_runs < 0:
        raise ValueError("warmup-runs must be >= 0")
    if args.measured_runs < 1:
        raise ValueError("measured-runs must be >= 1")

    raw_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    sentinel_checks: dict[str, float] = {}

    for scenario in SCENARIOS:
        for warmup_idx in range(args.warmup_runs):
            elapsed, sentinel_value = _run_once(scenario)
            raw_rows.append(
                {
                    "scenario_id": scenario.scenario_id,
                    "run_type": "warmup",
                    "run_index": warmup_idx + 1,
                    "seconds": elapsed,
                    "sentinel_value": sentinel_value,
                }
            )
            sentinel_checks[f"{scenario.scenario_id}.warmup"] = sentinel_value

        measured: list[float] = []
        for run_idx in range(args.measured_runs):
            elapsed, sentinel_value = _run_once(scenario)
            measured.append(elapsed)
            raw_rows.append(
                {
                    "scenario_id": scenario.scenario_id,
                    "run_type": "measured",
                    "run_index": run_idx + 1,
                    "seconds": elapsed,
                    "sentinel_value": sentinel_value,
                }
            )
            sentinel_checks[f"{scenario.scenario_id}.run{run_idx + 1}"] = sentinel_value

        summary_rows.append(
            {
                "scenario_id": scenario.scenario_id,
                "warmup_runs": args.warmup_runs,
                "measured_runs": args.measured_runs,
                "median_seconds": statistics.median(measured),
                "mean_seconds": statistics.fmean(measured),
                "min_seconds": min(measured),
                "max_seconds": max(measured),
                "std_seconds": statistics.stdev(measured) if len(measured) > 1 else 0.0,
                "relative_speedup_vs_baseline": 1.0,
            }
        )

    generated_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    metadata = {
        "generated_utc": generated_utc,
        "warmup_runs": args.warmup_runs,
        "measured_runs": args.measured_runs,
        "dataset_fingerprint_file": "baselines/example_fingerprint.json",
        "reference_backend_file": "reports/reference_backend.md",
        "python_environment_file": "reports/python_environment.md",
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "python_version": platform.python_version(),
        },
        "scenarios": [
            {
                "scenario_id": s.scenario_id,
                "function_name": s.function_name,
                "kwargs": s.kwargs,
                "sentinel_key": s.sentinel_key,
            }
            for s in SCENARIOS
        ],
        "sentinel_checks": sentinel_checks,
    }

    _write_csv(
        args.raw_output,
        raw_rows,
        fieldnames=["scenario_id", "run_type", "run_index", "seconds", "sentinel_value"],
    )
    _write_csv(
        args.summary_output,
        summary_rows,
        fieldnames=[
            "scenario_id",
            "warmup_runs",
            "measured_runs",
            "median_seconds",
            "mean_seconds",
            "min_seconds",
            "max_seconds",
            "std_seconds",
            "relative_speedup_vs_baseline",
        ],
    )
    _write_json(args.meta_output, metadata)
    _write_report(args.report_output, summary_rows, metadata)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
