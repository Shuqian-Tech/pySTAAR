#!/usr/bin/env python3
"""Measure cold-process versus warm-inprocess workflow latency."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import statistics
import subprocess
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


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _run_warm_once(scenario: Scenario) -> tuple[float, float]:
    fn = getattr(workflows, scenario.function_name)
    start = time.perf_counter()
    result = fn(**scenario.kwargs)
    elapsed = time.perf_counter() - start
    sentinel = float(result[scenario.sentinel_key])
    return elapsed, sentinel


def _run_cold_once(scenario: Scenario) -> tuple[float, float]:
    code = f"""
import json
import sys
from pathlib import Path
sys.path.insert(0, {str(SRC)!r})
from pystaar import workflows
fn = getattr(workflows, {scenario.function_name!r})
kwargs = json.loads({json.dumps(scenario.kwargs)!r})
result = fn(**kwargs)
print(float(result[{scenario.sentinel_key!r}]))
"""
    start = time.perf_counter()
    completed = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - start
    stdout = completed.stdout.strip().splitlines()
    if not stdout:
        raise RuntimeError(f"Cold-process run produced no output for {scenario.scenario_id}.")
    sentinel = float(stdout[-1])
    return elapsed, sentinel


def _write_report(
    path: Path,
    metadata: dict[str, Any],
    comparison_rows: list[dict[str, Any]],
) -> None:
    lines: list[str] = []
    lines.append("# Phase 3 Cold vs Warm Latency")
    lines.append("")
    lines.append(f"- Generated: {metadata['generated_utc']}")
    lines.append(f"- Dataset fingerprint: `{metadata['dataset_fingerprint_file']}`")
    lines.append(f"- Warm-up policy (warm mode): {metadata['warmup_runs']} run(s) discarded")
    lines.append(f"- Measured runs: {metadata['measured_runs']} per scenario per mode")
    lines.append("- Cold mode includes interpreter startup + imports + single workflow call.")
    lines.append("- Warm mode measures repeated calls in the same Python process.")
    lines.append("")
    lines.append("## Results")
    lines.append("")
    lines.append("| Scenario | Cold median (s) | Warm median (s) | Warm speedup vs cold |")
    lines.append("|---|---:|---:|---:|")
    for row in comparison_rows:
        lines.append(
            "| "
            + row["scenario_id"]
            + f" | {row['cold_median_seconds']:.6f}"
            + f" | {row['warm_median_seconds']:.6f}"
            + f" | {row['warm_speedup_vs_cold']:.2f}x |"
        )
    lines.append("")
    lines.append("## Artifacts")
    lines.append("")
    lines.append("- `benchmarks/phase3_cold_warm_raw.csv`")
    lines.append("- `benchmarks/phase3_cold_warm_summary.csv`")
    lines.append("- `benchmarks/phase3_cold_warm_comparison.csv`")
    lines.append("- `benchmarks/phase3_cold_warm_meta.json`")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--warmup-runs", type=int, default=1)
    parser.add_argument("--measured-runs", type=int, default=3)
    parser.add_argument(
        "--raw-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_cold_warm_raw.csv",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_cold_warm_summary.csv",
    )
    parser.add_argument(
        "--comparison-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_cold_warm_comparison.csv",
    )
    parser.add_argument(
        "--meta-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_cold_warm_meta.json",
    )
    parser.add_argument(
        "--report-output",
        type=Path,
        default=ROOT / "reports" / "performance_cold_warm.md",
    )
    args = parser.parse_args()

    if args.warmup_runs < 0:
        raise ValueError("warmup-runs must be >= 0")
    if args.measured_runs < 1:
        raise ValueError("measured-runs must be >= 1")

    raw_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    sentinel_checks: dict[str, float] = {}

    for scenario in SCENARIOS:
        for warmup_idx in range(args.warmup_runs):
            elapsed, sentinel = _run_warm_once(scenario)
            raw_rows.append(
                {
                    "mode": "warm_inprocess",
                    "scenario_id": scenario.scenario_id,
                    "run_type": "warmup",
                    "run_index": warmup_idx + 1,
                    "seconds": elapsed,
                    "sentinel_value": sentinel,
                }
            )
            sentinel_checks[f"warm.{scenario.scenario_id}.warmup{warmup_idx + 1}"] = sentinel

        warm_times: list[float] = []
        for run_idx in range(args.measured_runs):
            elapsed, sentinel = _run_warm_once(scenario)
            warm_times.append(elapsed)
            raw_rows.append(
                {
                    "mode": "warm_inprocess",
                    "scenario_id": scenario.scenario_id,
                    "run_type": "measured",
                    "run_index": run_idx + 1,
                    "seconds": elapsed,
                    "sentinel_value": sentinel,
                }
            )
            sentinel_checks[f"warm.{scenario.scenario_id}.run{run_idx + 1}"] = sentinel

        cold_times: list[float] = []
        for run_idx in range(args.measured_runs):
            elapsed, sentinel = _run_cold_once(scenario)
            cold_times.append(elapsed)
            raw_rows.append(
                {
                    "mode": "cold_process",
                    "scenario_id": scenario.scenario_id,
                    "run_type": "measured",
                    "run_index": run_idx + 1,
                    "seconds": elapsed,
                    "sentinel_value": sentinel,
                }
            )
            sentinel_checks[f"cold.{scenario.scenario_id}.run{run_idx + 1}"] = sentinel

        warm_median = statistics.median(warm_times)
        cold_median = statistics.median(cold_times)

        summary_rows.extend(
            [
                {
                    "mode": "warm_inprocess",
                    "scenario_id": scenario.scenario_id,
                    "warmup_runs": args.warmup_runs,
                    "measured_runs": args.measured_runs,
                    "median_seconds": warm_median,
                    "mean_seconds": statistics.fmean(warm_times),
                    "min_seconds": min(warm_times),
                    "max_seconds": max(warm_times),
                    "std_seconds": statistics.stdev(warm_times) if len(warm_times) > 1 else 0.0,
                },
                {
                    "mode": "cold_process",
                    "scenario_id": scenario.scenario_id,
                    "warmup_runs": 0,
                    "measured_runs": args.measured_runs,
                    "median_seconds": cold_median,
                    "mean_seconds": statistics.fmean(cold_times),
                    "min_seconds": min(cold_times),
                    "max_seconds": max(cold_times),
                    "std_seconds": statistics.stdev(cold_times) if len(cold_times) > 1 else 0.0,
                },
            ]
        )
        comparison_rows.append(
            {
                "scenario_id": scenario.scenario_id,
                "cold_median_seconds": cold_median,
                "warm_median_seconds": warm_median,
                "delta_seconds": warm_median - cold_median,
                "warm_speedup_vs_cold": cold_median / warm_median if warm_median > 0 else float("inf"),
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
        fieldnames=[
            "mode",
            "scenario_id",
            "run_type",
            "run_index",
            "seconds",
            "sentinel_value",
        ],
    )
    _write_csv(
        args.summary_output,
        summary_rows,
        fieldnames=[
            "mode",
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
    _write_csv(
        args.comparison_output,
        comparison_rows,
        fieldnames=[
            "scenario_id",
            "cold_median_seconds",
            "warm_median_seconds",
            "delta_seconds",
            "warm_speedup_vs_cold",
        ],
    )
    _write_json(args.meta_output, metadata)
    _write_report(args.report_output, metadata, comparison_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
