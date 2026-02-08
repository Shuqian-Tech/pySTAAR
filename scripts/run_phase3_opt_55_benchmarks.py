#!/usr/bin/env python3
"""Benchmark STAAR-55 optimization for related GLMM workflows."""

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

MODES = ("baseline_no_cache", "optimized_cache")


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _configure_mode(mode: str, original_fit_null_glmmkin):
    if mode == "baseline_no_cache":
        # workflows._fit_related_nullmodel bypasses caching when fit_null_glmmkin identity changes.
        def _passthrough_fit_null_glmmkin(*args, **kwargs):
            return original_fit_null_glmmkin(*args, **kwargs)

        workflows.fit_null_glmmkin = _passthrough_fit_null_glmmkin
        return
    if mode == "optimized_cache":
        workflows.fit_null_glmmkin = original_fit_null_glmmkin
        return
    raise ValueError(f"Unsupported mode: {mode}")


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


def _build_report(
    path: Path,
    metadata: dict[str, Any],
    summary_rows: list[dict[str, Any]],
    comparison_rows: list[dict[str, Any]],
) -> None:
    baseline_map = {
        row["scenario_id"]: row for row in summary_rows if row["mode"] == "baseline_no_cache"
    }
    optimized_map = {
        row["scenario_id"]: row for row in summary_rows if row["mode"] == "optimized_cache"
    }

    lines: list[str] = []
    lines.append("# Phase 3 Performance Snapshot (STAAR-55)")
    lines.append("")
    lines.append(f"- Generated: {metadata['generated_utc']}")
    lines.append(f"- Dataset fingerprint: `{metadata['dataset_fingerprint_file']}`")
    lines.append(f"- Warm-up policy: {metadata['warmup_runs']} warm-up run(s) discarded")
    lines.append(f"- Measured runs: {metadata['measured_runs']} per scenario per mode")
    lines.append("- Reported statistic: median seconds")
    lines.append("- Reference backend details: `reports/reference_backend.md`")
    lines.append("- Python environment details: `reports/python_environment.md`")
    lines.append("")
    lines.append("## Results by Mode")
    lines.append("")
    lines.append("| Scenario | Baseline no-cache median (s) | Optimized cache median (s) | Speedup |")
    lines.append("|---|---:|---:|---:|")
    for scenario in SCENARIOS:
        baseline = baseline_map[scenario.scenario_id]
        optimized = optimized_map[scenario.scenario_id]
        speedup = baseline["median_seconds"] / optimized["median_seconds"]
        lines.append(
            "| "
            + scenario.scenario_id
            + f" | {baseline['median_seconds']:.6f}"
            + f" | {optimized['median_seconds']:.6f}"
            + f" | {speedup:.2f}x |"
        )
    lines.append("")
    lines.append("## Comparison Artifact")
    lines.append("")
    lines.append("- `benchmarks/phase3_opt_55_comparison.csv`")
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("- Baseline mode disables related-null-model caching by replacing `workflows.fit_null_glmmkin` with a passthrough wrapper.")
    lines.append("- Optimized mode restores the original `fit_null_glmmkin` and uses the new cached related-null-model path.")
    lines.append("- Sentinels are expected to match between modes; values are recorded in `benchmarks/phase3_opt_55_meta.json`.")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--warmup-runs", type=int, default=1)
    parser.add_argument("--measured-runs", type=int, default=5)
    parser.add_argument(
        "--raw-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_opt_55_raw.csv",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_opt_55_summary.csv",
    )
    parser.add_argument(
        "--comparison-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_opt_55_comparison.csv",
    )
    parser.add_argument(
        "--meta-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_opt_55_meta.json",
    )
    parser.add_argument(
        "--report-output",
        type=Path,
        default=ROOT / "reports" / "performance_opt_55.md",
    )
    args = parser.parse_args()

    if args.warmup_runs < 0:
        raise ValueError("warmup-runs must be >= 0")
    if args.measured_runs < 1:
        raise ValueError("measured-runs must be >= 1")

    original_fit_null_glmmkin = workflows.fit_null_glmmkin
    raw_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    sentinel_checks: dict[str, float] = {}
    mode_scenario_sentinel: dict[tuple[str, str], float] = {}

    try:
        for scenario in SCENARIOS:
            for mode in MODES:
                _configure_mode(mode=mode, original_fit_null_glmmkin=original_fit_null_glmmkin)
                workflows._fit_related_nullmodel_cached.cache_clear()

                for warmup_idx in range(args.warmup_runs):
                    elapsed, sentinel_value = _run_once(scenario)
                    raw_rows.append(
                        {
                            "mode": mode,
                            "scenario_id": scenario.scenario_id,
                            "run_type": "warmup",
                            "run_index": warmup_idx + 1,
                            "seconds": elapsed,
                            "sentinel_value": sentinel_value,
                        }
                    )
                    sentinel_checks[
                        f"{mode}.{scenario.scenario_id}.warmup{warmup_idx + 1}"
                    ] = sentinel_value
                    mode_scenario_sentinel.setdefault((mode, scenario.scenario_id), sentinel_value)

                measured: list[float] = []
                for run_idx in range(args.measured_runs):
                    elapsed, sentinel_value = _run_once(scenario)
                    measured.append(elapsed)
                    raw_rows.append(
                        {
                            "mode": mode,
                            "scenario_id": scenario.scenario_id,
                            "run_type": "measured",
                            "run_index": run_idx + 1,
                            "seconds": elapsed,
                            "sentinel_value": sentinel_value,
                        }
                    )
                    sentinel_checks[
                        f"{mode}.{scenario.scenario_id}.run{run_idx + 1}"
                    ] = sentinel_value
                    mode_scenario_sentinel.setdefault((mode, scenario.scenario_id), sentinel_value)

                summary_rows.append(
                    {
                        "mode": mode,
                        "scenario_id": scenario.scenario_id,
                        "warmup_runs": args.warmup_runs,
                        "measured_runs": args.measured_runs,
                        "median_seconds": statistics.median(measured),
                        "mean_seconds": statistics.fmean(measured),
                        "min_seconds": min(measured),
                        "max_seconds": max(measured),
                        "std_seconds": statistics.stdev(measured) if len(measured) > 1 else 0.0,
                    }
                )
    finally:
        workflows.fit_null_glmmkin = original_fit_null_glmmkin

    summary_lookup = {
        (row["mode"], row["scenario_id"]): row for row in summary_rows
    }
    for scenario in SCENARIOS:
        baseline_sentinel = mode_scenario_sentinel[("baseline_no_cache", scenario.scenario_id)]
        optimized_sentinel = mode_scenario_sentinel[("optimized_cache", scenario.scenario_id)]
        if not abs(baseline_sentinel - optimized_sentinel) <= 1e-12:
            raise ValueError(
                "Sentinel mismatch between baseline and optimized mode for "
                f"{scenario.scenario_id}: {baseline_sentinel} vs {optimized_sentinel}"
            )
        baseline = summary_lookup[("baseline_no_cache", scenario.scenario_id)]
        optimized = summary_lookup[("optimized_cache", scenario.scenario_id)]
        comparison_rows.append(
            {
                "scenario_id": scenario.scenario_id,
                "baseline_no_cache_median_seconds": baseline["median_seconds"],
                "optimized_cache_median_seconds": optimized["median_seconds"],
                "delta_seconds": optimized["median_seconds"] - baseline["median_seconds"],
                "speedup_vs_baseline_no_cache": baseline["median_seconds"]
                / optimized["median_seconds"],
            }
        )

    generated_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    metadata = {
        "generated_utc": generated_utc,
        "warmup_runs": args.warmup_runs,
        "measured_runs": args.measured_runs,
        "modes": list(MODES),
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
            "baseline_no_cache_median_seconds",
            "optimized_cache_median_seconds",
            "delta_seconds",
            "speedup_vs_baseline_no_cache",
        ],
    )
    _write_json(args.meta_output, metadata)
    _build_report(args.report_output, metadata, summary_rows, comparison_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
