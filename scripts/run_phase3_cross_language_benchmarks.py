#!/usr/bin/env python3
"""Run Phase 3 cross-language benchmarks (Python vs R) and write comparison artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, cwd=ROOT, check=True)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_report(
    path: Path,
    rows: list[dict[str, Any]],
    warmup_runs: int,
    measured_runs: int,
    py_meta: dict[str, Any],
    r_meta: dict[str, Any],
) -> None:
    speedups = [row["python_vs_r_speedup"] for row in rows if row["python_vs_r_speedup"] > 0.0]
    geomean = statistics.geometric_mean(speedups) if speedups else float("nan")

    lines: list[str] = []
    lines.append("# Phase 3 Performance (Python vs R)")
    lines.append("")
    lines.append(f"- Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}")
    lines.append("- Dataset fingerprint: `baselines/example_fingerprint.json`")
    lines.append(f"- Warm-up policy: {warmup_runs} warm-up run(s) discarded")
    lines.append(f"- Measured runs: {measured_runs} per scenario")
    lines.append("- Reported statistic: median seconds")
    lines.append("- Reference backend details: `reports/reference_backend.md`")
    lines.append("- Python benchmark logs: `benchmarks/phase3_baseline_raw.csv`, `benchmarks/phase3_baseline_summary.csv`, `benchmarks/phase3_baseline_meta.json`")
    lines.append("- R benchmark logs: `benchmarks/phase3_baseline_r_raw.csv`, `benchmarks/phase3_baseline_r_summary.csv`, `benchmarks/phase3_baseline_r_meta.json`")
    lines.append("- Cross-language comparison: `benchmarks/phase3_cross_language_comparison.csv`")
    lines.append("")
    lines.append("## Environment")
    lines.append("")
    lines.append(f"- Python version: `{py_meta.get('platform', {}).get('python_version', 'unknown')}`")
    lines.append(f"- R version: `{r_meta.get('r_version', 'unknown')}`")
    lines.append(f"- STAAR package version: `{r_meta.get('staar_version', 'unknown')}`")
    lines.append("")
    lines.append("## Cross-Language Results")
    lines.append("")
    lines.append("| Scenario | Python median (s) | R median (s) | Python vs R speedup |")
    lines.append("|---|---:|---:|---:|")
    for row in rows:
        lines.append(
            f"| {row['scenario_id']}"
            f" | {row['python_median_seconds']:.6f}"
            f" | {row['r_median_seconds']:.6f}"
            f" | {row['python_vs_r_speedup']:.2f}x |"
        )
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(
        f"- Geometric mean Python speedup vs R across scenarios: `{geomean:.2f}x`."
        if math.isfinite(geomean)
        else "- Geometric mean Python speedup vs R: unavailable."
    )
    lines.append("- Speedup is computed as `R median / Python median` (values > 1 mean Python is faster).")
    lines.append("- This report is the Phase 3 baseline before Python optimizations.")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--warmup-runs", type=int, default=1)
    parser.add_argument("--measured-runs", type=int, default=5)
    parser.add_argument("--rscript-bin", default="Rscript")
    parser.add_argument(
        "--python-raw-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_baseline_raw.csv",
    )
    parser.add_argument(
        "--python-summary-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_baseline_summary.csv",
    )
    parser.add_argument(
        "--python-meta-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_baseline_meta.json",
    )
    parser.add_argument(
        "--r-raw-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_baseline_r_raw.csv",
    )
    parser.add_argument(
        "--r-summary-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_baseline_r_summary.csv",
    )
    parser.add_argument(
        "--r-meta-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_baseline_r_meta.json",
    )
    parser.add_argument(
        "--comparison-output",
        type=Path,
        default=ROOT / "benchmarks" / "phase3_cross_language_comparison.csv",
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

    _run(
        [
            sys.executable,
            str(ROOT / "scripts" / "run_phase3_benchmarks.py"),
            "--warmup-runs",
            str(args.warmup_runs),
            "--measured-runs",
            str(args.measured_runs),
            "--raw-output",
            str(args.python_raw_output),
            "--summary-output",
            str(args.python_summary_output),
            "--meta-output",
            str(args.python_meta_output),
            "--report-output",
            str(args.report_output),
        ]
    )

    _run(
        [
            args.rscript_bin,
            str(ROOT / "scripts" / "run_phase3_benchmarks_r.R"),
            "--warmup-runs",
            str(args.warmup_runs),
            "--measured-runs",
            str(args.measured_runs),
            "--raw-output",
            str(args.r_raw_output),
            "--summary-output",
            str(args.r_summary_output),
            "--meta-output",
            str(args.r_meta_output),
        ]
    )

    py_rows = _read_csv(args.python_summary_output)
    r_rows = _read_csv(args.r_summary_output)
    py_by_id = {row["scenario_id"]: row for row in py_rows}
    r_by_id = {row["scenario_id"]: row for row in r_rows}

    missing_in_r = sorted(set(py_by_id) - set(r_by_id))
    missing_in_py = sorted(set(r_by_id) - set(py_by_id))
    if missing_in_r or missing_in_py:
        raise ValueError(
            f"Scenario mismatch between Python and R summaries. "
            f"missing_in_r={missing_in_r}, missing_in_py={missing_in_py}"
        )

    comparison_rows: list[dict[str, Any]] = []
    for scenario_id in [row["scenario_id"] for row in py_rows]:
        py = py_by_id[scenario_id]
        r = r_by_id[scenario_id]
        py_median = float(py["median_seconds"])
        r_median = float(r["median_seconds"])
        py_vs_r = (r_median / py_median) if py_median > 0 else float("nan")
        r_vs_py = (py_median / r_median) if r_median > 0 else float("nan")
        comparison_rows.append(
            {
                "scenario_id": scenario_id,
                "warmup_runs": args.warmup_runs,
                "measured_runs": args.measured_runs,
                "python_median_seconds": py_median,
                "r_median_seconds": r_median,
                "python_min_seconds": float(py["min_seconds"]),
                "python_max_seconds": float(py["max_seconds"]),
                "r_min_seconds": float(r["min_seconds"]),
                "r_max_seconds": float(r["max_seconds"]),
                "python_vs_r_speedup": py_vs_r,
                "r_vs_python_speedup": r_vs_py,
            }
        )

    _write_csv(
        args.comparison_output,
        comparison_rows,
        fieldnames=[
            "scenario_id",
            "warmup_runs",
            "measured_runs",
            "python_median_seconds",
            "r_median_seconds",
            "python_min_seconds",
            "python_max_seconds",
            "r_min_seconds",
            "r_max_seconds",
            "python_vs_r_speedup",
            "r_vs_python_speedup",
        ],
    )

    py_meta = _read_json(args.python_meta_output)
    r_meta = _read_json(args.r_meta_output)
    _write_report(
        path=args.report_output,
        rows=comparison_rows,
        warmup_runs=args.warmup_runs,
        measured_runs=args.measured_runs,
        py_meta=py_meta,
        r_meta=r_meta,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
