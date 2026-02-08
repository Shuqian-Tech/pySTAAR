#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import statistics
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from pystaar import workflows


WORKFLOW_ORDER = [
    "staar_unrelated_glm",
    "staar_related_sparse_glmmkin",
    "staar_related_dense_glmmkin",
    "staar_unrelated_glm_cond",
    "staar_related_sparse_glmmkin_cond",
    "staar_related_dense_glmmkin_cond",
    "staar_unrelated_binary_spa",
    "staar_related_sparse_binary_spa",
    "staar_related_dense_binary_spa",
]


def _to_builtin(value: Any) -> Any:
    try:
        import numpy as np
    except Exception:  # pragma: no cover - optional defensive path
        np = None

    if isinstance(value, dict):
        return {k: _to_builtin(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_to_builtin(v) for v in value]
    if isinstance(value, tuple):
        return [_to_builtin(v) for v in value]
    if np is not None:
        if isinstance(value, np.generic):
            return value.item()
    return value


def _extract_metrics(name: str, payload: dict[str, Any]) -> dict[str, float]:
    if name in {
        "staar_unrelated_glm",
        "staar_related_sparse_glmmkin",
        "staar_related_dense_glmmkin",
    }:
        return {
            "num_variant": float(payload["num_variant"]),
            "results_STAAR_O": float(payload["results_STAAR_O"]),
        }

    if name in {
        "staar_unrelated_glm_cond",
        "staar_related_sparse_glmmkin_cond",
        "staar_related_dense_glmmkin_cond",
    }:
        return {
            "num_variant": float(payload["num_variant"]),
            "cMAC": float(payload["cMAC"]),
            "results_STAAR_O_cond": float(payload["results_STAAR_O_cond"]),
            "results_ACAT_O_cond": float(payload["results_ACAT_O_cond"]),
        }

    if name in {
        "staar_unrelated_binary_spa",
        "staar_related_sparse_binary_spa",
        "staar_related_dense_binary_spa",
    }:
        return {
            "num_variant": float(payload["num_variant"]),
            "cMAC": float(payload["cMAC"]),
            "case_count": float(payload["case_count"]),
            "results_STAAR_B": float(payload["results_STAAR_B"]),
        }

    raise ValueError(f"Unknown workflow for metric extraction: {name}")


def _run_once(
    dataset_dir: str,
    seed: int,
    rare_maf_cutoff: float,
    adj_variant_indices: tuple[int, ...],
) -> tuple[dict[str, dict[str, float]], dict[str, float], float]:
    results: dict[str, dict[str, float]] = {}
    timings: dict[str, float] = {}

    suite_start = time.perf_counter()

    for name in WORKFLOW_ORDER:
        fn = getattr(workflows, name)
        t0 = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            if name.endswith("_cond"):
                payload = fn(
                    dataset=dataset_dir,
                    seed=seed,
                    rare_maf_cutoff=rare_maf_cutoff,
                    method_cond="optimal",
                    adj_variant_indices=adj_variant_indices,
                )
            else:
                payload = fn(
                    dataset=dataset_dir,
                    seed=seed,
                    rare_maf_cutoff=rare_maf_cutoff,
                )
        timings[name] = float(time.perf_counter() - t0)
        results[name] = _extract_metrics(name, payload)

    total_elapsed = float(time.perf_counter() - suite_start)
    return results, timings, total_elapsed


def _aggregate_timings(
    timing_rows: list[dict[str, float]],
    total_rows: list[float],
    runs: int,
    warmup: int,
) -> dict[str, Any]:
    workflows_summary: dict[str, Any] = {}
    for name in WORKFLOW_ORDER:
        values = [float(row[name]) for row in timing_rows]
        workflows_summary[name] = {
            "timings_seconds": values,
            "median_seconds": float(statistics.median(values)),
            "mean_seconds": float(statistics.mean(values)),
            "min_seconds": float(min(values)),
            "max_seconds": float(max(values)),
        }

    return {
        "benchmark": {
            "language": "Python",
            "suite": "full_workflow_simulated",
            "runs": int(runs),
            "warmup": int(warmup),
            "total_timings_seconds": [float(v) for v in total_rows],
            "total_median_seconds": float(statistics.median(total_rows)),
            "total_mean_seconds": float(statistics.mean(total_rows)),
            "workflows": workflows_summary,
        }
    }


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


@dataclass
class Check:
    workflow: str
    metric: str
    r_value: float
    py_value: float
    abs_diff: float
    rel_diff: float
    atol: float
    rtol: float
    passed: bool


def _compare_results(
    r_results: dict[str, Any],
    py_results: dict[str, Any],
    atol: float,
    rtol: float,
) -> dict[str, Any]:
    checks: list[Check] = []
    failed: list[str] = []

    for workflow in WORKFLOW_ORDER:
        if workflow not in r_results:
            failed.append(f"{workflow}:missing_in_r")
            continue
        if workflow not in py_results:
            failed.append(f"{workflow}:missing_in_python")
            continue

        r_metrics = r_results[workflow]
        py_metrics = py_results[workflow]

        keys = sorted(set(r_metrics.keys()) | set(py_metrics.keys()))
        for key in keys:
            if key not in r_metrics or key not in py_metrics:
                failed.append(f"{workflow}.{key}:missing_key")
                continue

            r_value = float(r_metrics[key])
            py_value = float(py_metrics[key])
            abs_diff = abs(py_value - r_value)
            rel_diff = abs_diff / max(abs(r_value), 1e-12)
            passed = abs_diff <= (atol + rtol * abs(r_value))
            checks.append(
                Check(
                    workflow=workflow,
                    metric=key,
                    r_value=r_value,
                    py_value=py_value,
                    abs_diff=abs_diff,
                    rel_diff=rel_diff,
                    atol=atol,
                    rtol=rtol,
                    passed=passed,
                )
            )
            if not passed:
                failed.append(f"{workflow}.{key}")

    return {
        "passed": len(failed) == 0,
        "atol": float(atol),
        "rtol": float(rtol),
        "failed_checks": failed,
        "checks": [
            {
                "workflow": c.workflow,
                "metric": c.metric,
                "r_value": c.r_value,
                "py_value": c.py_value,
                "abs_diff": c.abs_diff,
                "rel_diff": c.rel_diff,
                "atol": c.atol,
                "rtol": c.rtol,
                "passed": c.passed,
            }
            for c in checks
        ],
    }


def _compare_performance(
    r_benchmark: dict[str, Any],
    py_benchmark: dict[str, Any],
) -> dict[str, Any]:
    r_data = r_benchmark["benchmark"]
    py_data = py_benchmark["benchmark"]

    workflows: dict[str, Any] = {}
    for name in WORKFLOW_ORDER:
        r_median = float(r_data["workflows"][name]["median_seconds"])
        py_median = float(py_data["workflows"][name]["median_seconds"])
        speedup = None if py_median <= 0 else r_median / py_median
        workflows[name] = {
            "r_median_seconds": r_median,
            "py_median_seconds": py_median,
            "speedup_r_div_py": speedup,
        }

    total_r = float(r_data["total_median_seconds"])
    total_py = float(py_data["total_median_seconds"])
    total_speedup = None if total_py <= 0 else total_r / total_py

    return {
        "suite": "full_workflow_simulated",
        "interpretation": "speedup_r_div_py > 1 means Python is faster on this benchmark target.",
        "total": {
            "r_median_seconds": total_r,
            "py_median_seconds": total_py,
            "speedup_r_div_py": total_speedup,
        },
        "workflows": workflows,
    }


def _render_markdown(report: dict[str, Any]) -> str:
    parity = report["parity"]
    performance = report["performance"]

    lines: list[str] = []
    lines.append("# 1KG Simulated Full-Workflow Report")
    lines.append("")
    lines.append("## Scope")
    lines.append("- Based on 1KG-derived genotype; `phred/pheno/kins` are simulated reproducibly.")
    lines.append("- Workflow families: STAAR, STAAR_cond, STAAR_Binary_SPA.")
    lines.append("- Coverage includes unrelated + related (sparse/dense) routes.")
    lines.append("")
    lines.append("## Correctness")
    lines.append(f"- Overall: {'PASS' if parity['passed'] else 'FAIL'}")
    lines.append(f"- Tolerance: atol={parity['atol']}, rtol={parity['rtol']}")
    if parity["failed_checks"]:
        lines.append(f"- Failed checks: {', '.join(parity['failed_checks'])}")
    else:
        lines.append("- Failed checks: none")
    lines.append("")
    lines.append("| Workflow | Metric | Status | Abs Diff | Rel Diff |")
    lines.append("|---|---|---:|---:|---:|")
    for item in parity["checks"]:
        status = "PASS" if item["passed"] else "FAIL"
        lines.append(
            f"| {item['workflow']} | {item['metric']} | {status} | "
            f"{item['abs_diff']:.3e} | {item['rel_diff']:.3e} |"
        )

    lines.append("")
    lines.append("## Performance")
    total = performance["total"]
    lines.append(f"- Total median (R): {total['r_median_seconds']:.4f}s")
    lines.append(f"- Total median (Python): {total['py_median_seconds']:.4f}s")
    if total["speedup_r_div_py"] is None:
        lines.append("- Total speedup (R/Python): unavailable")
    else:
        lines.append(f"- Total speedup (R/Python): {total['speedup_r_div_py']:.3f}x")
    lines.append(f"- Note: {performance['interpretation']}")
    lines.append("")
    lines.append("| Workflow | R median (s) | Python median (s) | R/Python |")
    lines.append("|---|---:|---:|---:|")
    for name in WORKFLOW_ORDER:
        item = performance["workflows"][name]
        speedup = item["speedup_r_div_py"]
        speedup_str = "-" if speedup is None else f"{speedup:.3f}x"
        lines.append(
            f"| {name} | {item['r_median_seconds']:.4f} | "
            f"{item['py_median_seconds']:.4f} | {speedup_str} |"
        )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run Python full workflow suite and compare to R outputs."
    )
    parser.add_argument("--dataset-dir", required=True)
    parser.add_argument("--r-results", required=True)
    parser.add_argument("--r-benchmark", required=True)
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--seed", type=int, default=600)
    parser.add_argument("--rare-maf-cutoff", type=float, default=0.05)
    parser.add_argument("--adj-variants", default="1", help="1-based indices, comma-separated")
    parser.add_argument("--atol", type=float, default=1e-6)
    parser.add_argument("--rtol", type=float, default=1e-3)
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()

    if args.runs < 1:
        raise ValueError("--runs must be >= 1")
    if args.warmup < 0:
        raise ValueError("--warmup must be >= 0")

    dataset_dir = str(Path(args.dataset_dir).expanduser().resolve())
    r_results_path = Path(args.r_results).expanduser().resolve()
    r_benchmark_path = Path(args.r_benchmark).expanduser().resolve()
    out_json = Path(args.out_json).expanduser().resolve()
    out_md = Path(args.out_md).expanduser().resolve()
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)

    r_results = _load_json(r_results_path)
    r_benchmark = _load_json(r_benchmark_path)

    benchmark_adj = r_benchmark.get("benchmark", {}).get("adj_variants_used")
    if isinstance(benchmark_adj, list) and len(benchmark_adj) > 0:
        adj_1_based = tuple(int(v) for v in benchmark_adj)
    else:
        adj_1_based = tuple(int(v.strip()) for v in args.adj_variants.split(",") if v.strip())
    if not adj_1_based:
        raise ValueError("--adj-variants must contain at least one index")
    adj_0_based = tuple(v - 1 for v in adj_1_based)
    if min(adj_0_based) < 0:
        raise ValueError("--adj-variants must be >= 1")

    for _ in range(args.warmup):
        _run_once(
            dataset_dir=dataset_dir,
            seed=args.seed,
            rare_maf_cutoff=args.rare_maf_cutoff,
            adj_variant_indices=adj_0_based,
        )

    timing_rows: list[dict[str, float]] = []
    total_rows: list[float] = []
    py_results: dict[str, dict[str, float]] | None = None
    for _ in range(args.runs):
        run_results, run_timing, total_elapsed = _run_once(
            dataset_dir=dataset_dir,
            seed=args.seed,
            rare_maf_cutoff=args.rare_maf_cutoff,
            adj_variant_indices=adj_0_based,
        )
        py_results = run_results
        timing_rows.append(run_timing)
        total_rows.append(total_elapsed)

    assert py_results is not None
    py_benchmark = _aggregate_timings(
        timing_rows=timing_rows,
        total_rows=total_rows,
        runs=args.runs,
        warmup=args.warmup,
    )

    parity = _compare_results(
        r_results=r_results,
        py_results=py_results,
        atol=args.atol,
        rtol=args.rtol,
    )
    performance = _compare_performance(
        r_benchmark=r_benchmark,
        py_benchmark=py_benchmark,
    )

    report = {
        "meta": {
            "scope": "1kg_simulated_full_workflow",
            "dataset_dir": dataset_dir,
            "seed": int(args.seed),
            "rare_maf_cutoff": float(args.rare_maf_cutoff),
            "adj_variants_1_based": list(adj_1_based),
            "r_results_file": str(r_results_path),
            "r_benchmark_file": str(r_benchmark_path),
        },
        "parity": parity,
        "performance": performance,
        "r_results": r_results,
        "python_results": py_results,
        "r_benchmark": r_benchmark,
        "python_benchmark": py_benchmark,
    }

    out_json.write_text(json.dumps(_to_builtin(report), indent=2, ensure_ascii=False), encoding="utf-8")
    out_md.write_text(_render_markdown(report), encoding="utf-8")
    print(f"Report written: {out_json}")
    print(f"Report written: {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
