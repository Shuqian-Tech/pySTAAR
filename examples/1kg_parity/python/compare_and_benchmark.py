#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import json
import statistics
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread


def _to_builtin(value: Any) -> Any:
    if isinstance(value, dict):
        return {k: _to_builtin(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_to_builtin(v) for v in value]
    if isinstance(value, tuple):
        return [_to_builtin(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return [_to_builtin(v) for v in value.tolist()]
    if isinstance(value, (pd.Series,)):
        return [_to_builtin(v) for v in value.tolist()]
    return value


def _load_matrix_market(path: Path) -> sp.csc_matrix:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as f:
            matrix = mmread(f)
    else:
        matrix = mmread(path)

    if not sp.issparse(matrix):
        return sp.csc_matrix(matrix)
    return matrix.tocsc()


def _load_and_summarize(matrix_path: Path, maf_path: Path, snploc_path: Path) -> dict[str, Any]:
    geno = _load_matrix_market(matrix_path)
    maf = pd.read_csv(maf_path)
    snploc = pd.read_csv(snploc_path)

    maf_values = maf["maf"].to_numpy(dtype=float)
    maf_q = np.quantile(maf_values, [0.0, 0.01, 0.5, 0.99, 1.0])

    variant_id_match = False
    if "variant_id" in maf.columns and "variant_id" in snploc.columns:
        variant_id_match = (
            maf["variant_id"].astype(str).reset_index(drop=True).equals(
                snploc["variant_id"].astype(str).reset_index(drop=True)
            )
        )

    return {
        "n_samples": int(geno.shape[0]),
        "n_variants": int(geno.shape[1]),
        "nnz": int(geno.nnz),
        "genotype_sum": float(np.asarray(geno.data, dtype=float).sum()),
        "maf_mean": float(np.mean(maf_values)),
        "maf_summary": {
            "min": float(maf_q[0]),
            "q01": float(maf_q[1]),
            "median": float(maf_q[2]),
            "q99": float(maf_q[3]),
            "max": float(maf_q[4]),
        },
        "maf_rows": int(maf.shape[0]),
        "snploc_rows": int(snploc.shape[0]),
        "variant_id_match": bool(variant_id_match),
    }


def _run_benchmark(
    matrix_path: Path,
    maf_path: Path,
    snploc_path: Path,
    runs: int,
    warmup: int,
) -> dict[str, Any]:
    if warmup < 0:
        raise ValueError("warmup must be >= 0")
    if runs < 1:
        raise ValueError("runs must be >= 1")

    for _ in range(warmup):
        _load_and_summarize(matrix_path, maf_path, snploc_path)

    timings: list[float] = []
    last_summary: dict[str, Any] | None = None
    for _ in range(runs):
        t0 = time.perf_counter()
        last_summary = _load_and_summarize(matrix_path, maf_path, snploc_path)
        timings.append(time.perf_counter() - t0)

    assert last_summary is not None
    return {
        "benchmark": {
            "language": "Python",
            "target": "matrix_market_load_and_summary",
            "runs": int(runs),
            "warmup": int(warmup),
            "timings_seconds": [float(x) for x in timings],
            "median_seconds": float(statistics.median(timings)),
            "mean_seconds": float(statistics.mean(timings)),
            "min_seconds": float(min(timings)),
            "max_seconds": float(max(timings)),
        },
        "summary": last_summary,
    }


def _make_check(
    name: str,
    r_value: Any,
    py_value: Any,
    tolerance: float | None = None,
) -> dict[str, Any]:
    if tolerance is None:
        passed = r_value == py_value
        return {
            "name": name,
            "passed": bool(passed),
            "r_value": _to_builtin(r_value),
            "py_value": _to_builtin(py_value),
            "tolerance": None,
            "abs_diff": None,
        }

    abs_diff = abs(float(r_value) - float(py_value))
    passed = abs_diff <= float(tolerance)
    return {
        "name": name,
        "passed": bool(passed),
        "r_value": float(r_value),
        "py_value": float(py_value),
        "tolerance": float(tolerance),
        "abs_diff": float(abs_diff),
    }


def _build_parity(
    r_summary: dict[str, Any],
    py_summary: dict[str, Any],
) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []

    checks.append(_make_check("n_samples", r_summary["n_samples"], py_summary["n_samples"]))
    checks.append(_make_check("n_variants", r_summary["n_variants"], py_summary["n_variants"]))
    checks.append(_make_check("nnz", r_summary["nnz"], py_summary["nnz"]))
    checks.append(_make_check("genotype_sum", r_summary["genotype_sum"], py_summary["genotype_sum"], 1e-8))
    checks.append(_make_check("maf_mean", r_summary["maf_mean"], py_summary["maf_mean"], 1e-12))

    for key in ("min", "q01", "median", "q99", "max"):
        checks.append(
            _make_check(
                f"maf_summary.{key}",
                r_summary["maf_summary"][key],
                py_summary["maf_summary"][key],
                1e-12,
            )
        )

    checks.append(_make_check("maf_rows_vs_n_variants", r_summary["n_variants"], py_summary["maf_rows"]))
    checks.append(_make_check("snploc_rows_vs_n_variants", r_summary["n_variants"], py_summary["snploc_rows"]))
    checks.append(_make_check("variant_id_match", True, py_summary["variant_id_match"]))

    passed = all(check["passed"] for check in checks)
    failed = [check["name"] for check in checks if not check["passed"]]
    return {
        "scope": "data-level parity (genotype/snploc/maf)",
        "passed": bool(passed),
        "failed_checks": failed,
        "checks": checks,
    }


def _load_json(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _build_performance(
    r_benchmark: dict[str, Any] | None,
    py_benchmark: dict[str, Any],
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "scope": "matrix_market_load_and_summary",
        "python_median_seconds": py_benchmark["benchmark"]["median_seconds"],
        "python_runs": py_benchmark["benchmark"]["runs"],
        "python_warmup": py_benchmark["benchmark"]["warmup"],
    }

    if r_benchmark is None:
        payload["note"] = "R benchmark not provided; speedup unavailable."
        return payload

    r_median = float(r_benchmark["benchmark"]["median_seconds"])
    py_median = float(py_benchmark["benchmark"]["median_seconds"])
    speedup = None
    if py_median > 0:
        speedup = r_median / py_median

    payload.update(
        {
            "r_median_seconds": r_median,
            "r_runs": int(r_benchmark["benchmark"]["runs"]),
            "r_warmup": int(r_benchmark["benchmark"]["warmup"]),
            "speedup_r_div_py": speedup,
            "interpretation": (
                "speedup_r_div_py > 1 means Python is faster on this benchmark target."
            ),
        }
    )
    return payload


def _render_markdown(report: dict[str, Any]) -> str:
    parity = report["parity"]
    perf = report["performance"]

    lines: list[str] = []
    lines.append("# 1KG Parity Report (Local)")
    lines.append("")
    lines.append("## Scope")
    lines.append("- Data-level parity only (`genotype`, `snploc`, `maf`).")
    lines.append("- Performance target: Matrix Market load + summary statistics.")
    lines.append("- This is not a full STAAR workflow parity benchmark.")
    lines.append("")
    lines.append("## Correctness")
    lines.append(f"- Overall: {'PASS' if parity['passed'] else 'FAIL'}")
    if parity["failed_checks"]:
        lines.append(f"- Failed checks: {', '.join(parity['failed_checks'])}")
    else:
        lines.append("- Failed checks: none")
    lines.append("")
    lines.append("| Check | Status | Abs Diff | Tolerance |")
    lines.append("|---|---:|---:|---:|")
    for check in parity["checks"]:
        status = "PASS" if check["passed"] else "FAIL"
        abs_diff = "-" if check["abs_diff"] is None else f"{check['abs_diff']:.3e}"
        tol = "-" if check["tolerance"] is None else f"{check['tolerance']:.3e}"
        lines.append(f"| {check['name']} | {status} | {abs_diff} | {tol} |")
    lines.append("")
    lines.append("## Performance")
    lines.append(
        f"- Python median: {perf['python_median_seconds']:.4f}s "
        f"(runs={perf['python_runs']}, warmup={perf['python_warmup']})"
    )
    if "r_median_seconds" in perf:
        lines.append(
            f"- R median: {perf['r_median_seconds']:.4f}s "
            f"(runs={perf['r_runs']}, warmup={perf['r_warmup']})"
        )
        speedup = perf.get("speedup_r_div_py")
        if speedup is None:
            lines.append("- Speedup (R/Python): unavailable")
        else:
            lines.append(f"- Speedup (R/Python): {speedup:.3f}x")
        lines.append(f"- Note: {perf['interpretation']}")
    else:
        lines.append(f"- Note: {perf['note']}")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare R vs Python summaries and benchmark.")
    parser.add_argument("--r-summary", required=True, help="Path to r_summary.json from R export script.")
    parser.add_argument("--matrix", required=True, help="Path to genotype Matrix Market file (.mtx or .mtx.gz).")
    parser.add_argument("--maf", required=True, help="Path to maf CSV(.gz).")
    parser.add_argument("--snploc", required=True, help="Path to snploc CSV(.gz).")
    parser.add_argument("--r-benchmark", default=None, help="Optional path to R benchmark JSON.")
    parser.add_argument("--runs", type=int, default=5, help="Python measured runs.")
    parser.add_argument("--warmup", type=int, default=1, help="Python warmup runs.")
    parser.add_argument("--out-json", required=True, help="Output JSON report path.")
    parser.add_argument("--out-md", required=True, help="Output Markdown report path.")
    args = parser.parse_args()

    r_summary_path = Path(args.r_summary).expanduser().resolve()
    matrix_path = Path(args.matrix).expanduser().resolve()
    maf_path = Path(args.maf).expanduser().resolve()
    snploc_path = Path(args.snploc).expanduser().resolve()
    r_benchmark_path = (
        None if args.r_benchmark is None else Path(args.r_benchmark).expanduser().resolve()
    )
    out_json = Path(args.out_json).expanduser().resolve()
    out_md = Path(args.out_md).expanduser().resolve()
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)

    r_summary_raw = _load_json(r_summary_path)
    assert r_summary_raw is not None
    r_summary = {
        "n_samples": int(r_summary_raw["n_samples"]),
        "n_variants": int(r_summary_raw["n_variants"]),
        "nnz": int(r_summary_raw["nnz"]),
        "genotype_sum": float(r_summary_raw["genotype_sum"]),
        "maf_mean": float(r_summary_raw["maf_mean"]),
        "maf_summary": {
            key: float(r_summary_raw["maf_summary"][key]) for key in ("min", "q01", "median", "q99", "max")
        },
    }

    py_benchmark = _run_benchmark(
        matrix_path=matrix_path,
        maf_path=maf_path,
        snploc_path=snploc_path,
        runs=args.runs,
        warmup=args.warmup,
    )
    parity = _build_parity(r_summary=r_summary, py_summary=py_benchmark["summary"])

    r_benchmark = _load_json(r_benchmark_path) if r_benchmark_path is not None else None
    performance = _build_performance(r_benchmark=r_benchmark, py_benchmark=py_benchmark)

    report = {
        "meta": {
            "scope": "1kg_parity_local",
            "matrix_file": str(matrix_path),
            "maf_file": str(maf_path),
            "snploc_file": str(snploc_path),
            "r_summary_file": str(r_summary_path),
            "r_benchmark_file": None if r_benchmark_path is None else str(r_benchmark_path),
        },
        "parity": parity,
        "performance": performance,
        "python_benchmark": py_benchmark["benchmark"],
        "python_summary": py_benchmark["summary"],
        "r_summary": r_summary,
    }

    out_json.write_text(json.dumps(_to_builtin(report), indent=2, ensure_ascii=False), encoding="utf-8")
    out_md.write_text(_render_markdown(report), encoding="utf-8")
    print(f"Report written: {out_json}")
    print(f"Report written: {out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
