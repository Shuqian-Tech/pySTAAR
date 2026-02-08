#!/usr/bin/env python3
"""Run wheel build + fresh-venv install smoke checks for release readiness."""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import tempfile
import textwrap
import time
import venv
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]


def _run(
    command: list[str],
    *,
    cwd: Path | None = None,
) -> dict[str, Any]:
    started = time.perf_counter()
    proc = subprocess.run(
        command,
        cwd=str(cwd) if cwd is not None else None,
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - started
    return {
        "command": command,
        "cwd": str(cwd) if cwd is not None else None,
        "returncode": proc.returncode,
        "seconds": elapsed,
        "stdout": proc.stdout.strip(),
        "stderr": proc.stderr.strip(),
    }


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _write_report(path: Path, payload: dict[str, Any]) -> None:
    lines: list[str] = []
    lines.append("# Release Readiness Smoke Check")
    lines.append("")
    lines.append(f"- Generated: {payload['generated_utc']}")
    lines.append(f"- Platform: {payload['platform']['system']} {payload['platform']['release']} ({payload['platform']['machine']})")
    lines.append(f"- Python: {payload['platform']['python_version']}")
    lines.append(f"- Overall status: {'PASS' if payload['ok'] else 'FAIL'}")
    lines.append("")
    lines.append("## Steps")
    lines.append("")
    for step in payload["steps"]:
        status = "PASS" if step["returncode"] == 0 else "FAIL"
        lines.append(
            f"- `{status}` `{step['label']}` in `{step['seconds']:.2f}s`: "
            + " ".join(step["command"])
        )
    lines.append("")
    lines.append("## Artifacts")
    lines.append("")
    lines.append("- JSON details: `reports/release_smoke.json`")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--json-output",
        type=Path,
        default=ROOT / "reports" / "release_smoke.json",
    )
    parser.add_argument(
        "--report-output",
        type=Path,
        default=ROOT / "reports" / "release_readiness.md",
    )
    args = parser.parse_args()

    generated_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    payload: dict[str, Any] = {
        "generated_utc": generated_utc,
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "python_version": platform.python_version(),
        },
        "steps": [],
        "ok": False,
    }

    with tempfile.TemporaryDirectory(prefix="pystaar_release_smoke_") as tmp:
        tmp_path = Path(tmp)
        wheel_dir = tmp_path / "wheelhouse"
        wheel_dir.mkdir(parents=True, exist_ok=True)

        step = _run(
            [sys.executable, "-m", "pip", "wheel", ".", "--no-deps", "-w", str(wheel_dir)],
            cwd=ROOT,
        )
        step["label"] = "build_wheel"
        payload["steps"].append(step)
        if step["returncode"] != 0:
            _write_json(args.json_output, payload)
            _write_report(args.report_output, payload)
            return 1

        wheels = sorted(wheel_dir.glob("pystaar-*.whl"))
        if not wheels:
            payload["steps"].append(
                {
                    "label": "locate_wheel",
                    "command": ["glob", "pystaar-*.whl"],
                    "cwd": str(wheel_dir),
                    "returncode": 1,
                    "seconds": 0.0,
                    "stdout": "",
                    "stderr": "No wheel generated.",
                }
            )
            _write_json(args.json_output, payload)
            _write_report(args.report_output, payload)
            return 1
        wheel_path = wheels[-1]

        venv_dir = tmp_path / "venv"
        venv.EnvBuilder(with_pip=True, clear=True).create(venv_dir)
        venv_python = venv_dir / "bin" / "python"
        if not venv_python.exists():
            venv_python = venv_dir / "Scripts" / "python.exe"

        step = _run([str(venv_python), "-m", "pip", "install", "--upgrade", "pip"])
        step["label"] = "upgrade_pip"
        payload["steps"].append(step)
        if step["returncode"] != 0:
            _write_json(args.json_output, payload)
            _write_report(args.report_output, payload)
            return 1

        step = _run([str(venv_python), "-m", "pip", "install", str(wheel_path)])
        step["label"] = "install_wheel"
        payload["steps"].append(step)
        if step["returncode"] != 0:
            _write_json(args.json_output, payload)
            _write_report(args.report_output, payload)
            return 1

        smoke_code = textwrap.dedent(
            """
            import pystaar
            info = pystaar.get_runtime_cache_info()
            assert "workflows" in info and "data" in info
            result = pystaar.staar_unrelated_glm(dataset="example", seed=600, rare_maf_cutoff=0.05)
            assert "results_STAAR_O" in result
            assert 0.0 <= float(result["results_STAAR_O"]) <= 1.0
            after = pystaar.clear_runtime_caches(include_dataset_cache=True)
            assert after["after"]["workflows"]["related_nullmodel"]["currsize"] == 0
            print("SMOKE_OK", float(result["results_STAAR_O"]))
            """
        ).strip()
        step = _run([str(venv_python), "-c", smoke_code])
        step["label"] = "import_and_workflow_smoke"
        payload["steps"].append(step)
        if step["returncode"] != 0:
            _write_json(args.json_output, payload)
            _write_report(args.report_output, payload)
            return 1

    payload["ok"] = all(step["returncode"] == 0 for step in payload["steps"])
    _write_json(args.json_output, payload)
    _write_report(args.report_output, payload)
    return 0 if payload["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
