import json
import math
from pathlib import Path

import pytest
import yaml

from pystaar import workflows

ROOT = Path(__file__).resolve().parents[2]
SPECS_DIR = ROOT / "specs"
BASELINES_DIR = ROOT / "baselines"

ISSUE_BY_SCENARIO = {
    "ai_staar_related_dense_glmmkin_find_weight": "STAAR-25",
    "ai_staar_related_dense_glmmkin": "STAAR-22",
    "ai_staar_related_sparse_glmmkin_find_weight": "STAAR-24",
    "ai_staar_related_sparse_glmmkin": "STAAR-21",
    "ai_staar_unrelated_glm_find_weight": "STAAR-23",
    "ai_staar_unrelated_glm": "STAAR-20",
    "indiv_score_unrelated_glm": "STAAR-14",
    "indiv_score_related_sparse_glmmkin": "STAAR-15",
    "indiv_score_related_dense_glmmkin": "STAAR-16",
    "indiv_score_unrelated_glm_cond": "STAAR-17",
    "indiv_score_related_sparse_glmmkin_cond": "STAAR-18",
    "indiv_score_related_dense_glmmkin_cond": "STAAR-19",
    "staar_unrelated_glm": "STAAR-1",
    "staar_unrelated_glm_rare_maf_0_01": "STAAR-26",
    "staar_unrelated_glm_cond": "STAAR-4",
    "staar_unrelated_binary_spa": "STAAR-7",
    "staar_unrelated_binary_spa_case_q90": "STAAR-27",
    "staar_unrelated_binary_spa_filter": "STAAR-11",
    "staar_related_sparse_binary_spa": "STAAR-8",
    "staar_related_sparse_binary_spa_filter": "STAAR-12",
    "staar_related_dense_binary_spa": "STAAR-9",
    "staar_related_dense_binary_spa_filter": "STAAR-13",
    "staar_related_sparse_glmmkin": "STAAR-2",
    "staar_related_sparse_glmmkin_rare_maf_0_01": "STAAR-28",
    "staar_related_dense_glmmkin": "STAAR-3",
    "staar_related_dense_glmmkin_rare_maf_0_01": "STAAR-29",
    "staar_related_sparse_glmmkin_cond": "STAAR-5",
    "staar_related_dense_glmmkin_cond": "STAAR-6",
}


def _load_specs():
    return sorted(SPECS_DIR.glob("*.yaml"))


def _parse_entry_point(entry_point: str) -> str:
    # Expect format like: func_name(arg1, arg2=...)
    name = entry_point.split("(", 1)[0].strip()
    if not name:
        raise ValueError(f"Could not parse python_entry_point: {entry_point}")
    return name


def _load_json_path(file_path: Path, path: str):
    if not file_path.is_absolute():
        if file_path.parts and file_path.parts[0] == "baselines":
            file_path = ROOT / file_path
        else:
            file_path = BASELINES_DIR / file_path
    data = json.loads(file_path.read_text())
    if not path.startswith("$."):
        raise ValueError(f"Unsupported JSONPath: {path}")
    keys = [segment for segment in path[2:].split(".") if segment]
    for key in keys:
        if not isinstance(data, dict):
            raise KeyError(f"Path {path} not found in {file_path}")
        data = data[key]
    return data


def _assert_scalar_close(label: str, actual, expected, rtol: float, atol: float):
    if isinstance(expected, list):
        if len(expected) == 0:
            if actual != []:
                raise AssertionError(f"{label} expected empty list, got {actual}")
            return
        raise AssertionError(f"{label} expected scalar, got list: {expected}")
    if not math.isclose(actual, expected, rel_tol=rtol, abs_tol=atol):
        raise AssertionError(
            f"{label} mismatch: actual={actual} expected={expected} "
            f"(rtol={rtol}, atol={atol})"
        )


def _assert_mapping_close(label: str, actual, expected, rtol: float, atol: float):
    if set(actual.keys()) != set(expected.keys()):
        raise AssertionError(
            f"{label} keys mismatch: actual={sorted(actual.keys())} "
            f"expected={sorted(expected.keys())}"
        )
    for key in sorted(expected.keys()):
        _assert_scalar_close(f"{label}.{key}", actual[key], expected[key], rtol, atol)


@pytest.mark.parametrize("spec_path", _load_specs())
def test_scenario_parity(spec_path: Path):
    spec = yaml.safe_load(spec_path.read_text())
    scenario_id = spec["scenario_id"]
    issue = ISSUE_BY_SCENARIO.get(scenario_id, "STAAR-UNKNOWN")
    xfail_cfg = spec.get("xfail")
    xfail_issue = issue
    xfail_reason = ""
    if isinstance(xfail_cfg, dict):
        xfail_issue = str(xfail_cfg.get("issue", issue))
        xfail_reason = str(xfail_cfg.get("reason", "")).strip()

    func_name = _parse_entry_point(spec["python_entry_point"])
    func = getattr(workflows, func_name)

    dataset = spec["inputs"]["dataset"]
    params = spec.get("parameters", {})

    try:
        try:
            results = func(dataset, **params)
        except NotImplementedError as exc:
            if xfail_cfg:
                message = f"{scenario_id} expected-fail ({xfail_issue}): {exc}"
                if xfail_reason:
                    message = f"{message}; {xfail_reason}"
                pytest.xfail(message)
            pytest.xfail(f"{scenario_id} not implemented ({issue}): {exc}")

        if results is None:
            raise AssertionError(
                f"{scenario_id} returned None. Expected a dict of sentinel values."
            )

        for sentinel in spec["sentinels"]:
            expected = _load_json_path(Path(sentinel["file"]), sentinel["path"])
            actual = results.get(sentinel["name"]) if isinstance(results, dict) else None
            if actual is None:
                raise AssertionError(
                    f"{scenario_id} missing sentinel '{sentinel['name']}' in results"
                )

            comparison = sentinel["comparison"]
            rtol = float(sentinel.get("rtol", 0.0))
            atol = float(sentinel.get("atol", 0.0))

            if comparison == "scalar_approx":
                _assert_scalar_close(sentinel["name"], actual, expected, rtol, atol)
            elif comparison == "mapping_approx":
                _assert_mapping_close(sentinel["name"], actual, expected, rtol, atol)
            else:
                raise ValueError(
                    f"Unsupported comparison type: {comparison} in {spec_path}"
                )
    except AssertionError as exc:
        if xfail_cfg:
            message = f"{scenario_id} expected-fail ({xfail_issue}): {exc}"
            if xfail_reason:
                message = f"{message}; {xfail_reason}"
            pytest.xfail(message)
        raise
