import math
import shutil
from types import SimpleNamespace

import pytest

from pystaar import workflows
from pystaar import staar_core
from pystaar.data import load_example_dataset
from pystaar.models import NullModelGLMMKin


def _load_runtime_ai_metadata_from_example_files():
    groups_df = workflows.pd.read_csv(workflows.DATA_DIR / workflows.AI_POP_GROUPS_FILE)
    pop_groups = groups_df["pop_group"].astype(str).to_numpy()
    pop_levels = list(dict.fromkeys(pop_groups.tolist()))

    w11_df = workflows.pd.read_csv(workflows.DATA_DIR / workflows.AI_POP_WEIGHTS_1_1_FILE)
    w125_df = workflows.pd.read_csv(workflows.DATA_DIR / workflows.AI_POP_WEIGHTS_1_25_FILE)
    weight_cols = [c for c in w11_df.columns if c != "population"]
    w11 = w11_df.set_index("population").loc[pop_levels, weight_cols].to_numpy(dtype=float)
    w125 = w125_df.set_index("population").loc[pop_levels, weight_cols].to_numpy(dtype=float)
    return pop_groups, pop_levels, w11, w125


def _copy_example_dataset_to_directory(target_dir):
    file_map = {
        "example_geno.mtx": "geno.mtx",
        "example_phred.csv": "phred.csv",
        "example_pheno_unrelated.csv": "pheno_unrelated.csv",
        "example_pheno_related.csv": "pheno_related.csv",
        "example_kins_sparse.mtx": "kins_sparse.mtx",
        "example_kins_dense.mtx": "kins_dense.mtx",
    }
    target_dir.mkdir(parents=True, exist_ok=True)
    for src_name, dst_name in file_map.items():
        shutil.copy2(workflows.DATA_DIR / src_name, target_dir / dst_name)
    return target_dir


def test_related_workflow_skips_precomputed_cov_for_nonbaseline_cutoff(monkeypatch):
    original_read_csv = workflows.pd.read_csv
    blocked_files = {"example_glmmkin_cov.csv", "example_glmmkin_scaled_residuals.csv"}

    def _guard_precomputed_reads(*args, **kwargs):
        path = str(args[0]) if args else ""
        if any(name in path for name in blocked_files):
            raise AssertionError("precomputed artifacts should not be loaded")
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(workflows.pd, "read_csv", _guard_precomputed_reads)

    results = workflows.staar_related_sparse_glmmkin(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.01,
    )

    assert results["num_variant"] > 0
    assert results["num_variant"] != 163.0


def test_unrelated_workflow_runs_with_directory_dataset(tmp_path):
    dataset_dir = _copy_example_dataset_to_directory(tmp_path / "dataset_copy")

    results = workflows.staar_unrelated_glm(
        dataset=str(dataset_dir),
        seed=600,
        rare_maf_cutoff=0.05,
    )

    assert results["num_variant"] == 163.0
    assert 0.0 <= results["results_STAAR_O"] <= 1.0


def test_related_binary_workflow_runs_with_dataset_object():
    dataset_obj = load_example_dataset()

    results = workflows.staar_related_sparse_binary_spa(
        dataset=dataset_obj,
        seed=600,
        rare_maf_cutoff=0.05,
        case_quantile=0.90,
        use_precomputed_artifacts=False,
    )

    assert results["num_variant"] == 163.0
    assert results["case_count"] > 0.0
    assert 0.0 < results["results_STAAR_B"] <= 1.0


def test_related_workflow_precomputed_mode_does_not_load_glmm_cov_artifacts(monkeypatch):
    original_read_csv = workflows.pd.read_csv

    def _guard_precomputed_reads(*args, **kwargs):
        path = str(args[0]) if args else ""
        if path.endswith("example_glmmkin_cov.csv") or path.endswith(
            "example_glmmkin_cov_rare_maf_0_01.csv"
        ):
            raise AssertionError("related GLMM precomputed mode should not load covariance artifacts")
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(workflows.pd, "read_csv", _guard_precomputed_reads)

    results = workflows.staar_related_sparse_glmmkin(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.01,
        use_precomputed_artifacts=True,
    )

    assert results["num_variant"] == 153.0


def test_related_nullmodel_cache_reuses_fit_for_string_dataset(monkeypatch):
    dummy_data = SimpleNamespace(
        geno=workflows.np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 1.0]], dtype=float),
        phred=workflows.np.zeros((3, 1), dtype=float),
        pheno_related=workflows.pd.DataFrame(
            {
                "Y": [0.0, 1.0],
                "X1": [0.0, 1.0],
                "X2": [1.0, 0.0],
            }
        ),
        kins_sparse=workflows.sp.identity(2, format="csc"),
        kins_dense=workflows.sp.identity(2, format="csc"),
    )

    fit_calls = {"count": 0}

    def _fake_fit_null_glmmkin(*args, **kwargs):
        fit_calls["count"] += 1
        zeros = workflows.np.zeros(2, dtype=float)
        return NullModelGLMMKin(
            X=workflows.np.zeros((2, 3), dtype=float),
            y=zeros,
            fitted=zeros,
            residuals=zeros,
            scaled_residuals=zeros,
            beta=workflows.np.zeros(3, dtype=float),
            theta=workflows.np.array([1.0, 1.0], dtype=float),
            cov=workflows.np.eye(3, dtype=float),
            sigma_solver=None,
            Sigma_iX=workflows.np.zeros((2, 3), dtype=float),
            sparse_kins=True,
        )

    def _fake_staar(*args, **kwargs):
        return {
            "num_variant": 2.0,
            "results_STAAR_O": 0.11,
            "results_STAAR_S_1_25": {"STAAR-S(1,25)": 0.21},
            "results_STAAR_S_1_1": {"STAAR-S(1,1)": 0.22},
            "results_STAAR_B_1_25": {"STAAR-B(1,25)": 0.31},
            "results_STAAR_B_1_1": {"STAAR-B(1,1)": 0.32},
            "results_STAAR_A_1_25": {"STAAR-A(1,25)": 0.41},
            "results_STAAR_A_1_1": {"STAAR-A(1,1)": 0.42},
        }

    monkeypatch.setattr(workflows, "load_dataset", lambda dataset: dummy_data)
    monkeypatch.setattr(workflows, "fit_null_glmmkin", _fake_fit_null_glmmkin)
    monkeypatch.setattr(workflows, "_ORIGINAL_FIT_NULL_GLMMKIN", _fake_fit_null_glmmkin)
    monkeypatch.setattr(workflows, "staar", _fake_staar)

    workflows._fit_related_nullmodel_cached.cache_clear()
    try:
        first = workflows.staar_related_sparse_glmmkin(
            dataset="example",
            use_precomputed_artifacts=False,
        )
        second = workflows.staar_related_sparse_glmmkin(
            dataset="example",
            use_precomputed_artifacts=False,
        )
    finally:
        workflows._fit_related_nullmodel_cached.cache_clear()

    assert fit_calls["count"] == 1
    assert first["num_variant"] == 2.0
    assert second["results_STAAR_O"] == pytest.approx(0.11)


def test_related_ai_cache_reuses_ai_staar_for_string_dataset(monkeypatch):
    dummy_data = SimpleNamespace(
        geno=workflows.np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 1.0]], dtype=float),
        phred=workflows.np.zeros((3, 1), dtype=float),
        pheno_related=workflows.pd.DataFrame(
            {
                "Y": [0.0, 1.0],
                "X1": [0.0, 1.0],
                "X2": [1.0, 0.0],
            }
        ),
        kins_sparse=workflows.sp.identity(2, format="csc"),
        kins_dense=workflows.sp.identity(2, format="csc"),
    )

    fit_calls = {"count": 0}
    ai_calls = {"count": 0}

    def _fake_fit_null_glmmkin(*args, **kwargs):
        fit_calls["count"] += 1
        zeros = workflows.np.zeros(2, dtype=float)
        return NullModelGLMMKin(
            X=workflows.np.zeros((2, 3), dtype=float),
            y=zeros,
            fitted=zeros,
            residuals=zeros,
            scaled_residuals=zeros,
            beta=workflows.np.zeros(3, dtype=float),
            theta=workflows.np.array([1.0, 1.0], dtype=float),
            cov=workflows.np.eye(3, dtype=float),
            sigma_solver=None,
            Sigma_iX=workflows.np.zeros((2, 3), dtype=float),
            sparse_kins=True,
        )

    def _fake_ai_staar(*args, **kwargs):
        ai_calls["count"] += 1
        return {
            "num_variant": 2.0,
            "cMAC": 3.0,
            "results_STAAR_O": 0.12,
            "results_ACAT_O": 0.13,
            "results_STAAR_S_1_25": {"STAAR-S(1,25)": 0.22},
            "results_STAAR_S_1_1": {"STAAR-S(1,1)": 0.23},
            "results_STAAR_B_1_25": {"STAAR-B(1,25)": 0.32},
            "results_STAAR_B_1_1": {"STAAR-B(1,1)": 0.33},
            "results_STAAR_A_1_25": {"STAAR-A(1,25)": 0.42},
            "results_STAAR_A_1_1": {"STAAR-A(1,1)": 0.43},
        }

    monkeypatch.setattr(workflows, "load_dataset", lambda dataset: dummy_data)
    monkeypatch.setattr(workflows, "fit_null_glmmkin", _fake_fit_null_glmmkin)
    monkeypatch.setattr(workflows, "_ORIGINAL_FIT_NULL_GLMMKIN", _fake_fit_null_glmmkin)
    monkeypatch.setattr(workflows, "ai_staar", _fake_ai_staar)
    monkeypatch.setattr(workflows, "_ORIGINAL_AI_STAAR", _fake_ai_staar)
    monkeypatch.setattr(
        workflows,
        "_resolve_ai_metadata",
        lambda **kwargs: (
            workflows.np.array(["EUR", "EUR"]),
            ["EUR"],
            workflows.np.ones((1, 2), dtype=float),
            workflows.np.ones((1, 2), dtype=float),
        ),
    )

    workflows._fit_related_nullmodel_cached.cache_clear()
    workflows._related_ai_common_cached.cache_clear()
    try:
        first = workflows.ai_staar_related_sparse_glmmkin(
            dataset="example",
            use_precomputed_artifacts=False,
        )
        second = workflows.ai_staar_related_sparse_glmmkin(
            dataset="example",
            use_precomputed_artifacts=False,
        )
    finally:
        workflows._fit_related_nullmodel_cached.cache_clear()
        workflows._related_ai_common_cached.cache_clear()

    assert fit_calls["count"] == 1
    assert ai_calls["count"] == 1
    assert first["results_STAAR_O"] == pytest.approx(0.12)
    assert second["results_STAAR_O"] == pytest.approx(0.12)


def test_related_workflow_precomputed_mode_uses_scaled_without_cov_for_baseline_cutoff(
    monkeypatch,
):
    original_fit_null_glmmkin = workflows.fit_null_glmmkin
    captured = {}

    def _track_precomputed_inputs(*args, **kwargs):
        captured["precomputed_cov"] = kwargs.get("precomputed_cov")
        captured["precomputed_scaled_residuals"] = kwargs.get("precomputed_scaled_residuals")
        return original_fit_null_glmmkin(*args, **kwargs)

    monkeypatch.setattr(workflows, "fit_null_glmmkin", _track_precomputed_inputs)

    results = workflows.staar_related_sparse_glmmkin(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        use_precomputed_artifacts=True,
    )

    assert results["num_variant"] == 163.0
    assert captured["precomputed_cov"] is None
    assert captured["precomputed_scaled_residuals"] is not None


def test_related_workflow_precomputed_mode_omits_cov_for_nonbaseline_cutoff(monkeypatch):
    original_fit_null_glmmkin = workflows.fit_null_glmmkin
    captured = {}

    def _track_precomputed_inputs(*args, **kwargs):
        captured["precomputed_cov"] = kwargs.get("precomputed_cov")
        captured["precomputed_scaled_residuals"] = kwargs.get("precomputed_scaled_residuals")
        return original_fit_null_glmmkin(*args, **kwargs)

    monkeypatch.setattr(workflows, "fit_null_glmmkin", _track_precomputed_inputs)

    results = workflows.staar_related_sparse_glmmkin(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.01,
        use_precomputed_artifacts=True,
    )

    assert results["num_variant"] == 153.0
    assert captured["precomputed_cov"] is None
    assert captured["precomputed_scaled_residuals"] is not None


def test_related_workflow_precomputed_mode_does_not_load_cov_for_baseline_cutoff(monkeypatch):
    original_read_csv = workflows.pd.read_csv

    def _guard_cov_reads(*args, **kwargs):
        path = str(args[0]) if args else ""
        if path.endswith("example_glmmkin_cov.csv"):
            raise AssertionError("baseline related GLMM precomputed mode should not load covariance artifacts")
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(workflows.pd, "read_csv", _guard_cov_reads)

    results = workflows.staar_related_sparse_glmmkin(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        use_precomputed_artifacts=True,
    )

    assert results["num_variant"] == 163.0


def test_related_workflow_defaults_to_pure_path_for_baseline_cutoff(monkeypatch):
    original_fit_null_glmmkin = workflows.fit_null_glmmkin
    captured = {}

    def _track_precomputed_inputs(*args, **kwargs):
        captured["precomputed_cov"] = kwargs.get("precomputed_cov")
        captured["precomputed_scaled_residuals"] = kwargs.get("precomputed_scaled_residuals")
        return original_fit_null_glmmkin(*args, **kwargs)

    monkeypatch.setattr(workflows, "fit_null_glmmkin", _track_precomputed_inputs)

    results = workflows.staar_related_sparse_glmmkin(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
    )

    assert results["num_variant"] == 163.0
    assert captured["precomputed_cov"] is None
    assert captured["precomputed_scaled_residuals"] is None


def test_related_ai_workflow_passes_precomputed_theta_in_precomputed_mode(monkeypatch):
    original_fit_null_glmmkin = workflows.fit_null_glmmkin
    captured = {}

    def _track_precomputed_inputs(*args, **kwargs):
        captured["precomputed_theta"] = kwargs.get("precomputed_theta")
        return original_fit_null_glmmkin(*args, **kwargs)

    monkeypatch.setattr(workflows, "fit_null_glmmkin", _track_precomputed_inputs)

    workflows.ai_staar_related_sparse_glmmkin(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        use_precomputed_artifacts=True,
    )

    assert captured["precomputed_theta"] is not None
    assert captured["precomputed_theta"][0] == pytest.approx(
        workflows.BASELINE_GLMMKIN_THETA_SPARSE[0]
    )
    assert captured["precomputed_theta"][1] == pytest.approx(
        workflows.BASELINE_GLMMKIN_THETA_SPARSE[1]
    )


def test_related_workflow_omits_precomputed_theta_in_pure_mode(monkeypatch):
    original_fit_null_glmmkin = workflows.fit_null_glmmkin
    captured = {}

    def _track_precomputed_inputs(*args, **kwargs):
        captured["precomputed_theta"] = kwargs.get("precomputed_theta")
        return original_fit_null_glmmkin(*args, **kwargs)

    monkeypatch.setattr(workflows, "fit_null_glmmkin", _track_precomputed_inputs)

    workflows.staar_related_sparse_glmmkin(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        use_precomputed_artifacts=False,
    )

    assert captured["precomputed_theta"] is None


def test_unrelated_binary_spa_filter_cutoff_one_matches_full_spa():
    full_spa = workflows.staar_unrelated_binary_spa(
        dataset="example",
        seed=600,
        SPA_p_filter=False,
    )
    filtered = workflows.staar_unrelated_binary_spa(
        dataset="example",
        seed=600,
        SPA_p_filter=True,
        p_filter_cutoff=1.0,
    )

    assert full_spa["num_variant"] == filtered["num_variant"]
    assert full_spa["cMAC"] == filtered["cMAC"]
    assert full_spa["case_count"] == filtered["case_count"]
    assert math.isclose(full_spa["results_STAAR_B"], filtered["results_STAAR_B"], rel_tol=1e-10, abs_tol=1e-12)
    for key, expected in full_spa["results_STAAR_B_1_25"].items():
        assert math.isclose(filtered["results_STAAR_B_1_25"][key], expected, rel_tol=1e-10, abs_tol=1e-12)
    for key, expected in full_spa["results_STAAR_B_1_1"].items():
        assert math.isclose(filtered["results_STAAR_B_1_1"][key], expected, rel_tol=1e-10, abs_tol=1e-12)


def test_unrelated_binary_spa_invalid_filter_cutoff():
    with pytest.raises(ValueError, match="p_filter_cutoff"):
        workflows.staar_unrelated_binary_spa(
            dataset="example",
            seed=600,
            SPA_p_filter=True,
            p_filter_cutoff=0.0,
        )


def test_related_binary_spa_filter_cutoff_one_matches_full_spa():
    full_spa = workflows.staar_related_sparse_binary_spa(
        dataset="example",
        seed=600,
        SPA_p_filter=False,
    )
    filtered = workflows.staar_related_sparse_binary_spa(
        dataset="example",
        seed=600,
        SPA_p_filter=True,
        p_filter_cutoff=1.0,
    )

    assert full_spa["num_variant"] == filtered["num_variant"]
    assert full_spa["cMAC"] == filtered["cMAC"]
    assert full_spa["case_count"] == filtered["case_count"]
    assert math.isclose(full_spa["results_STAAR_B"], filtered["results_STAAR_B"], rel_tol=1e-10, abs_tol=1e-12)
    for key, expected in full_spa["results_STAAR_B_1_25"].items():
        assert math.isclose(filtered["results_STAAR_B_1_25"][key], expected, rel_tol=1e-10, abs_tol=1e-12)
    for key, expected in full_spa["results_STAAR_B_1_1"].items():
        assert math.isclose(filtered["results_STAAR_B_1_1"][key], expected, rel_tol=1e-10, abs_tol=1e-12)


@pytest.mark.parametrize(
    "workflow",
    [
        workflows.staar_related_sparse_binary_spa,
        workflows.staar_related_dense_binary_spa,
    ],
)
def test_related_binary_spa_supports_nonbaseline_case_quantile(workflow):
    results = workflow(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        case_quantile=0.90,
        use_precomputed_artifacts=False,
    )

    assert results["case_count"] > 0.0
    assert 0.0 < results["results_STAAR_B"] <= 1.0


@pytest.mark.parametrize(
    "workflow",
    [
        workflows.staar_related_sparse_binary_spa,
        workflows.staar_related_dense_binary_spa,
    ],
)
def test_related_binary_spa_nonbaseline_case_quantile_skips_precomputed_scaled_residuals(
    monkeypatch,
    workflow,
):
    original_read_csv = workflows.pd.read_csv

    def _guard_scaled_residual_reads(*args, **kwargs):
        path = str(args[0]) if args else ""
        if path.endswith("example_glmmkin_binary_spa_sparse_scaled_residuals.csv"):
            raise AssertionError(
                "nonbaseline case_quantile should not use precomputed scaled residual artifact"
            )
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(workflows.pd, "read_csv", _guard_scaled_residual_reads)

    results = workflow(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        case_quantile=0.90,
        use_precomputed_artifacts=True,
    )

    assert results["case_count"] > 0.0
    assert 0.0 < results["results_STAAR_B"] <= 1.0


def test_related_binary_spa_prefilter_loads_precomputed_cov(monkeypatch):
    original_staar_binary_spa = workflows.staar_binary_spa
    captured = {}

    def _track_staar_binary_spa(*args, **kwargs):
        obj = kwargs.get("obj_nullmodel")
        captured["has_precomputed_cov_filter"] = (
            obj is not None and getattr(obj, "precomputed_cov_filter", None) is not None
        )
        return original_staar_binary_spa(*args, **kwargs)

    monkeypatch.setattr(workflows, "staar_binary_spa", _track_staar_binary_spa)

    workflows.staar_related_sparse_binary_spa(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        SPA_p_filter=True,
        p_filter_cutoff=0.05,
        use_precomputed_artifacts=True,
    )

    assert captured["has_precomputed_cov_filter"] is True


@pytest.mark.parametrize(
    "workflow,blocked_file",
    [
        (
            workflows.staar_related_sparse_binary_spa,
            "example_glmmkin_binary_spa_sparse_cov_filter.csv",
        ),
        (
            workflows.staar_related_dense_binary_spa,
            "example_glmmkin_binary_spa_dense_cov_filter.csv",
        ),
    ],
)
def test_related_binary_spa_prefilter_does_not_load_precomputed_cov_file(
    monkeypatch,
    workflow,
    blocked_file,
):
    original_read_csv = workflows.pd.read_csv

    def _guard_cov_filter_reads(*args, **kwargs):
        path = str(args[0]) if args else ""
        if path.endswith(blocked_file):
            raise AssertionError("related SPA prefilter should not load precomputed cov file")
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(workflows.pd, "read_csv", _guard_cov_filter_reads)

    results = workflow(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        SPA_p_filter=True,
        p_filter_cutoff=0.05,
        use_precomputed_artifacts=True,
    )

    assert results["num_variant"] == 163.0
    assert 0.0 < results["results_STAAR_B"] <= 1.0


@pytest.mark.parametrize(
    "workflow,blocked_files",
    [
        (
            workflows.staar_related_sparse_binary_spa,
            (
                "example_glmmkin_binary_spa_sparse_XW.csv",
                "example_glmmkin_binary_spa_sparse_XXWX_inv.csv",
            ),
        ),
        (
            workflows.staar_related_dense_binary_spa,
            (
                "example_glmmkin_binary_spa_dense_XW.csv",
                "example_glmmkin_binary_spa_dense_XXWX_inv.csv",
            ),
        ),
    ],
)
def test_related_binary_spa_precomputed_path_reconstructs_nullmodel_components(
    monkeypatch,
    workflow,
    blocked_files,
):
    original_read_csv = workflows.pd.read_csv

    def _guard_component_reads(*args, **kwargs):
        path = str(args[0]) if args else ""
        if any(path.endswith(name) for name in blocked_files):
            raise AssertionError("related SPA precomputed path should reconstruct components")
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(workflows.pd, "read_csv", _guard_component_reads)

    results = workflow(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        SPA_p_filter=True,
        p_filter_cutoff=0.05,
        use_precomputed_artifacts=True,
    )

    assert results["num_variant"] == 163.0
    assert 0.0 < results["results_STAAR_B"] <= 1.0


def test_related_binary_spa_dense_precomputed_path_uses_shared_scaled_residual_artifact(
    monkeypatch,
):
    original_read_csv = workflows.pd.read_csv

    def _guard_dense_scaled_residual_reads(*args, **kwargs):
        path = str(args[0]) if args else ""
        if path.endswith("example_glmmkin_binary_spa_dense_scaled_residuals.csv"):
            raise AssertionError("dense related SPA path should use shared scaled residual artifact")
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(workflows.pd, "read_csv", _guard_dense_scaled_residual_reads)

    results = workflows.staar_related_dense_binary_spa(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        SPA_p_filter=True,
        p_filter_cutoff=0.05,
        use_precomputed_artifacts=True,
    )

    assert results["num_variant"] == 163.0
    assert 0.0 < results["results_STAAR_B"] <= 1.0


@pytest.mark.parametrize(
    "workflow",
    [
        workflows.staar_related_sparse_binary_spa,
        workflows.staar_related_dense_binary_spa,
    ],
)
def test_related_binary_spa_precomputed_path_does_not_load_fitted_artifacts(
    monkeypatch,
    workflow,
):
    original_read_csv = workflows.pd.read_csv

    def _guard_fitted_reads(*args, **kwargs):
        path = str(args[0]) if args else ""
        if path.endswith("example_glmmkin_binary_spa_sparse_fitted.csv") or path.endswith(
            "example_glmmkin_binary_spa_dense_fitted.csv"
        ):
            raise AssertionError("related SPA precomputed path should not load fitted artifacts")
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(workflows.pd, "read_csv", _guard_fitted_reads)

    results = workflow(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        SPA_p_filter=True,
        p_filter_cutoff=0.05,
        use_precomputed_artifacts=True,
    )

    assert results["num_variant"] == 163.0
    assert 0.0 < results["results_STAAR_B"] <= 1.0


def test_related_binary_spa_prefilter_defaults_to_pure_covariance(monkeypatch):
    original_staar_binary_spa = workflows.staar_binary_spa
    captured = {}

    def _track_staar_binary_spa(*args, **kwargs):
        obj = kwargs.get("obj_nullmodel")
        captured["has_precomputed_cov_filter"] = (
            obj is not None and getattr(obj, "precomputed_cov_filter", None) is not None
        )
        return original_staar_binary_spa(*args, **kwargs)

    monkeypatch.setattr(workflows, "staar_binary_spa", _track_staar_binary_spa)

    workflows.staar_related_sparse_binary_spa(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        SPA_p_filter=True,
        p_filter_cutoff=0.05,
    )

    assert captured["has_precomputed_cov_filter"] is False


@pytest.mark.parametrize(
    "workflow",
    [
        workflows.staar_related_sparse_binary_spa,
        workflows.staar_related_dense_binary_spa,
    ],
)
def test_related_binary_spa_pure_path_is_not_degenerate(workflow):
    results = workflow(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        use_precomputed_artifacts=False,
    )

    assert 0.0 < results["results_STAAR_B"] <= 1.0


@pytest.mark.parametrize(
    "workflow,spa_filter",
    [
        (workflows.staar_related_sparse_binary_spa, False),
        (workflows.staar_related_dense_binary_spa, False),
        (workflows.staar_related_sparse_binary_spa, True),
        (workflows.staar_related_dense_binary_spa, True),
    ],
)
def test_related_binary_spa_pure_path_tracks_precomputed(workflow, spa_filter):
    kwargs = {"SPA_p_filter": spa_filter}
    if spa_filter:
        kwargs["p_filter_cutoff"] = 0.05

    pure = workflow(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        use_precomputed_artifacts=False,
        **kwargs,
    )
    precomputed = workflow(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        use_precomputed_artifacts=True,
        **kwargs,
    )

    assert abs(pure["results_STAAR_B"] - precomputed["results_STAAR_B"]) < 2.5e-6
    assert (
        abs(
            pure["results_STAAR_B_1_25"]["STAAR-B(1,25)"]
            - precomputed["results_STAAR_B_1_25"]["STAAR-B(1,25)"]
        )
        < 2.5e-6
    )


def test_unrelated_binary_spa_prefilter_does_not_load_precomputed_cov(monkeypatch):
    original_read_csv = workflows.pd.read_csv

    def _guard_cov_filter_reads(*args, **kwargs):
        path = str(args[0]) if args else ""
        if path.endswith("example_glm_binary_spa_cov_filter.csv"):
            raise AssertionError("unrelated SPA prefilter should not load precomputed cov")
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(workflows.pd, "read_csv", _guard_cov_filter_reads)

    results = workflows.staar_unrelated_binary_spa(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        SPA_p_filter=True,
        p_filter_cutoff=0.05,
    )

    assert results["num_variant"] == 163.0
    assert 0.0 < results["results_STAAR_B"] <= 1.0


def test_indiv_score_unrelated_workflow_runs():
    results = workflows.indiv_score_unrelated_glm(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
    )
    assert results["num_variant"] > 0
    assert results["num_tested"] > 0
    assert len(results["pvalue_samples"]) > 0


def test_indiv_score_related_cond_workflow_runs():
    results = workflows.indiv_score_related_sparse_glmmkin_cond(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
    )
    assert results["num_variant"] > 0
    assert results["num_tested"] > 0
    assert len(results["pvalue_cond_samples"]) > 0


@pytest.mark.parametrize(
    "workflow",
    [
        workflows.staar_related_sparse_glmmkin_cond,
        workflows.staar_related_dense_glmmkin_cond,
    ],
)
def test_related_cond_workflow_does_not_load_precomputed_cond_cov(
    monkeypatch,
    workflow,
):
    original_read_csv = workflows.pd.read_csv

    def _guard_cond_cov_reads(*args, **kwargs):
        path = str(args[0]) if args else ""
        if path.endswith("example_glmmkin_cov_cond_sparse.csv") or path.endswith(
            "example_glmmkin_cov_cond_dense.csv"
        ):
            raise AssertionError("related conditional workflow should not load cond covariance artifacts")
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(workflows.pd, "read_csv", _guard_cond_cov_reads)

    results = workflow(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        use_precomputed_artifacts=True,
    )
    assert results["num_variant"] == 163.0
    assert 0.0 <= results["results_STAAR_O_cond"] <= 1.0


@pytest.mark.parametrize(
    "workflow",
    [
        workflows.indiv_score_related_sparse_glmmkin_cond,
        workflows.indiv_score_related_dense_glmmkin_cond,
    ],
)
def test_related_indiv_cond_workflow_does_not_load_precomputed_cond_cov(
    monkeypatch,
    workflow,
):
    original_read_csv = workflows.pd.read_csv

    def _guard_cond_cov_reads(*args, **kwargs):
        path = str(args[0]) if args else ""
        if path.endswith("example_glmmkin_cov_cond_sparse.csv") or path.endswith(
            "example_glmmkin_cov_cond_dense.csv"
        ):
            raise AssertionError("related indiv conditional workflow should not load cond covariance artifacts")
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(workflows.pd, "read_csv", _guard_cond_cov_reads)

    results = workflow(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        use_precomputed_artifacts=True,
    )
    assert results["num_variant"] == 163.0
    assert results["num_tested"] > 0


def test_ai_staar_unrelated_workflow_runs():
    results = workflows.ai_staar_unrelated_glm(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
    )
    assert results["num_variant"] > 0
    assert results["cMAC"] > 0
    assert 0.0 <= results["results_STAAR_O"] <= 1.0


def test_ai_staar_unrelated_uses_group_stat_covariance_fast_path(monkeypatch):
    original = staar_core._ai_unrelated_gaussian_cov_from_group_stats
    called = {"count": 0}

    def _wrapped(*args, **kwargs):
        called["count"] += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(staar_core, "_ai_unrelated_gaussian_cov_from_group_stats", _wrapped)

    results = workflows.ai_staar_unrelated_glm(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
    )

    assert called["count"] > 0
    assert 0.0 <= results["results_STAAR_O"] <= 1.0


def test_ai_staar_unrelated_requires_runtime_metadata_for_nonexample_dataset():
    dataset_obj = load_example_dataset()

    with pytest.raises(ValueError, match="AI metadata must be provided for non-example datasets"):
        workflows.ai_staar_unrelated_glm(
            dataset=dataset_obj,
            seed=600,
            rare_maf_cutoff=0.05,
        )


def test_ai_staar_unrelated_accepts_runtime_metadata_for_nonexample_dataset():
    dataset_obj = load_example_dataset()
    pop_groups, pop_levels, w11, w125 = _load_runtime_ai_metadata_from_example_files()

    results = workflows.ai_staar_unrelated_glm(
        dataset=dataset_obj,
        seed=600,
        rare_maf_cutoff=0.05,
        pop_groups=pop_groups,
        pop_levels=pop_levels,
        pop_weights_1_1=w11,
        pop_weights_1_25=w125,
    )

    assert results["num_variant"] > 0
    assert results["cMAC"] > 0
    assert 0.0 <= results["results_STAAR_O"] <= 1.0


def test_ai_staar_related_sparse_workflow_runs():
    results = workflows.ai_staar_related_sparse_glmmkin(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
    )
    assert results["num_variant"] > 0
    assert results["cMAC"] > 0
    assert 0.0 <= results["results_STAAR_O"] <= 1.0


def test_ai_staar_related_accepts_runtime_metadata_for_nonexample_dataset():
    dataset_obj = load_example_dataset()
    pop_groups, pop_levels, w11, w125 = _load_runtime_ai_metadata_from_example_files()

    results = workflows.ai_staar_related_sparse_glmmkin(
        dataset=dataset_obj,
        seed=600,
        rare_maf_cutoff=0.05,
        pop_groups=pop_groups,
        pop_levels=pop_levels,
        pop_weights_1_1=w11,
        pop_weights_1_25=w125,
    )

    assert results["num_variant"] > 0
    assert results["cMAC"] > 0
    assert 0.0 <= results["results_STAAR_O"] <= 1.0


def test_ai_staar_related_dense_workflow_runs():
    results = workflows.ai_staar_related_dense_glmmkin(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
    )
    assert results["num_variant"] > 0
    assert results["cMAC"] > 0
    assert 0.0 <= results["results_STAAR_O"] <= 1.0


@pytest.mark.parametrize(
    "workflow",
    [
        workflows.ai_staar_related_sparse_glmmkin,
        workflows.ai_staar_related_dense_glmmkin,
    ],
)
def test_ai_staar_related_workflow_does_not_load_precomputed_ai_cov_artifacts(
    monkeypatch,
    workflow,
):
    original_read_csv = workflows.pd.read_csv

    def _guard_ai_cov_reads(*args, **kwargs):
        path = str(args[0]) if args else ""
        if "example_ai_cov_" in path:
            raise AssertionError("related AI workflow should not load precomputed AI covariance artifacts")
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(workflows.pd, "read_csv", _guard_ai_cov_reads)

    results = workflow(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        use_precomputed_artifacts=True,
    )
    assert results["num_variant"] == 163.0
    assert 0.0 <= results["results_STAAR_O"] <= 1.0


def test_ai_staar_unrelated_find_weight_workflow_runs():
    results = workflows.ai_staar_unrelated_glm_find_weight(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
    )
    assert results["num_variant"] > 0
    assert results["cMAC"] > 0
    assert 0.0 <= results["results_STAAR_O"] <= 1.0
    assert len(results["weight_all_1"]) > 0
    assert len(results["results_weight_staar_o"]) > 0


def test_ai_staar_related_sparse_find_weight_workflow_runs():
    results = workflows.ai_staar_related_sparse_glmmkin_find_weight(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
    )
    assert results["num_variant"] > 0
    assert results["cMAC"] > 0
    assert 0.0 <= results["results_STAAR_O"] <= 1.0
    assert len(results["weight_all_1"]) > 0
    assert len(results["results_weight_staar_o"]) > 0


def test_ai_staar_related_dense_find_weight_workflow_runs():
    results = workflows.ai_staar_related_dense_glmmkin_find_weight(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
    )
    assert results["num_variant"] > 0
    assert results["cMAC"] > 0
    assert 0.0 <= results["results_STAAR_O"] <= 1.0
    assert len(results["weight_all_1"]) > 0
    assert len(results["results_weight_staar_o"]) > 0


@pytest.mark.parametrize(
    "workflow",
    [
        workflows.ai_staar_related_sparse_glmmkin_find_weight,
        workflows.ai_staar_related_dense_glmmkin_find_weight,
    ],
)
def test_ai_staar_related_find_weight_workflow_does_not_load_precomputed_ai_cov_artifacts(
    monkeypatch,
    workflow,
):
    original_read_csv = workflows.pd.read_csv

    def _guard_ai_cov_reads(*args, **kwargs):
        path = str(args[0]) if args else ""
        if "example_ai_cov_" in path:
            raise AssertionError(
                "related AI find-weight workflow should not load precomputed AI covariance artifacts"
            )
        return original_read_csv(*args, **kwargs)

    monkeypatch.setattr(workflows.pd, "read_csv", _guard_ai_cov_reads)

    results = workflow(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        use_precomputed_artifacts=True,
    )
    assert results["num_variant"] == 163.0
    assert 0.0 <= results["results_STAAR_O"] <= 1.0
