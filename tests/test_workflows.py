import math

import pytest

from pystaar import workflows


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


def test_related_workflow_uses_precomputed_cov_for_baseline_cutoff(monkeypatch):
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
    assert captured["precomputed_cov"] is not None
    assert captured["precomputed_cov"].shape == (163, 163)
    assert captured["precomputed_scaled_residuals"] is not None


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

    assert abs(pure["results_STAAR_B"] - precomputed["results_STAAR_B"]) < 1e-3


def test_unrelated_binary_spa_prefilter_loads_precomputed_cov(monkeypatch):
    original_staar_binary_spa = workflows.staar_binary_spa
    captured = {}

    def _track_staar_binary_spa(*args, **kwargs):
        obj = kwargs.get("obj_nullmodel")
        captured["has_precomputed_cov_filter"] = (
            obj is not None and getattr(obj, "precomputed_cov_filter", None) is not None
        )
        return original_staar_binary_spa(*args, **kwargs)

    monkeypatch.setattr(workflows, "staar_binary_spa", _track_staar_binary_spa)

    workflows.staar_unrelated_binary_spa(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
        SPA_p_filter=True,
        p_filter_cutoff=0.05,
    )

    assert captured["has_precomputed_cov_filter"] is True


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


def test_ai_staar_unrelated_workflow_runs():
    results = workflows.ai_staar_unrelated_glm(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
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


def test_ai_staar_related_dense_workflow_runs():
    results = workflows.ai_staar_related_dense_glmmkin(
        dataset="example",
        seed=600,
        rare_maf_cutoff=0.05,
    )
    assert results["num_variant"] > 0
    assert results["cMAC"] > 0
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
