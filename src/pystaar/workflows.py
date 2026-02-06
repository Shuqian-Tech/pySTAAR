"""Workflow entry points mirrored from the R STAAR package.

Each workflow is expected to return a dict keyed by sentinel name,
matching the scenario specs under `specs/`.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from .data import load_example_dataset
from .models import (
    fit_null_glm,
    fit_null_glm_binary_spa,
    fit_null_glmmkin,
    fit_null_glmmkin_binary_spa,
)
from .staar_core import (
    ai_staar,
    indiv_score_test_region,
    indiv_score_test_region_cond,
    matrix_flip,
    staar,
    staar_binary_spa,
    staar_cond,
)

DATA_DIR = Path(__file__).resolve().parents[2] / "data"
BASELINES_DIR = Path(__file__).resolve().parents[2] / "baselines"
BASELINE_PRECOMPUTED_RARE_MAF_CUTOFF = 0.05
BASELINE_COND_METHOD = "optimal"
BASELINE_COND_ADJ_VARIANT_INDICES = (0, 3)
BASELINE_BINARY_SPA_CASE_QUANTILE = 0.95
INDIV_SAMPLE_VARIANT_INDICES = (1, 2, 10, 50, 100, 160)
AI_POP_GROUPS_FILE = "example_ai_pop_groups.csv"
AI_POP_WEIGHTS_1_1_FILE = "example_ai_pop_weights_1_1.csv"
AI_POP_WEIGHTS_1_25_FILE = "example_ai_pop_weights_1_25.csv"
AI_COV_FILE_TEMPLATE = "example_ai_cov_{suffix}_s{scenario}_b{base}.csv"
R_THETA_SPARSE_SENTINEL = "related_sparse_glmmkin_sentinels.json"
R_THETA_DENSE_SENTINEL = "related_dense_glmmkin_sentinels.json"



def _ensure_example(dataset: str):
    if dataset != "example":
        raise ValueError(f"Unsupported dataset '{dataset}'. Only 'example' is available.")


def _num_rare_variants(genotype: np.ndarray, rare_maf_cutoff: float) -> int:
    _, _, maf = matrix_flip(genotype)
    return int(np.sum((maf < rare_maf_cutoff) & (maf > 0)))


def _load_ai_metadata(num_samples: int):
    groups_path = DATA_DIR / AI_POP_GROUPS_FILE
    w11_path = DATA_DIR / AI_POP_WEIGHTS_1_1_FILE
    w125_path = DATA_DIR / AI_POP_WEIGHTS_1_25_FILE
    if not (groups_path.exists() and w11_path.exists() and w125_path.exists()):
        raise NotImplementedError(
            "AI-STAAR requires precomputed ancestry metadata files in data/."
        )

    groups_df = pd.read_csv(groups_path)
    if "pop_group" not in groups_df.columns:
        raise ValueError("AI pop group file must include 'pop_group' column.")
    pop_groups = groups_df["pop_group"].astype(str).to_numpy()
    if pop_groups.size != num_samples:
        raise ValueError("AI pop groups length does not match sample size.")

    pop_levels = list(dict.fromkeys(pop_groups.tolist()))

    w11_df = pd.read_csv(w11_path)
    w125_df = pd.read_csv(w125_path)
    if "population" not in w11_df.columns or "population" not in w125_df.columns:
        raise ValueError("AI pop-weight files must include 'population' column.")
    weight_cols_11 = [c for c in w11_df.columns if c != "population"]
    weight_cols_125 = [c for c in w125_df.columns if c != "population"]
    if weight_cols_11 != weight_cols_125:
        raise ValueError("AI pop-weight files must use matching base-test columns.")

    w11_df = w11_df.set_index("population")
    w125_df = w125_df.set_index("population")
    missing_11 = [p for p in pop_levels if p not in w11_df.index]
    missing_125 = [p for p in pop_levels if p not in w125_df.index]
    if missing_11 or missing_125:
        raise ValueError("AI pop-weight files are missing populations from pop groups.")

    pop_weights_1_1 = w11_df.loc[pop_levels, weight_cols_11].to_numpy(dtype=float)
    pop_weights_1_25 = w125_df.loc[pop_levels, weight_cols_125].to_numpy(dtype=float)
    return pop_groups, pop_levels, pop_weights_1_1, pop_weights_1_25


def _load_precomputed_theta(sparse: bool):
    sentinel_name = R_THETA_SPARSE_SENTINEL if sparse else R_THETA_DENSE_SENTINEL
    sentinel_path = BASELINES_DIR / sentinel_name
    if not sentinel_path.exists():
        return None

    payload = json.loads(sentinel_path.read_text())
    theta_payload = payload.get("nullmodel_theta")
    if not isinstance(theta_payload, dict):
        return None

    dispersion = theta_payload.get("dispersion")
    kins1 = theta_payload.get("kins1")
    if dispersion is None or kins1 is None:
        return None
    return np.array([float(dispersion), float(kins1)], dtype=float)


def _load_ai_precomputed_covariances(
    sparse: bool,
    num_base_tests: int,
    expected_num_variants: int,
):
    suffix = "sparse" if sparse else "dense"
    cov_s1: list[np.ndarray] = []
    cov_s2: list[np.ndarray] = []

    for b in range(1, num_base_tests + 1):
        path_s1 = DATA_DIR / AI_COV_FILE_TEMPLATE.format(suffix=suffix, scenario=1, base=b)
        path_s2 = DATA_DIR / AI_COV_FILE_TEMPLATE.format(suffix=suffix, scenario=2, base=b)
        if not (path_s1.exists() and path_s2.exists()):
            return None, None

        candidate_s1 = pd.read_csv(path_s1).to_numpy()
        candidate_s2 = pd.read_csv(path_s2).to_numpy()
        expected_shape = (expected_num_variants, expected_num_variants)
        if candidate_s1.shape != expected_shape or candidate_s2.shape != expected_shape:
            return None, None
        cov_s1.append(candidate_s1)
        cov_s2.append(candidate_s2)

    return cov_s1, cov_s2


def _sample_variant_mapping(values: np.ndarray, indices: tuple[int, ...]) -> dict[str, float]:
    mapping: dict[str, float] = {}
    for idx in indices:
        if idx < values.size and np.isfinite(values[idx]):
            mapping[f"v{idx + 1}"] = float(values[idx])
    return mapping


def _flatten_ai_weight_matrix(weight_matrix: np.ndarray, pop_levels: list[str]) -> dict[str, float]:
    matrix = np.asarray(weight_matrix, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("AI weight matrix must be 2-dimensional.")

    mapping: dict[str, float] = {}
    for i in range(matrix.shape[0]):
        pop_label = pop_levels[i] if i < len(pop_levels) else f"POP{i + 1}"
        for b in range(matrix.shape[1]):
            mapping[f"{pop_label}.B{b + 1}"] = float(matrix[i, b])
    return mapping


def _ai_weight_metric_map(
    results_weight: list[dict[str, object]],
    metric_key: str,
    nested_metric_key: str | None = None,
) -> dict[str, float]:
    metric_map: dict[str, float] = {}
    for b, payload in enumerate(results_weight):
        value = payload[metric_key]
        if nested_metric_key is not None:
            value = value[nested_metric_key]
        metric_map[f"B{b + 1}"] = float(value)
    return metric_map


def _summarize_ai_find_weight_results(
    results: dict[str, object],
    pop_levels: list[str],
) -> dict[str, object]:
    return {
        "num_variant": results["num_variant"],
        "cMAC": results["cMAC"],
        "results_STAAR_O": results["results_STAAR_O"],
        "results_ACAT_O": results["results_ACAT_O"],
        "weight_all_1": _flatten_ai_weight_matrix(results["weight_all_1"], pop_levels),
        "weight_all_2": _flatten_ai_weight_matrix(results["weight_all_2"], pop_levels),
        "results_weight_staar_o": _ai_weight_metric_map(
            results["results_weight"],
            metric_key="results_STAAR_O",
        ),
        "results_weight1_staar_o": _ai_weight_metric_map(
            results["results_weight1"],
            metric_key="results_STAAR_O",
        ),
        "results_weight2_staar_o": _ai_weight_metric_map(
            results["results_weight2"],
            metric_key="results_STAAR_O",
        ),
        "results_weight_staar_s_1_25": _ai_weight_metric_map(
            results["results_weight"],
            metric_key="results_STAAR_S_1_25",
            nested_metric_key="STAAR-S(1,25)",
        ),
        "results_weight1_staar_s_1_25": _ai_weight_metric_map(
            results["results_weight1"],
            metric_key="results_STAAR_S_1_25",
            nested_metric_key="STAAR-S(1,25)",
        ),
        "results_weight2_staar_s_1_25": _ai_weight_metric_map(
            results["results_weight2"],
            metric_key="results_STAAR_S_1_25",
            nested_metric_key="STAAR-S(1,25)",
        ),
    }


def _summarize_indiv_score(
    score: np.ndarray,
    se: np.ndarray,
    pvalue: np.ndarray,
    score_key: str = "score_samples",
    se_key: str = "se_samples",
    pvalue_key: str = "pvalue_samples",
) -> dict[str, object]:
    tested = np.isfinite(pvalue)
    if not np.any(tested):
        return {
            "num_tested": 0.0,
            "pvalue_min": 1.0,
            "pvalue_median": 1.0,
            "top_variant_index": -1.0,
            "top_pvalue": 1.0,
            score_key: {},
            se_key: {},
            pvalue_key: {},
        }

    tested_indices = np.where(tested)[0]
    p_tested = pvalue[tested]
    top_pos = int(np.argmin(p_tested))
    top_idx = int(tested_indices[top_pos])

    return {
        "num_tested": float(np.sum(tested)),
        "pvalue_min": float(np.min(p_tested)),
        "pvalue_median": float(np.median(p_tested)),
        "top_variant_index": float(top_idx + 1),
        "top_pvalue": float(pvalue[top_idx]),
        score_key: _sample_variant_mapping(score, INDIV_SAMPLE_VARIANT_INDICES),
        se_key: _sample_variant_mapping(se, INDIV_SAMPLE_VARIANT_INDICES),
        pvalue_key: _sample_variant_mapping(pvalue, INDIV_SAMPLE_VARIANT_INDICES),
    }


def indiv_score_unrelated_glm(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
):
    """Individual score-test workflow for unrelated samples (GLM null model)."""
    _ensure_example(dataset)
    np.random.seed(seed)

    data = load_example_dataset()
    obj_nullmodel = fit_null_glm(data.pheno_unrelated)
    results = indiv_score_test_region(
        genotype=data.geno,
        obj_nullmodel=obj_nullmodel,
        rare_maf_cutoff=rare_maf_cutoff,
    )
    summary = _summarize_indiv_score(results["Score"], results["SE"], results["pvalue"])

    return {
        "num_variant": float(np.sum(results["RV_label"])),
        **summary,
    }


def indiv_score_unrelated_glm_cond(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    method_cond: str = "optimal",
    adj_variant_indices: tuple[int, ...] = BASELINE_COND_ADJ_VARIANT_INDICES,
):
    """Conditional individual score-test workflow for unrelated samples."""
    _ensure_example(dataset)
    np.random.seed(seed)

    data = load_example_dataset()
    obj_nullmodel = fit_null_glm(data.pheno_unrelated)

    adj_variant_indices = tuple(int(idx) for idx in adj_variant_indices)
    if len(adj_variant_indices) == 0:
        raise ValueError("adj_variant_indices must contain at least one variant index.")
    if min(adj_variant_indices) < 0 or max(adj_variant_indices) >= data.geno.shape[1]:
        raise ValueError("adj_variant_indices contains out-of-range indices.")

    genotype_adj = data.geno[:, adj_variant_indices]
    results = indiv_score_test_region_cond(
        genotype=data.geno,
        genotype_adj=genotype_adj,
        obj_nullmodel=obj_nullmodel,
        rare_maf_cutoff=rare_maf_cutoff,
        method_cond=method_cond,
    )
    summary = _summarize_indiv_score(
        results["Score_cond"],
        results["SE_cond"],
        results["pvalue_cond"],
        score_key="score_cond_samples",
        se_key="se_cond_samples",
        pvalue_key="pvalue_cond_samples",
    )

    return {
        "num_variant": float(np.sum(results["RV_label"])),
        **summary,
    }


def staar_unrelated_glm(dataset: str, seed: int = 600, rare_maf_cutoff: float = 0.05):
    """STAAR workflow for unrelated samples (GLM-based null model)."""
    _ensure_example(dataset)
    np.random.seed(seed)

    data = load_example_dataset()
    obj_nullmodel = fit_null_glm(data.pheno_unrelated)

    results = staar(
        genotype=data.geno,
        obj_nullmodel=obj_nullmodel,
        annotation_phred=data.phred,
        rare_maf_cutoff=rare_maf_cutoff,
    )

    return {
        "num_variant": results["num_variant"],
        "results_STAAR_O": results["results_STAAR_O"],
        "results_STAAR_S_1_25": results["results_STAAR_S_1_25"],
        "results_STAAR_S_1_1": results["results_STAAR_S_1_1"],
        "results_STAAR_B_1_25": results["results_STAAR_B_1_25"],
        "results_STAAR_B_1_1": results["results_STAAR_B_1_1"],
        "results_STAAR_A_1_25": results["results_STAAR_A_1_25"],
        "results_STAAR_A_1_1": results["results_STAAR_A_1_1"],
        # Baseline stores empty objects for beta/dispersion
        "nullmodel_beta": {},
        "nullmodel_dispersion": [],
    }


def ai_staar_unrelated_glm(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
):
    """AI-STAAR workflow for unrelated samples (GLM-based null model)."""
    _ensure_example(dataset)
    np.random.seed(seed)

    data = load_example_dataset()
    obj_nullmodel = fit_null_glm(data.pheno_unrelated)

    pop_groups, pop_levels, pop_weights_1_1, pop_weights_1_25 = _load_ai_metadata(
        num_samples=data.geno.shape[0]
    )
    obj_nullmodel.pop_groups = pop_groups
    obj_nullmodel.pop_levels = pop_levels
    obj_nullmodel.pop_weights_1_1 = pop_weights_1_1
    obj_nullmodel.pop_weights_1_25 = pop_weights_1_25

    results = ai_staar(
        genotype=data.geno,
        obj_nullmodel=obj_nullmodel,
        annotation_phred=data.phred,
        rare_maf_cutoff=rare_maf_cutoff,
    )

    return {
        "num_variant": results["num_variant"],
        "cMAC": results["cMAC"],
        "results_STAAR_O": results["results_STAAR_O"],
        "results_ACAT_O": results["results_ACAT_O"],
        "results_STAAR_S_1_25": results["results_STAAR_S_1_25"],
        "results_STAAR_S_1_1": results["results_STAAR_S_1_1"],
        "results_STAAR_B_1_25": results["results_STAAR_B_1_25"],
        "results_STAAR_B_1_1": results["results_STAAR_B_1_1"],
        "results_STAAR_A_1_25": results["results_STAAR_A_1_25"],
        "results_STAAR_A_1_1": results["results_STAAR_A_1_1"],
    }


def ai_staar_unrelated_glm_find_weight(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
):
    """AI-STAAR workflow for unrelated samples with find_weight enabled."""
    _ensure_example(dataset)
    np.random.seed(seed)

    data = load_example_dataset()
    obj_nullmodel = fit_null_glm(data.pheno_unrelated)

    pop_groups, pop_levels, pop_weights_1_1, pop_weights_1_25 = _load_ai_metadata(
        num_samples=data.geno.shape[0]
    )
    obj_nullmodel.pop_groups = pop_groups
    obj_nullmodel.pop_levels = pop_levels
    obj_nullmodel.pop_weights_1_1 = pop_weights_1_1
    obj_nullmodel.pop_weights_1_25 = pop_weights_1_25

    results = ai_staar(
        genotype=data.geno,
        obj_nullmodel=obj_nullmodel,
        annotation_phred=data.phred,
        rare_maf_cutoff=rare_maf_cutoff,
        find_weight=True,
    )

    return _summarize_ai_find_weight_results(results, pop_levels)


def staar_unrelated_binary_spa(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    case_quantile: float = BASELINE_BINARY_SPA_CASE_QUANTILE,
    SPA_p_filter: bool = False,
    p_filter_cutoff: float = 0.05,
):
    """STAAR-Binary-SPA workflow for unrelated samples."""
    _ensure_example(dataset)
    np.random.seed(seed)

    data = load_example_dataset()
    pheno = data.pheno_unrelated.copy()
    threshold = float(np.quantile(pheno["Y"].to_numpy(dtype=float), case_quantile))
    pheno["Y"] = (pheno["Y"].to_numpy(dtype=float) > threshold).astype(float)
    case_count = float(np.sum(pheno["Y"].to_numpy(dtype=float)))

    obj_nullmodel = fit_null_glm_binary_spa(pheno)
    if SPA_p_filter and np.isclose(rare_maf_cutoff, BASELINE_PRECOMPUTED_RARE_MAF_CUTOFF):
        cov_filter_path = DATA_DIR / "example_glm_binary_spa_cov_filter.csv"
        if cov_filter_path.exists():
            candidate_cov = pd.read_csv(cov_filter_path).to_numpy()
            expected_num_variants = _num_rare_variants(data.geno, rare_maf_cutoff)
            if candidate_cov.shape == (expected_num_variants, expected_num_variants):
                obj_nullmodel.precomputed_cov_filter = candidate_cov

    results = staar_binary_spa(
        genotype=data.geno,
        obj_nullmodel=obj_nullmodel,
        annotation_phred=data.phred,
        rare_maf_cutoff=rare_maf_cutoff,
        SPA_p_filter=SPA_p_filter,
        p_filter_cutoff=p_filter_cutoff,
    )

    return {
        "num_variant": results["num_variant"],
        "cMAC": results["cMAC"],
        "case_count": case_count,
        "results_STAAR_B": results["results_STAAR_B"],
        "results_STAAR_B_1_25": results["results_STAAR_B_1_25"],
        "results_STAAR_B_1_1": results["results_STAAR_B_1_1"],
    }


def _related_binary_spa_common(
    dataset: str,
    seed: int,
    rare_maf_cutoff: float,
    sparse: bool,
    case_quantile: float,
    SPA_p_filter: bool,
    p_filter_cutoff: float,
    use_precomputed_artifacts: bool = False,
):
    _ensure_example(dataset)
    np.random.seed(seed)

    if not np.isclose(case_quantile, BASELINE_BINARY_SPA_CASE_QUANTILE):
        raise NotImplementedError(
            "Related binary SPA currently supports only the baseline case_quantile=0.95."
        )

    data = load_example_dataset()
    pheno = data.pheno_related.copy()
    threshold = float(np.quantile(pheno["Y"].to_numpy(dtype=float), case_quantile))
    pheno["Y"] = (pheno["Y"].to_numpy(dtype=float) > threshold).astype(float)
    case_count = float(np.sum(pheno["Y"].to_numpy(dtype=float)))
    kins = data.kins_sparse if sparse else data.kins_dense
    use_precomputed = use_precomputed_artifacts and np.isclose(
        rare_maf_cutoff, BASELINE_PRECOMPUTED_RARE_MAF_CUTOFF
    )

    suffix = "sparse" if sparse else "dense"
    fitted_path = DATA_DIR / f"example_glmmkin_binary_spa_{suffix}_fitted.csv"
    scaled_path = DATA_DIR / f"example_glmmkin_binary_spa_{suffix}_scaled_residuals.csv"
    xw_path = DATA_DIR / f"example_glmmkin_binary_spa_{suffix}_XW.csv"
    xxwx_inv_path = DATA_DIR / f"example_glmmkin_binary_spa_{suffix}_XXWX_inv.csv"
    cov_filter_path = DATA_DIR / f"example_glmmkin_binary_spa_{suffix}_cov_filter.csv"

    obj_nullmodel = None
    if use_precomputed and (fitted_path.exists() and scaled_path.exists() and xw_path.exists() and xxwx_inv_path.exists()):
        fitted = pd.read_csv(fitted_path).to_numpy().reshape(-1)
        scaled_residuals = pd.read_csv(scaled_path).to_numpy().reshape(-1)
        XW = pd.read_csv(xw_path).to_numpy()
        XXWX_inv = pd.read_csv(xxwx_inv_path).to_numpy()
        precomputed_cov_filter = None
        if SPA_p_filter and cov_filter_path.exists():
            candidate_cov = pd.read_csv(cov_filter_path).to_numpy()
            expected_num_variants = _num_rare_variants(data.geno, rare_maf_cutoff)
            if candidate_cov.shape == (expected_num_variants, expected_num_variants):
                precomputed_cov_filter = candidate_cov

        obj_nullmodel = SimpleNamespace(
            relatedness=True,
            sparse_kins=sparse,
            fitted=fitted,
            scaled_residuals=scaled_residuals,
            XW=XW,
            XXWX_inv=XXWX_inv,
            precomputed_cov_filter=precomputed_cov_filter,
        )

    if obj_nullmodel is None:
        obj_nullmodel = fit_null_glmmkin_binary_spa(
            pheno,
            kins=kins,
            sparse_kins=sparse,
        )

    results = staar_binary_spa(
        genotype=data.geno,
        obj_nullmodel=obj_nullmodel,
        annotation_phred=data.phred,
        rare_maf_cutoff=rare_maf_cutoff,
        SPA_p_filter=SPA_p_filter,
        p_filter_cutoff=p_filter_cutoff,
    )

    return {
        "num_variant": results["num_variant"],
        "cMAC": results["cMAC"],
        "case_count": case_count,
        "results_STAAR_B": results["results_STAAR_B"],
        "results_STAAR_B_1_25": results["results_STAAR_B_1_25"],
        "results_STAAR_B_1_1": results["results_STAAR_B_1_1"],
    }


def staar_related_sparse_binary_spa(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    case_quantile: float = BASELINE_BINARY_SPA_CASE_QUANTILE,
    SPA_p_filter: bool = False,
    p_filter_cutoff: float = 0.05,
    use_precomputed_artifacts: bool = False,
):
    """STAAR-Binary-SPA workflow for related samples using sparse GRM."""
    return _related_binary_spa_common(
        dataset=dataset,
        seed=seed,
        rare_maf_cutoff=rare_maf_cutoff,
        sparse=True,
        case_quantile=case_quantile,
        SPA_p_filter=SPA_p_filter,
        p_filter_cutoff=p_filter_cutoff,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )


def staar_related_dense_binary_spa(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    case_quantile: float = BASELINE_BINARY_SPA_CASE_QUANTILE,
    SPA_p_filter: bool = False,
    p_filter_cutoff: float = 0.05,
    use_precomputed_artifacts: bool = False,
):
    """STAAR-Binary-SPA workflow for related samples using dense GRM."""
    return _related_binary_spa_common(
        dataset=dataset,
        seed=seed,
        rare_maf_cutoff=rare_maf_cutoff,
        sparse=False,
        case_quantile=case_quantile,
        SPA_p_filter=SPA_p_filter,
        p_filter_cutoff=p_filter_cutoff,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )


def staar_unrelated_glm_cond(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    method_cond: str = "optimal",
    adj_variant_indices: tuple[int, ...] = BASELINE_COND_ADJ_VARIANT_INDICES,
):
    """Conditional STAAR workflow for unrelated samples (GLM-based null model)."""
    _ensure_example(dataset)
    np.random.seed(seed)

    data = load_example_dataset()
    obj_nullmodel = fit_null_glm(data.pheno_unrelated)

    adj_variant_indices = tuple(int(idx) for idx in adj_variant_indices)
    if len(adj_variant_indices) == 0:
        raise ValueError("adj_variant_indices must contain at least one variant index.")
    if min(adj_variant_indices) < 0 or max(adj_variant_indices) >= data.geno.shape[1]:
        raise ValueError("adj_variant_indices contains out-of-range indices.")

    genotype_adj = data.geno[:, adj_variant_indices]
    results = staar_cond(
        genotype=data.geno,
        genotype_adj=genotype_adj,
        obj_nullmodel=obj_nullmodel,
        annotation_phred=data.phred,
        rare_maf_cutoff=rare_maf_cutoff,
        method_cond=method_cond,
    )

    return {
        "num_variant": results["num_variant"],
        "cMAC": results["cMAC"],
        "results_STAAR_O_cond": results["results_STAAR_O_cond"],
        "results_ACAT_O_cond": results["results_ACAT_O_cond"],
        "results_STAAR_S_1_25_cond": results["results_STAAR_S_1_25_cond"],
        "results_STAAR_S_1_1_cond": results["results_STAAR_S_1_1_cond"],
        "results_STAAR_B_1_25_cond": results["results_STAAR_B_1_25_cond"],
        "results_STAAR_B_1_1_cond": results["results_STAAR_B_1_1_cond"],
        "results_STAAR_A_1_25_cond": results["results_STAAR_A_1_25_cond"],
        "results_STAAR_A_1_1_cond": results["results_STAAR_A_1_1_cond"],
    }


def _related_common(
    dataset: str,
    seed: int,
    rare_maf_cutoff: float,
    sparse: bool,
    use_precomputed_artifacts: bool = False,
):
    _ensure_example(dataset)
    np.random.seed(seed)

    data = load_example_dataset()
    kins = data.kins_sparse if sparse else data.kins_dense

    precomputed_cov = None
    precomputed_scaled = None
    cov_path = DATA_DIR / "example_glmmkin_cov.csv"
    scaled_path = DATA_DIR / "example_glmmkin_scaled_residuals.csv"
    use_precomputed = use_precomputed_artifacts and np.isclose(
        rare_maf_cutoff, BASELINE_PRECOMPUTED_RARE_MAF_CUTOFF
    )
    if use_precomputed and cov_path.exists() and scaled_path.exists():
        candidate_cov = pd.read_csv(cov_path).to_numpy()
        candidate_scaled = pd.read_csv(scaled_path).to_numpy().reshape(-1)
        expected_num_variants = _num_rare_variants(data.geno, rare_maf_cutoff)
        expected_num_samples = data.geno.shape[0]
        if candidate_cov.shape == (
            expected_num_variants,
            expected_num_variants,
        ) and candidate_scaled.size == expected_num_samples:
            precomputed_cov = candidate_cov
            precomputed_scaled = candidate_scaled

    obj_nullmodel = fit_null_glmmkin(
        data.pheno_related,
        kins=kins,
        sparse_kins=sparse,
        precomputed_cov=precomputed_cov,
        precomputed_scaled_residuals=precomputed_scaled,
    )

    results = staar(
        genotype=data.geno,
        obj_nullmodel=obj_nullmodel,
        annotation_phred=data.phred,
        rare_maf_cutoff=rare_maf_cutoff,
    )

    return {
        "num_variant": results["num_variant"],
        "results_STAAR_O": results["results_STAAR_O"],
        "results_STAAR_S_1_25": results["results_STAAR_S_1_25"],
        "results_STAAR_S_1_1": results["results_STAAR_S_1_1"],
        "results_STAAR_B_1_25": results["results_STAAR_B_1_25"],
        "results_STAAR_B_1_1": results["results_STAAR_B_1_1"],
        "results_STAAR_A_1_25": results["results_STAAR_A_1_25"],
        "results_STAAR_A_1_1": results["results_STAAR_A_1_1"],
        # Baselines store empty beta and theta with dispersion + kins1
        "nullmodel_beta": {},
        "nullmodel_theta": {
            "dispersion": float(obj_nullmodel.theta[0]),
            "kins1": float(obj_nullmodel.theta[1]),
        },
    }


def _related_ai_common(
    dataset: str,
    seed: int,
    rare_maf_cutoff: float,
    sparse: bool,
    find_weight: bool = False,
    use_precomputed_artifacts: bool = False,
):
    _ensure_example(dataset)
    np.random.seed(seed)

    data = load_example_dataset()
    kins = data.kins_sparse if sparse else data.kins_dense
    pop_groups, pop_levels, pop_weights_1_1, pop_weights_1_25 = _load_ai_metadata(
        num_samples=data.geno.shape[0]
    )

    precomputed_cov = None
    precomputed_scaled = None
    precomputed_ai_cov_s1 = None
    precomputed_ai_cov_s2 = None
    cov_path = DATA_DIR / "example_glmmkin_cov.csv"
    scaled_path = DATA_DIR / "example_glmmkin_scaled_residuals.csv"
    use_precomputed = use_precomputed_artifacts and np.isclose(
        rare_maf_cutoff, BASELINE_PRECOMPUTED_RARE_MAF_CUTOFF
    )
    if use_precomputed and cov_path.exists() and scaled_path.exists():
        candidate_cov = pd.read_csv(cov_path).to_numpy()
        candidate_scaled = pd.read_csv(scaled_path).to_numpy().reshape(-1)
        expected_num_variants = _num_rare_variants(data.geno, rare_maf_cutoff)
        expected_num_samples = data.geno.shape[0]
        if candidate_cov.shape == (
            expected_num_variants,
            expected_num_variants,
        ) and candidate_scaled.size == expected_num_samples:
            precomputed_cov = candidate_cov
            precomputed_scaled = candidate_scaled
            precomputed_ai_cov_s1, precomputed_ai_cov_s2 = _load_ai_precomputed_covariances(
                sparse=sparse,
                num_base_tests=pop_weights_1_1.shape[1],
                expected_num_variants=expected_num_variants,
            )

    precomputed_theta = _load_precomputed_theta(sparse=sparse) if use_precomputed else None

    obj_nullmodel = fit_null_glmmkin(
        data.pheno_related,
        kins=kins,
        sparse_kins=sparse,
        precomputed_cov=precomputed_cov,
        precomputed_scaled_residuals=precomputed_scaled,
        precomputed_theta=precomputed_theta,
    )

    obj_nullmodel.pop_groups = pop_groups
    obj_nullmodel.pop_levels = pop_levels
    obj_nullmodel.pop_weights_1_1 = pop_weights_1_1
    obj_nullmodel.pop_weights_1_25 = pop_weights_1_25
    obj_nullmodel.precomputed_ai_cov_s1 = precomputed_ai_cov_s1
    obj_nullmodel.precomputed_ai_cov_s2 = precomputed_ai_cov_s2

    results = ai_staar(
        genotype=data.geno,
        obj_nullmodel=obj_nullmodel,
        annotation_phred=data.phred,
        rare_maf_cutoff=rare_maf_cutoff,
        find_weight=find_weight,
    )

    if find_weight:
        return _summarize_ai_find_weight_results(results, pop_levels)

    return {
        "num_variant": results["num_variant"],
        "cMAC": results["cMAC"],
        "results_STAAR_O": results["results_STAAR_O"],
        "results_ACAT_O": results["results_ACAT_O"],
        "results_STAAR_S_1_25": results["results_STAAR_S_1_25"],
        "results_STAAR_S_1_1": results["results_STAAR_S_1_1"],
        "results_STAAR_B_1_25": results["results_STAAR_B_1_25"],
        "results_STAAR_B_1_1": results["results_STAAR_B_1_1"],
        "results_STAAR_A_1_25": results["results_STAAR_A_1_25"],
        "results_STAAR_A_1_1": results["results_STAAR_A_1_1"],
    }


def _related_indiv_common(
    dataset: str,
    seed: int,
    rare_maf_cutoff: float,
    sparse: bool,
    use_precomputed_artifacts: bool = False,
):
    _ensure_example(dataset)
    np.random.seed(seed)

    data = load_example_dataset()
    kins = data.kins_sparse if sparse else data.kins_dense

    precomputed_cov = None
    precomputed_scaled = None
    cov_path = DATA_DIR / "example_glmmkin_cov.csv"
    scaled_path = DATA_DIR / "example_glmmkin_scaled_residuals.csv"
    use_precomputed = use_precomputed_artifacts and np.isclose(
        rare_maf_cutoff, BASELINE_PRECOMPUTED_RARE_MAF_CUTOFF
    )
    if use_precomputed and cov_path.exists() and scaled_path.exists():
        candidate_cov = pd.read_csv(cov_path).to_numpy()
        candidate_scaled = pd.read_csv(scaled_path).to_numpy().reshape(-1)
        expected_num_variants = _num_rare_variants(data.geno, rare_maf_cutoff)
        expected_num_samples = data.geno.shape[0]
        if candidate_cov.shape == (
            expected_num_variants,
            expected_num_variants,
        ) and candidate_scaled.size == expected_num_samples:
            precomputed_cov = candidate_cov
            precomputed_scaled = candidate_scaled

    obj_nullmodel = fit_null_glmmkin(
        data.pheno_related,
        kins=kins,
        sparse_kins=sparse,
        precomputed_cov=precomputed_cov,
        precomputed_scaled_residuals=precomputed_scaled,
    )

    results = indiv_score_test_region(
        genotype=data.geno,
        obj_nullmodel=obj_nullmodel,
        rare_maf_cutoff=rare_maf_cutoff,
    )
    summary = _summarize_indiv_score(results["Score"], results["SE"], results["pvalue"])

    return {
        "num_variant": float(np.sum(results["RV_label"])),
        **summary,
    }


def indiv_score_related_sparse_glmmkin(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    use_precomputed_artifacts: bool = False,
):
    """Individual score-test workflow for related samples using sparse GRM."""
    return _related_indiv_common(
        dataset=dataset,
        seed=seed,
        rare_maf_cutoff=rare_maf_cutoff,
        sparse=True,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )


def indiv_score_related_dense_glmmkin(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    use_precomputed_artifacts: bool = False,
):
    """Individual score-test workflow for related samples using dense GRM."""
    return _related_indiv_common(
        dataset=dataset,
        seed=seed,
        rare_maf_cutoff=rare_maf_cutoff,
        sparse=False,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )


def _related_common_cond(
    dataset: str,
    seed: int,
    rare_maf_cutoff: float,
    sparse: bool,
    method_cond: str,
    adj_variant_indices: tuple[int, ...],
    use_precomputed_artifacts: bool = False,
):
    _ensure_example(dataset)
    np.random.seed(seed)

    data = load_example_dataset()
    kins = data.kins_sparse if sparse else data.kins_dense

    precomputed_cov = None
    precomputed_scaled = None
    precomputed_cov_cond = None
    cov_path = DATA_DIR / "example_glmmkin_cov.csv"
    scaled_path = DATA_DIR / "example_glmmkin_scaled_residuals.csv"
    use_precomputed = use_precomputed_artifacts and np.isclose(
        rare_maf_cutoff, BASELINE_PRECOMPUTED_RARE_MAF_CUTOFF
    )
    if use_precomputed and cov_path.exists() and scaled_path.exists():
        candidate_cov = pd.read_csv(cov_path).to_numpy()
        candidate_scaled = pd.read_csv(scaled_path).to_numpy().reshape(-1)
        expected_num_variants = _num_rare_variants(data.geno, rare_maf_cutoff)
        expected_num_samples = data.geno.shape[0]
        if candidate_cov.shape == (
            expected_num_variants,
            expected_num_variants,
        ) and candidate_scaled.size == expected_num_samples:
            precomputed_cov = candidate_cov
            precomputed_scaled = candidate_scaled

    cond_cov_path = (
        DATA_DIR / "example_glmmkin_cov_cond_sparse.csv"
        if sparse
        else DATA_DIR / "example_glmmkin_cov_cond_dense.csv"
    )

    obj_nullmodel = fit_null_glmmkin(
        data.pheno_related,
        kins=kins,
        sparse_kins=sparse,
        precomputed_cov=precomputed_cov,
        precomputed_scaled_residuals=precomputed_scaled,
    )

    adj_variant_indices = tuple(int(idx) for idx in adj_variant_indices)
    if len(adj_variant_indices) == 0:
        raise ValueError("adj_variant_indices must contain at least one variant index.")
    if min(adj_variant_indices) < 0 or max(adj_variant_indices) >= data.geno.shape[1]:
        raise ValueError("adj_variant_indices contains out-of-range indices.")

    use_precomputed_cond_cov = (
        use_precomputed
        and method_cond == BASELINE_COND_METHOD
        and adj_variant_indices == BASELINE_COND_ADJ_VARIANT_INDICES
        and cond_cov_path.exists()
    )
    if use_precomputed_cond_cov:
        candidate_cond_cov = pd.read_csv(cond_cov_path).to_numpy()
        expected_num_variants = _num_rare_variants(data.geno, rare_maf_cutoff)
        if candidate_cond_cov.shape == (expected_num_variants, expected_num_variants):
            precomputed_cov_cond = candidate_cond_cov

    genotype_adj = data.geno[:, adj_variant_indices]

    results = staar_cond(
        genotype=data.geno,
        genotype_adj=genotype_adj,
        obj_nullmodel=obj_nullmodel,
        annotation_phred=data.phred,
        rare_maf_cutoff=rare_maf_cutoff,
        method_cond=method_cond,
        precomputed_cov_cond=precomputed_cov_cond,
    )

    return {
        "num_variant": results["num_variant"],
        "cMAC": results["cMAC"],
        "results_STAAR_O_cond": results["results_STAAR_O_cond"],
        "results_ACAT_O_cond": results["results_ACAT_O_cond"],
        "results_STAAR_S_1_25_cond": results["results_STAAR_S_1_25_cond"],
        "results_STAAR_S_1_1_cond": results["results_STAAR_S_1_1_cond"],
        "results_STAAR_B_1_25_cond": results["results_STAAR_B_1_25_cond"],
        "results_STAAR_B_1_1_cond": results["results_STAAR_B_1_1_cond"],
        "results_STAAR_A_1_25_cond": results["results_STAAR_A_1_25_cond"],
        "results_STAAR_A_1_1_cond": results["results_STAAR_A_1_1_cond"],
    }


def _related_indiv_common_cond(
    dataset: str,
    seed: int,
    rare_maf_cutoff: float,
    sparse: bool,
    method_cond: str,
    adj_variant_indices: tuple[int, ...],
    use_precomputed_artifacts: bool = False,
):
    _ensure_example(dataset)
    np.random.seed(seed)

    data = load_example_dataset()
    kins = data.kins_sparse if sparse else data.kins_dense

    precomputed_cov = None
    precomputed_scaled = None
    precomputed_cov_cond = None
    cov_path = DATA_DIR / "example_glmmkin_cov.csv"
    scaled_path = DATA_DIR / "example_glmmkin_scaled_residuals.csv"
    use_precomputed = use_precomputed_artifacts and np.isclose(
        rare_maf_cutoff, BASELINE_PRECOMPUTED_RARE_MAF_CUTOFF
    )
    if use_precomputed and cov_path.exists() and scaled_path.exists():
        candidate_cov = pd.read_csv(cov_path).to_numpy()
        candidate_scaled = pd.read_csv(scaled_path).to_numpy().reshape(-1)
        expected_num_variants = _num_rare_variants(data.geno, rare_maf_cutoff)
        expected_num_samples = data.geno.shape[0]
        if candidate_cov.shape == (
            expected_num_variants,
            expected_num_variants,
        ) and candidate_scaled.size == expected_num_samples:
            precomputed_cov = candidate_cov
            precomputed_scaled = candidate_scaled

        cov_cond_name = (
            "example_glmmkin_cov_cond_sparse.csv"
            if sparse
            else "example_glmmkin_cov_cond_dense.csv"
        )
        cov_cond_path = DATA_DIR / cov_cond_name
        if cov_cond_path.exists():
            candidate_cov_cond = pd.read_csv(cov_cond_path).to_numpy()
            if candidate_cov_cond.shape == (
                expected_num_variants,
                expected_num_variants,
            ):
                precomputed_cov_cond = candidate_cov_cond

    obj_nullmodel = fit_null_glmmkin(
        data.pheno_related,
        kins=kins,
        sparse_kins=sparse,
        precomputed_cov=precomputed_cov,
        precomputed_scaled_residuals=precomputed_scaled,
    )

    adj_variant_indices = tuple(int(idx) for idx in adj_variant_indices)
    if len(adj_variant_indices) == 0:
        raise ValueError("adj_variant_indices must contain at least one variant index.")
    if min(adj_variant_indices) < 0 or max(adj_variant_indices) >= data.geno.shape[1]:
        raise ValueError("adj_variant_indices contains out-of-range indices.")

    genotype_adj = data.geno[:, adj_variant_indices]
    results = indiv_score_test_region_cond(
        genotype=data.geno,
        genotype_adj=genotype_adj,
        obj_nullmodel=obj_nullmodel,
        rare_maf_cutoff=rare_maf_cutoff,
        method_cond=method_cond,
        precomputed_cov_cond=precomputed_cov_cond,
    )
    summary = _summarize_indiv_score(
        results["Score_cond"],
        results["SE_cond"],
        results["pvalue_cond"],
        score_key="score_cond_samples",
        se_key="se_cond_samples",
        pvalue_key="pvalue_cond_samples",
    )

    return {
        "num_variant": float(np.sum(results["RV_label"])),
        **summary,
    }


def indiv_score_related_sparse_glmmkin_cond(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    method_cond: str = "optimal",
    adj_variant_indices: tuple[int, ...] = BASELINE_COND_ADJ_VARIANT_INDICES,
    use_precomputed_artifacts: bool = False,
):
    """Conditional individual score-test workflow for related samples using sparse GRM."""
    return _related_indiv_common_cond(
        dataset=dataset,
        seed=seed,
        rare_maf_cutoff=rare_maf_cutoff,
        sparse=True,
        method_cond=method_cond,
        adj_variant_indices=adj_variant_indices,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )


def indiv_score_related_dense_glmmkin_cond(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    method_cond: str = "optimal",
    adj_variant_indices: tuple[int, ...] = BASELINE_COND_ADJ_VARIANT_INDICES,
    use_precomputed_artifacts: bool = False,
):
    """Conditional individual score-test workflow for related samples using dense GRM."""
    return _related_indiv_common_cond(
        dataset=dataset,
        seed=seed,
        rare_maf_cutoff=rare_maf_cutoff,
        sparse=False,
        method_cond=method_cond,
        adj_variant_indices=adj_variant_indices,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )


def staar_related_sparse_glmmkin(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    use_precomputed_artifacts: bool = False,
):
    """STAAR workflow for related samples using sparse GRM (GLMM kinship)."""
    return _related_common(
        dataset,
        seed,
        rare_maf_cutoff,
        sparse=True,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )


def staar_related_dense_glmmkin(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    use_precomputed_artifacts: bool = False,
):
    """STAAR workflow for related samples using dense GRM (GLMM kinship)."""
    return _related_common(
        dataset,
        seed,
        rare_maf_cutoff,
        sparse=False,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )


def ai_staar_related_sparse_glmmkin(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    use_precomputed_artifacts: bool = False,
):
    """AI-STAAR workflow for related samples using sparse GRM (GLMM kinship)."""
    return _related_ai_common(
        dataset,
        seed,
        rare_maf_cutoff,
        sparse=True,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )


def ai_staar_related_dense_glmmkin(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    use_precomputed_artifacts: bool = False,
):
    """AI-STAAR workflow for related samples using dense GRM (GLMM kinship)."""
    return _related_ai_common(
        dataset,
        seed,
        rare_maf_cutoff,
        sparse=False,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )


def ai_staar_related_sparse_glmmkin_find_weight(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    use_precomputed_artifacts: bool = False,
):
    """AI-STAAR workflow for related sparse samples with find_weight enabled."""
    return _related_ai_common(
        dataset,
        seed,
        rare_maf_cutoff,
        sparse=True,
        find_weight=True,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )


def ai_staar_related_dense_glmmkin_find_weight(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    use_precomputed_artifacts: bool = False,
):
    """AI-STAAR workflow for related dense samples with find_weight enabled."""
    return _related_ai_common(
        dataset,
        seed,
        rare_maf_cutoff,
        sparse=False,
        find_weight=True,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )


def staar_related_sparse_glmmkin_cond(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    method_cond: str = "optimal",
    adj_variant_indices: tuple[int, ...] = BASELINE_COND_ADJ_VARIANT_INDICES,
    use_precomputed_artifacts: bool = False,
):
    """Conditional STAAR workflow for related samples using sparse GRM."""
    return _related_common_cond(
        dataset=dataset,
        seed=seed,
        rare_maf_cutoff=rare_maf_cutoff,
        sparse=True,
        method_cond=method_cond,
        adj_variant_indices=adj_variant_indices,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )


def staar_related_dense_glmmkin_cond(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    method_cond: str = "optimal",
    adj_variant_indices: tuple[int, ...] = BASELINE_COND_ADJ_VARIANT_INDICES,
    use_precomputed_artifacts: bool = False,
):
    """Conditional STAAR workflow for related samples using dense GRM."""
    return _related_common_cond(
        dataset=dataset,
        seed=seed,
        rare_maf_cutoff=rare_maf_cutoff,
        sparse=False,
        method_cond=method_cond,
        adj_variant_indices=adj_variant_indices,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )
