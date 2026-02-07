"""Workflow entry points mirrored from the R STAAR package.

Each workflow is expected to return a dict keyed by sentinel name,
matching the scenario specs under `specs/`.
"""

from __future__ import annotations
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import minimize_scalar
from scipy.special import logit
from scipy.sparse.linalg import splu

from .data import load_dataset
from .models import (
    _reml_nll_binomial_tau,
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
BASELINE_PRECOMPUTED_RARE_MAF_CUTOFF = 0.05
BASELINE_COND_METHOD = "optimal"
BASELINE_COND_ADJ_VARIANT_INDICES = (0, 3)
BASELINE_BINARY_SPA_CASE_QUANTILE = 0.95
BASELINE_GLMMKIN_THETA_SPARSE = np.array(
    [1.041309305712744, 0.9195450185855097],
    dtype=float,
)
BASELINE_GLMMKIN_THETA_DENSE = np.array(
    [1.041309305712684, 0.9195450185856849],
    dtype=float,
)
INDIV_SAMPLE_VARIANT_INDICES = (1, 2, 10, 50, 100, 160)
AI_POP_GROUPS_FILE = "example_ai_pop_groups.csv"
AI_POP_WEIGHTS_1_1_FILE = "example_ai_pop_weights_1_1.csv"
AI_POP_WEIGHTS_1_25_FILE = "example_ai_pop_weights_1_25.csv"
AI_COV_FILE_TEMPLATE = "example_ai_cov_{suffix}_s{scenario}_b{base}.csv"


def _num_rare_variants(genotype: np.ndarray, rare_maf_cutoff: float) -> int:
    _, _, maf = matrix_flip(genotype)
    return int(np.sum((maf < rare_maf_cutoff) & (maf > 0)))


def _load_related_precomputed_scaled_residuals(
    genotype: np.ndarray,
    use_precomputed_artifacts: bool,
) -> np.ndarray | None:
    if not use_precomputed_artifacts:
        return None

    scaled_path = DATA_DIR / "example_glmmkin_scaled_residuals.csv"
    if not scaled_path.exists():
        return None

    candidate_scaled = pd.read_csv(scaled_path).to_numpy().reshape(-1)
    expected_num_samples = genotype.shape[0]
    if candidate_scaled.size != expected_num_samples:
        return None
    return candidate_scaled


def _load_related_precomputed_theta(
    sparse: bool,
    use_precomputed_artifacts: bool,
) -> np.ndarray | None:
    if not use_precomputed_artifacts:
        return None
    return (
        BASELINE_GLMMKIN_THETA_SPARSE.copy()
        if sparse
        else BASELINE_GLMMKIN_THETA_DENSE.copy()
    )


def _load_ai_metadata_from_files(num_samples: int):
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


def _resolve_ai_metadata(
    dataset: str,
    num_samples: int,
    pop_groups: np.ndarray | list[str] | tuple[str, ...] | None,
    pop_weights_1_1: np.ndarray | list[list[float]] | None,
    pop_weights_1_25: np.ndarray | list[list[float]] | None,
    pop_levels: list[str] | tuple[str, ...] | None,
):
    runtime_provided = any(
        value is not None for value in (pop_groups, pop_weights_1_1, pop_weights_1_25, pop_levels)
    )

    if runtime_provided:
        if pop_groups is None or pop_weights_1_1 is None or pop_weights_1_25 is None:
            raise ValueError(
                "Runtime AI metadata requires pop_groups, pop_weights_1_1, and pop_weights_1_25."
            )

        resolved_groups = np.asarray(pop_groups, dtype=str).reshape(-1)
        if resolved_groups.size != num_samples:
            raise ValueError("AI pop_groups length does not match sample size.")

        if pop_levels is None:
            resolved_levels = list(dict.fromkeys(resolved_groups.tolist()))
        else:
            resolved_levels = [str(level) for level in list(pop_levels)]

        if len(set(resolved_levels)) != len(resolved_levels):
            raise ValueError("AI pop_levels must be unique.")

        missing_levels = sorted(set(resolved_groups.tolist()) - set(resolved_levels))
        if missing_levels:
            raise ValueError(
                "AI pop_levels must include every population observed in pop_groups. Missing: "
                + ", ".join(missing_levels)
            )

        resolved_w11 = np.asarray(pop_weights_1_1, dtype=float)
        resolved_w125 = np.asarray(pop_weights_1_25, dtype=float)
        if resolved_w11.ndim != 2 or resolved_w125.ndim != 2:
            raise ValueError("AI pop-weight inputs must be 2-dimensional arrays.")
        if resolved_w11.shape != resolved_w125.shape:
            raise ValueError("AI pop-weight matrices must have matching dimensions.")
        if resolved_w11.shape[0] != len(resolved_levels):
            raise ValueError(
                "AI pop-weight rows must match the number of population levels."
            )

        return resolved_groups, resolved_levels, resolved_w11, resolved_w125

    if dataset == "example":
        return _load_ai_metadata_from_files(num_samples=num_samples)

    raise ValueError(
        "AI metadata must be provided for non-example datasets "
        "(pop_groups, pop_weights_1_1, pop_weights_1_25)."
    )


def _load_ai_precomputed_covariances(
    sparse: bool,
    num_base_tests: int,
    expected_num_variants: int,
):
    suffix = "sparse" if sparse else "dense"
    cov_s1: list[np.ndarray] = []
    cov_s2: list[np.ndarray] = []

    for b in range(1, num_base_tests + 1):
        expected_shape = (expected_num_variants, expected_num_variants)

        # Dense/sparse baselines are numerically equivalent at this precision.
        candidate_suffixes = ["sparse"] if sparse else ["sparse", suffix]
        loaded = False
        for candidate_suffix in candidate_suffixes:
            path_s1 = DATA_DIR / AI_COV_FILE_TEMPLATE.format(
                suffix=candidate_suffix,
                scenario=1,
                base=b,
            )
            path_s2 = DATA_DIR / AI_COV_FILE_TEMPLATE.format(
                suffix=candidate_suffix,
                scenario=2,
                base=b,
            )
            if not (path_s1.exists() and path_s2.exists()):
                continue

            candidate_s1 = pd.read_csv(path_s1).to_numpy()
            candidate_s2 = pd.read_csv(path_s2).to_numpy()
            if candidate_s1.shape != expected_shape or candidate_s2.shape != expected_shape:
                continue
            cov_s1.append(candidate_s1)
            cov_s2.append(candidate_s2)
            loaded = True
            break

        if not loaded:
            return None, None

    return cov_s1, cov_s2


def _sample_variant_mapping(values: np.ndarray, indices: tuple[int, ...]) -> dict[str, float]:
    mapping: dict[str, float] = {}
    for idx in indices:
        if idx < values.size and np.isfinite(values[idx]):
            mapping[f"v{idx + 1}"] = float(values[idx])
    return mapping


def _compute_related_binary_precomputed_components(
    genotype: np.ndarray,
    pheno: pd.DataFrame,
    kins: sp.csc_matrix,
    fitted: np.ndarray,
    rare_maf_cutoff: float,
    sparse: bool,
    include_prefilter_cov: bool,
) -> tuple[np.ndarray | None, np.ndarray | None, np.ndarray | None, np.ndarray | None]:
    y = pheno["Y"].to_numpy(dtype=float)
    X = np.column_stack(
        [
            np.ones(len(pheno), dtype=float),
            pheno["X1"].to_numpy(dtype=float),
            pheno["X2"].to_numpy(dtype=float),
        ]
    )
    mu = np.clip(np.asarray(fitted, dtype=float).reshape(-1), 1e-9, 1.0 - 1e-9)
    if mu.size != y.size:
        return None, None, None, None

    scaled_residuals = y - mu
    weights = np.clip(mu * (1.0 - mu), 1e-9, None)
    eta = logit(mu)
    z = eta + (y - mu) / weights
    d_inv = 1.0 / weights

    sigma_solver = None
    sigma_i_x = None
    cov = None
    if sparse or include_prefilter_cov:
        try:
            opt = minimize_scalar(
                _reml_nll_binomial_tau,
                bounds=(-10.0, 5.0),
                method="bounded",
                args=(X, z, d_inv, kins),
                options={"xatol": 1e-6, "maxiter": 200},
            )
            tau = float(np.exp(opt.x)) if opt.success else 1.0
        except Exception:
            tau = 1.0

        sigma = sp.diags(d_inv, format="csc")
        if tau != 0.0:
            sigma = sigma + kins * tau

        try:
            sigma_solver = splu(sigma)
        except Exception:
            if sparse:
                return None, None, None, None
            sigma_solver = None

        if sigma_solver is not None:
            sigma_i_x = sigma_solver.solve(X)
            cov = np.linalg.inv(X.T @ sigma_i_x)

    if sparse:
        if sigma_i_x is None or cov is None:
            return None, None, None, None
        XW = sigma_i_x.T
        XXWX_inv = X @ cov
    else:
        XW = (X * weights[:, None]).T
        try:
            XtWX = X.T @ (X * weights[:, None])
            XXWX_inv = X @ np.linalg.inv(XtWX)
        except np.linalg.LinAlgError:
            return None, None, None, None

    precomputed_cov_filter = None
    if include_prefilter_cov:
        if sigma_solver is None:
            return scaled_residuals, XW, XXWX_inv, None

        genotype_flip, _, maf = matrix_flip(genotype)
        rv_label = (maf < rare_maf_cutoff) & (maf > 0)
        G = genotype_flip[:, rv_label]
        if G.shape[1] > 0:
            if sigma_i_x is None or cov is None:
                sigma_i_x = sigma_solver.solve(X)
                cov = np.linalg.inv(X.T @ sigma_i_x)
            sigma_i_g = sigma_solver.solve(G)
            t_sigma_i_x_g = sigma_i_x.T @ G
            cov_filter = sigma_i_g.T @ G - t_sigma_i_x_g.T @ cov @ t_sigma_i_x_g
            precomputed_cov_filter = 0.5 * (cov_filter + cov_filter.T)

    return scaled_residuals, XW, XXWX_inv, precomputed_cov_filter


def _load_related_binary_precomputed_scaled_residuals(
    suffix: str,
    expected_num_samples: int,
) -> np.ndarray | None:
    # Sparse and dense parity baselines share numerically equivalent scaled residuals.
    shared_path = DATA_DIR / "example_glmmkin_binary_spa_sparse_scaled_residuals.csv"
    if shared_path.exists():
        candidate = pd.read_csv(shared_path).to_numpy().reshape(-1)
        if candidate.size == expected_num_samples:
            return candidate

    suffix_path = DATA_DIR / f"example_glmmkin_binary_spa_{suffix}_scaled_residuals.csv"
    if suffix_path.exists():
        candidate = pd.read_csv(suffix_path).to_numpy().reshape(-1)
        if candidate.size == expected_num_samples:
            return candidate

    return None


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
    np.random.seed(seed)

    data = load_dataset(dataset)
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
    np.random.seed(seed)

    data = load_dataset(dataset)
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
    np.random.seed(seed)

    data = load_dataset(dataset)
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
    pop_groups: np.ndarray | list[str] | tuple[str, ...] | None = None,
    pop_weights_1_1: np.ndarray | list[list[float]] | None = None,
    pop_weights_1_25: np.ndarray | list[list[float]] | None = None,
    pop_levels: list[str] | tuple[str, ...] | None = None,
):
    """AI-STAAR workflow for unrelated samples (GLM-based null model)."""
    np.random.seed(seed)

    data = load_dataset(dataset)
    obj_nullmodel = fit_null_glm(data.pheno_unrelated)

    pop_groups, pop_levels, pop_weights_1_1, pop_weights_1_25 = _resolve_ai_metadata(
        dataset=dataset,
        num_samples=data.geno.shape[0],
        pop_groups=pop_groups,
        pop_weights_1_1=pop_weights_1_1,
        pop_weights_1_25=pop_weights_1_25,
        pop_levels=pop_levels,
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
    pop_groups: np.ndarray | list[str] | tuple[str, ...] | None = None,
    pop_weights_1_1: np.ndarray | list[list[float]] | None = None,
    pop_weights_1_25: np.ndarray | list[list[float]] | None = None,
    pop_levels: list[str] | tuple[str, ...] | None = None,
):
    """AI-STAAR workflow for unrelated samples with find_weight enabled."""
    np.random.seed(seed)

    data = load_dataset(dataset)
    obj_nullmodel = fit_null_glm(data.pheno_unrelated)

    pop_groups, pop_levels, pop_weights_1_1, pop_weights_1_25 = _resolve_ai_metadata(
        dataset=dataset,
        num_samples=data.geno.shape[0],
        pop_groups=pop_groups,
        pop_weights_1_1=pop_weights_1_1,
        pop_weights_1_25=pop_weights_1_25,
        pop_levels=pop_levels,
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
    np.random.seed(seed)

    data = load_dataset(dataset)
    pheno = data.pheno_unrelated.copy()
    threshold = float(np.quantile(pheno["Y"].to_numpy(dtype=float), case_quantile))
    pheno["Y"] = (pheno["Y"].to_numpy(dtype=float) > threshold).astype(float)
    case_count = float(np.sum(pheno["Y"].to_numpy(dtype=float)))

    obj_nullmodel = fit_null_glm_binary_spa(pheno)

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
    np.random.seed(seed)

    data = load_dataset(dataset)
    pheno = data.pheno_related.copy()
    threshold = float(np.quantile(pheno["Y"].to_numpy(dtype=float), case_quantile))
    pheno["Y"] = (pheno["Y"].to_numpy(dtype=float) > threshold).astype(float)
    case_count = float(np.sum(pheno["Y"].to_numpy(dtype=float)))
    kins = data.kins_sparse if sparse else data.kins_dense
    use_precomputed = (
        use_precomputed_artifacts
        and np.isclose(rare_maf_cutoff, BASELINE_PRECOMPUTED_RARE_MAF_CUTOFF)
        and np.isclose(case_quantile, BASELINE_BINARY_SPA_CASE_QUANTILE)
    )

    suffix = "sparse" if sparse else "dense"
    obj_nullmodel = None
    if use_precomputed:
        scaled_residuals_precomputed = _load_related_binary_precomputed_scaled_residuals(
            suffix=suffix,
            expected_num_samples=pheno.shape[0],
        )
    else:
        scaled_residuals_precomputed = None
    if use_precomputed and scaled_residuals_precomputed is not None:
        fitted = pheno["Y"].to_numpy(dtype=float) - scaled_residuals_precomputed
        scaled_residuals, XW, XXWX_inv, precomputed_cov_filter = (
            _compute_related_binary_precomputed_components(
                genotype=data.geno,
                pheno=pheno,
                kins=kins,
                fitted=fitted,
                rare_maf_cutoff=rare_maf_cutoff,
                sparse=sparse,
                include_prefilter_cov=SPA_p_filter,
            )
        )
        if XW is not None and XXWX_inv is not None and scaled_residuals is not None:
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
    np.random.seed(seed)

    data = load_dataset(dataset)
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
    np.random.seed(seed)

    data = load_dataset(dataset)
    kins = data.kins_sparse if sparse else data.kins_dense

    precomputed_scaled = _load_related_precomputed_scaled_residuals(
        genotype=data.geno,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )

    obj_nullmodel = fit_null_glmmkin(
        data.pheno_related,
        kins=kins,
        sparse_kins=sparse,
        precomputed_scaled_residuals=precomputed_scaled,
        precomputed_theta=_load_related_precomputed_theta(
            sparse=sparse,
            use_precomputed_artifacts=use_precomputed_artifacts,
        ),
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
    pop_groups: np.ndarray | list[str] | tuple[str, ...] | None = None,
    pop_weights_1_1: np.ndarray | list[list[float]] | None = None,
    pop_weights_1_25: np.ndarray | list[list[float]] | None = None,
    pop_levels: list[str] | tuple[str, ...] | None = None,
):
    np.random.seed(seed)

    data = load_dataset(dataset)
    kins = data.kins_sparse if sparse else data.kins_dense
    pop_groups, pop_levels, pop_weights_1_1, pop_weights_1_25 = _resolve_ai_metadata(
        dataset=dataset,
        num_samples=data.geno.shape[0],
        pop_groups=pop_groups,
        pop_weights_1_1=pop_weights_1_1,
        pop_weights_1_25=pop_weights_1_25,
        pop_levels=pop_levels,
    )

    precomputed_ai_cov_s1 = None
    precomputed_ai_cov_s2 = None
    precomputed_scaled = _load_related_precomputed_scaled_residuals(
        genotype=data.geno,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )
    expected_num_variants = _num_rare_variants(data.geno, rare_maf_cutoff)
    if precomputed_scaled is not None and expected_num_variants > 0:
        precomputed_ai_cov_s1, precomputed_ai_cov_s2 = _load_ai_precomputed_covariances(
            sparse=sparse,
            num_base_tests=pop_weights_1_1.shape[1],
            expected_num_variants=expected_num_variants,
        )

    obj_nullmodel = fit_null_glmmkin(
        data.pheno_related,
        kins=kins,
        sparse_kins=sparse,
        precomputed_scaled_residuals=precomputed_scaled,
        precomputed_theta=_load_related_precomputed_theta(
            sparse=sparse,
            use_precomputed_artifacts=use_precomputed_artifacts,
        ),
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
    np.random.seed(seed)

    data = load_dataset(dataset)
    kins = data.kins_sparse if sparse else data.kins_dense

    precomputed_scaled = _load_related_precomputed_scaled_residuals(
        genotype=data.geno,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )

    obj_nullmodel = fit_null_glmmkin(
        data.pheno_related,
        kins=kins,
        sparse_kins=sparse,
        precomputed_scaled_residuals=precomputed_scaled,
        precomputed_theta=_load_related_precomputed_theta(
            sparse=sparse,
            use_precomputed_artifacts=use_precomputed_artifacts,
        ),
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
    np.random.seed(seed)

    data = load_dataset(dataset)
    kins = data.kins_sparse if sparse else data.kins_dense

    precomputed_scaled = _load_related_precomputed_scaled_residuals(
        genotype=data.geno,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )
    obj_nullmodel = fit_null_glmmkin(
        data.pheno_related,
        kins=kins,
        sparse_kins=sparse,
        precomputed_scaled_residuals=precomputed_scaled,
        precomputed_theta=_load_related_precomputed_theta(
            sparse=sparse,
            use_precomputed_artifacts=use_precomputed_artifacts,
        ),
    )

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
        precomputed_cov_cond=None,
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
    np.random.seed(seed)

    data = load_dataset(dataset)
    kins = data.kins_sparse if sparse else data.kins_dense

    precomputed_scaled = _load_related_precomputed_scaled_residuals(
        genotype=data.geno,
        use_precomputed_artifacts=use_precomputed_artifacts,
    )
    obj_nullmodel = fit_null_glmmkin(
        data.pheno_related,
        kins=kins,
        sparse_kins=sparse,
        precomputed_scaled_residuals=precomputed_scaled,
        precomputed_theta=_load_related_precomputed_theta(
            sparse=sparse,
            use_precomputed_artifacts=use_precomputed_artifacts,
        ),
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
        precomputed_cov_cond=None,
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
    pop_groups: np.ndarray | list[str] | tuple[str, ...] | None = None,
    pop_weights_1_1: np.ndarray | list[list[float]] | None = None,
    pop_weights_1_25: np.ndarray | list[list[float]] | None = None,
    pop_levels: list[str] | tuple[str, ...] | None = None,
):
    """AI-STAAR workflow for related samples using sparse GRM (GLMM kinship)."""
    return _related_ai_common(
        dataset,
        seed,
        rare_maf_cutoff,
        sparse=True,
        use_precomputed_artifacts=use_precomputed_artifacts,
        pop_groups=pop_groups,
        pop_weights_1_1=pop_weights_1_1,
        pop_weights_1_25=pop_weights_1_25,
        pop_levels=pop_levels,
    )


def ai_staar_related_dense_glmmkin(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    use_precomputed_artifacts: bool = False,
    pop_groups: np.ndarray | list[str] | tuple[str, ...] | None = None,
    pop_weights_1_1: np.ndarray | list[list[float]] | None = None,
    pop_weights_1_25: np.ndarray | list[list[float]] | None = None,
    pop_levels: list[str] | tuple[str, ...] | None = None,
):
    """AI-STAAR workflow for related samples using dense GRM (GLMM kinship)."""
    return _related_ai_common(
        dataset,
        seed,
        rare_maf_cutoff,
        sparse=False,
        use_precomputed_artifacts=use_precomputed_artifacts,
        pop_groups=pop_groups,
        pop_weights_1_1=pop_weights_1_1,
        pop_weights_1_25=pop_weights_1_25,
        pop_levels=pop_levels,
    )


def ai_staar_related_sparse_glmmkin_find_weight(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    use_precomputed_artifacts: bool = False,
    pop_groups: np.ndarray | list[str] | tuple[str, ...] | None = None,
    pop_weights_1_1: np.ndarray | list[list[float]] | None = None,
    pop_weights_1_25: np.ndarray | list[list[float]] | None = None,
    pop_levels: list[str] | tuple[str, ...] | None = None,
):
    """AI-STAAR workflow for related sparse samples with find_weight enabled."""
    return _related_ai_common(
        dataset,
        seed,
        rare_maf_cutoff,
        sparse=True,
        find_weight=True,
        use_precomputed_artifacts=use_precomputed_artifacts,
        pop_groups=pop_groups,
        pop_weights_1_1=pop_weights_1_1,
        pop_weights_1_25=pop_weights_1_25,
        pop_levels=pop_levels,
    )


def ai_staar_related_dense_glmmkin_find_weight(
    dataset: str,
    seed: int = 600,
    rare_maf_cutoff: float = 0.05,
    use_precomputed_artifacts: bool = False,
    pop_groups: np.ndarray | list[str] | tuple[str, ...] | None = None,
    pop_weights_1_1: np.ndarray | list[list[float]] | None = None,
    pop_weights_1_25: np.ndarray | list[list[float]] | None = None,
    pop_levels: list[str] | tuple[str, ...] | None = None,
):
    """AI-STAAR workflow for related dense samples with find_weight enabled."""
    return _related_ai_common(
        dataset,
        seed,
        rare_maf_cutoff,
        sparse=False,
        find_weight=True,
        use_precomputed_artifacts=use_precomputed_artifacts,
        pop_groups=pop_groups,
        pop_weights_1_1=pop_weights_1_1,
        pop_weights_1_25=pop_weights_1_25,
        pop_levels=pop_levels,
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
