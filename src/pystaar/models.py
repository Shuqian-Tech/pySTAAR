"""Null model fitting utilities for STAAR migration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.optimize import minimize, minimize_scalar
from scipy.sparse.linalg import splu
from scipy.special import expit


@dataclass
class NullModelGLM:
    X: np.ndarray
    y: np.ndarray
    fitted: np.ndarray
    residuals: np.ndarray
    beta: np.ndarray
    dispersion: float
    weights: np.ndarray
    family: str
    relatedness: bool = False
    use_SPA: bool = False
    XW: np.ndarray | None = None
    XXWX_inv: np.ndarray | None = None


@dataclass
class NullModelGLMMKin:
    X: np.ndarray
    y: np.ndarray
    fitted: np.ndarray
    residuals: np.ndarray
    scaled_residuals: np.ndarray
    beta: np.ndarray
    theta: np.ndarray
    cov: np.ndarray
    sigma_solver: splu
    Sigma_iX: np.ndarray
    sparse_kins: bool
    relatedness: bool = True
    precomputed_cov: np.ndarray | None = None


@dataclass
class NullModelGLMMKinBinarySPA:
    X: np.ndarray
    y: np.ndarray
    fitted: np.ndarray
    residuals: np.ndarray
    scaled_residuals: np.ndarray
    beta: np.ndarray
    theta: np.ndarray
    cov: np.ndarray
    sigma_solver: splu
    Sigma_iX: np.ndarray
    sparse_kins: bool
    weights: np.ndarray
    XW: np.ndarray
    XXWX_inv: np.ndarray
    relatedness: bool = True
    use_SPA: bool = True
    family: str = "binomial"
    precomputed_cov_filter: np.ndarray | None = None


def _normalize_covariate_cols(
    covariate_cols: str | Iterable[str] | None,
) -> tuple[str, ...]:
    if covariate_cols is None:
        return ()
    if isinstance(covariate_cols, str):
        return (covariate_cols,)
    return tuple(covariate_cols)


def _design_matrix(
    df: pd.DataFrame,
    outcome_col: str = "Y",
    covariate_cols: str | Iterable[str] | None = ("X1", "X2"),
    add_intercept: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    if outcome_col not in df.columns:
        raise ValueError(f"Outcome column '{outcome_col}' not found in data frame.")

    covariate_names = _normalize_covariate_cols(covariate_cols)
    missing_covariates = [col for col in covariate_names if col not in df.columns]
    if missing_covariates:
        raise ValueError(
            "Covariate columns not found in data frame: " + ", ".join(missing_covariates)
        )

    y = df[outcome_col].to_numpy(dtype=float)
    n = len(df)
    parts: list[np.ndarray] = []
    if add_intercept:
        parts.append(np.ones((n, 1), dtype=float))
    if covariate_names:
        parts.append(df.loc[:, list(covariate_names)].to_numpy(dtype=float))
    if not parts:
        raise ValueError("Design matrix must include at least one column.")

    X = parts[0] if len(parts) == 1 else np.column_stack(parts)
    return X, y


def _fit_binomial_irls(
    X: np.ndarray,
    y: np.ndarray,
    max_iter: int = 100,
    tol: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    beta = np.zeros(X.shape[1], dtype=float)
    eps = 1e-12

    for _ in range(max_iter):
        eta = X @ beta
        mu = np.clip(expit(eta), eps, 1.0 - eps)
        w = np.clip(mu * (1.0 - mu), eps, None)
        z = eta + (y - mu) / w

        WX = X * w[:, None]
        XtWX = X.T @ WX
        XtWz = X.T @ (w * z)

        beta_new = np.linalg.solve(XtWX, XtWz)
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break
        beta = beta_new

    fitted = np.clip(expit(X @ beta), eps, 1.0 - eps)
    weights = np.clip(fitted * (1.0 - fitted), eps, None)
    return beta, fitted, weights


def _fit_binomial_glm_rstyle(
    X: np.ndarray,
    y: np.ndarray,
    max_iter: int = 100,
    tol: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Mirror R glm() binomial fit outputs used by fit_null_glm_Binary_SPA.

    R's glm object stores working weights from the last IWLS update step. Those
    can differ slightly from fitted*(1-fitted), and STAAR's SPA prefilter
    covariance depends on these exact stored weights.
    """

    eps = 1e-12
    prior_weights = np.ones_like(y, dtype=float)
    mu = np.clip((prior_weights * y + 0.5) / (prior_weights + 1.0), eps, 1.0 - eps)
    eta = np.log(mu / (1.0 - mu))

    beta = np.zeros(X.shape[1], dtype=float)
    dev_old = np.inf
    w_last = np.ones_like(y, dtype=float)

    for _ in range(max_iter):
        var_mu = np.clip(mu * (1.0 - mu), eps, None)
        mu_eta = np.clip(mu * (1.0 - mu), eps, None)
        z = eta + (y - mu) / mu_eta
        w = np.sqrt(prior_weights * (mu_eta**2) / var_mu)

        Xw = X * w[:, None]
        zw = z * w
        beta = np.linalg.lstsq(Xw, zw, rcond=None)[0]

        eta = X @ beta
        mu = np.clip(expit(eta), eps, 1.0 - eps)
        dev = -2.0 * np.sum(y * np.log(mu) + (1.0 - y) * np.log(1.0 - mu))
        w_last = w

        if np.abs(dev - dev_old) / (0.1 + np.abs(dev)) < tol:
            break
        dev_old = dev

    weights = w_last**2
    return beta, mu, weights


def fit_null_glm(
    df: pd.DataFrame,
    *,
    outcome_col: str = "Y",
    covariate_cols: str | Iterable[str] | None = ("X1", "X2"),
    add_intercept: bool = True,
) -> NullModelGLM:
    X, y = _design_matrix(
        df,
        outcome_col=outcome_col,
        covariate_cols=covariate_cols,
        add_intercept=add_intercept,
    )
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    fitted = X @ beta
    residuals = y - fitted
    df_resid = X.shape[0] - X.shape[1]
    if df_resid <= 0:
        raise ValueError(
            "Residual degrees of freedom must be positive (n_samples must exceed model parameters)."
        )
    dispersion = float(np.sum(residuals**2) / df_resid)
    weights = np.ones_like(y)
    return NullModelGLM(
        X=X,
        y=y,
        fitted=fitted,
        residuals=residuals,
        beta=beta,
        dispersion=dispersion,
        weights=weights,
        family="gaussian",
    )


def fit_null_glm_binary_spa(
    df: pd.DataFrame,
    max_iter: int = 100,
    tol: float = 1e-8,
    *,
    outcome_col: str = "Y",
    covariate_cols: str | Iterable[str] | None = ("X1", "X2"),
    add_intercept: bool = True,
) -> NullModelGLM:
    X, y = _design_matrix(
        df,
        outcome_col=outcome_col,
        covariate_cols=covariate_cols,
        add_intercept=add_intercept,
    )
    y = y.astype(float)
    if not np.all(np.isin(y, [0.0, 1.0])):
        raise ValueError("Binary SPA null model requires Y coded as 0/1.")

    beta, fitted, weights = _fit_binomial_glm_rstyle(
        X=X,
        y=y,
        max_iter=max_iter,
        tol=tol,
    )
    residuals = y - fitted

    XtWX = X.T @ (X * weights[:, None])
    XXWX_inv = X @ np.linalg.inv(XtWX)
    XW = (X * weights[:, None]).T

    return NullModelGLM(
        X=X,
        y=y,
        fitted=fitted,
        residuals=residuals,
        beta=beta,
        dispersion=1.0,
        weights=weights,
        family="binomial",
        use_SPA=True,
        XW=XW,
        XXWX_inv=XXWX_inv,
    )


def _reml_nll(log_tau: np.ndarray, X: np.ndarray, y: np.ndarray, kins: sp.csc_matrix) -> float:
    tau1 = float(np.exp(log_tau[0]))
    tau2 = float(np.exp(log_tau[1]))
    n = X.shape[0]
    sigma = sp.identity(n, format="csc")
    sigma.data *= tau1
    if tau2 != 0.0:
        sigma = sigma + kins * tau2
    lu = splu(sigma)
    diag_u = lu.U.diagonal()
    logdet = float(np.sum(np.log(np.abs(diag_u))))

    sigma_i_x = lu.solve(X)
    xt_sigma_i_x = X.T @ sigma_i_x
    sign, logdet_x = np.linalg.slogdet(xt_sigma_i_x)
    if sign <= 0:
        return np.inf
    cov = np.linalg.inv(xt_sigma_i_x)
    sigma_i_y = lu.solve(y)
    beta = cov @ (X.T @ sigma_i_y)
    residual = y - X @ beta
    sigma_i_r = lu.solve(residual)
    quad = float(residual @ sigma_i_r)

    p = X.shape[1]
    ll = -0.5 * (logdet + logdet_x + quad + (n - p) * np.log(2 * np.pi))
    return -ll


def _reml_nll_binomial_tau(
    log_tau: float,
    X: np.ndarray,
    z: np.ndarray,
    d_inv: np.ndarray,
    kins: sp.csc_matrix,
) -> float:
    tau = float(np.exp(log_tau))
    sigma = sp.diags(d_inv, format="csc")
    if tau != 0.0:
        sigma = sigma + kins * tau
    lu = splu(sigma)

    sigma_i_x = lu.solve(X)
    xt_sigma_i_x = X.T @ sigma_i_x
    sign, logdet_x = np.linalg.slogdet(xt_sigma_i_x)
    if sign <= 0:
        return np.inf

    cov = np.linalg.inv(xt_sigma_i_x)
    sigma_i_z = lu.solve(z)
    beta = cov @ (X.T @ sigma_i_z)
    residual = z - X @ beta
    sigma_i_r = lu.solve(residual)
    quad = float(residual @ sigma_i_r)

    diag_u = lu.U.diagonal()
    logdet = float(np.sum(np.log(np.abs(diag_u))))
    return 0.5 * (logdet + logdet_x + quad)


def estimate_tau_reml(X: np.ndarray, y: np.ndarray, kins: sp.csc_matrix) -> np.ndarray:
    init = np.log([1.0, 1.0])
    result = minimize(
        _reml_nll,
        init,
        args=(X, y, kins),
        method="Nelder-Mead",
        options={"maxiter": 200, "xatol": 1e-8, "fatol": 1e-8},
    )
    return np.exp(result.x)


def fit_null_glmmkin(
    df: pd.DataFrame,
    kins: sp.csc_matrix,
    sparse_kins: bool,
    precomputed_cov: np.ndarray | None = None,
    precomputed_scaled_residuals: np.ndarray | None = None,
    precomputed_theta: np.ndarray | None = None,
    *,
    outcome_col: str = "Y",
    covariate_cols: str | Iterable[str] | None = ("X1", "X2"),
    add_intercept: bool = True,
) -> NullModelGLMMKin:
    X, y = _design_matrix(
        df,
        outcome_col=outcome_col,
        covariate_cols=covariate_cols,
        add_intercept=add_intercept,
    )
    if precomputed_theta is None:
        theta = estimate_tau_reml(X, y, kins)
    else:
        theta = np.asarray(precomputed_theta, dtype=float).reshape(-1)
        if theta.size != 2:
            raise ValueError("precomputed_theta must contain exactly two components.")
        if np.any(theta <= 0):
            raise ValueError("precomputed_theta components must be positive.")

    n = X.shape[0]
    sigma = sp.identity(n, format="csc")
    sigma.data *= theta[0]
    sigma = sigma + kins * theta[1]
    lu = splu(sigma)

    sigma_i_x = lu.solve(X)
    xt_sigma_i_x = X.T @ sigma_i_x
    cov = np.linalg.inv(xt_sigma_i_x)
    sigma_i_y = lu.solve(y)
    beta = cov @ (X.T @ sigma_i_y)

    # Match GMMAT's eta computation in R_fitglmm_ai
    diag_sigma = theta[0]
    fitted = y - diag_sigma * (sigma_i_y - sigma_i_x @ beta)
    residuals = y - fitted
    if precomputed_scaled_residuals is None:
        scaled_residuals = residuals / theta[0]
    else:
        scaled_residuals = precomputed_scaled_residuals

    return NullModelGLMMKin(
        X=X,
        y=y,
        fitted=fitted,
        residuals=residuals,
        scaled_residuals=scaled_residuals,
        beta=beta,
        theta=theta,
        cov=cov,
        sigma_solver=lu,
        Sigma_iX=sigma_i_x,
        sparse_kins=sparse_kins,
        precomputed_cov=precomputed_cov,
    )


def fit_null_glmmkin_binary_spa(
    df: pd.DataFrame,
    kins: sp.csc_matrix,
    sparse_kins: bool,
    max_iter: int = 80,
    tol: float = 1e-8,
    *,
    outcome_col: str = "Y",
    covariate_cols: str | Iterable[str] | None = ("X1", "X2"),
    add_intercept: bool = True,
) -> NullModelGLMMKinBinarySPA:
    X, y = _design_matrix(
        df,
        outcome_col=outcome_col,
        covariate_cols=covariate_cols,
        add_intercept=add_intercept,
    )
    y = y.astype(float)
    if not np.all(np.isin(y, [0.0, 1.0])):
        raise ValueError("Binary SPA null model requires Y coded as 0/1.")

    eps = 1e-9
    # R-style binomial initialization matches the null-model start used by STAAR's
    # binary workflows more closely than the plain IRLS start.
    beta, _, _ = _fit_binomial_glm_rstyle(
        X=X,
        y=y,
        max_iter=100,
        tol=1e-8,
    )
    eta = X @ beta
    tau = 1.0

    for _ in range(max_iter):
        fitted = np.clip(expit(eta), eps, 1.0 - eps)
        weights = np.clip(fitted * (1.0 - fitted), eps, None)
        z = eta + (y - fitted) / weights
        d_inv = 1.0 / weights

        try:
            opt = minimize_scalar(
                _reml_nll_binomial_tau,
                bounds=(-10.0, 5.0),
                method="bounded",
                args=(X, z, d_inv, kins),
                options={"xatol": 1e-6, "maxiter": 200},
            )
            if opt.success:
                tau_new = float(np.exp(opt.x))
            else:
                tau_new = tau
        except Exception:
            tau_new = tau

        sigma = sp.diags(d_inv, format="csc")
        if tau_new != 0.0:
            sigma = sigma + kins * tau_new
        lu = splu(sigma)

        sigma_i_x = lu.solve(X)
        xt_sigma_i_x = X.T @ sigma_i_x
        cov = np.linalg.inv(xt_sigma_i_x)
        sigma_i_z = lu.solve(z)
        beta_new = cov @ (X.T @ sigma_i_z)
        residual_working = z - X @ beta_new
        u = tau_new * (kins @ lu.solve(residual_working))
        eta_new = X @ beta_new + u

        delta = float(np.max(np.abs(eta_new - eta)))
        eta = eta_new
        beta = beta_new
        tau = tau_new
        if delta < tol:
            break

    fitted = np.clip(expit(eta), eps, 1.0 - eps)
    residuals = y - fitted
    weights = np.clip(fitted * (1.0 - fitted), eps, None)
    d_inv = 1.0 / weights
    sigma = sp.diags(d_inv, format="csc")
    if tau != 0.0:
        sigma = sigma + kins * tau
    lu = splu(sigma)
    sigma_i_x = lu.solve(X)
    xt_sigma_i_x = X.T @ sigma_i_x
    cov = np.linalg.inv(xt_sigma_i_x)

    # Keep residuals on the Y-mu scale for SPA score construction.
    scaled_residuals = residuals

    if sparse_kins:
        XW = sigma_i_x.T
        XXWX_inv = X @ cov
    else:
        XtWX = X.T @ (X * weights[:, None])
        XXWX_inv = X @ np.linalg.inv(XtWX)
        XW = (X * weights[:, None]).T

    theta = np.array([1.0, tau], dtype=float)

    return NullModelGLMMKinBinarySPA(
        X=X,
        y=y,
        fitted=fitted,
        residuals=residuals,
        scaled_residuals=scaled_residuals,
        beta=beta,
        theta=theta,
        cov=cov,
        sigma_solver=lu,
        Sigma_iX=sigma_i_x,
        sparse_kins=sparse_kins,
        weights=weights,
        XW=XW,
        XXWX_inv=XXWX_inv,
    )
