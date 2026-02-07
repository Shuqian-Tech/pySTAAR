"""Core STAAR computations translated from the R/C++ implementation."""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np
from scipy import special, stats

from .staar_stats import cct, cct_pval, saddle


def matrix_flip(G: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Impute missing genotypes, flip to minor allele, and compute AF/MAF.

    Mirrors matrix_flip.cpp.
    """
    G = G.copy().astype(float)
    observed = G > -1
    num_observed = np.sum(observed, axis=0)

    with np.errstate(invalid="ignore", divide="ignore"):
        AF = np.sum(np.where(observed, G, 0.0), axis=0) / (2.0 * num_observed)

    low_af = AF <= 0.5
    missing_idx = np.where(~observed)
    if missing_idx[0].size > 0:
        impute_values = np.where(low_af, 0.0, 2.0)
        G[missing_idx] = impute_values[missing_idx[1]]

    flip_mask = ~low_af
    if np.any(flip_mask):
        G[:, flip_mask] = 2.0 - G[:, flip_mask]

    MAF = np.where(low_af, AF, 1.0 - AF)

    return G, AF, MAF


def _compute_weights(maf: np.ndarray, annotation_phred: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    annotation_rank = 1.0 - 10.0 ** (-annotation_phred / 10.0)

    w_1 = stats.beta.pdf(maf, 1, 25)
    w_2 = stats.beta.pdf(maf, 1, 1)
    base = stats.beta.pdf(maf, 0.5, 0.5)

    if annotation_phred.shape[1] == 0:
        w_B = np.column_stack([w_1, w_2])
        w_S = np.column_stack([w_1, w_2])
        w_A = np.column_stack([
            w_1**2 / base**2,
            w_2**2 / base**2,
        ])
    else:
        w_B_1 = np.column_stack([w_1, annotation_rank * w_1[:, None]])
        w_B_2 = np.column_stack([w_2, annotation_rank * w_2[:, None]])
        w_B = np.column_stack([w_B_1, w_B_2])

        w_S_1 = np.column_stack([w_1, np.sqrt(annotation_rank) * w_1[:, None]])
        w_S_2 = np.column_stack([w_2, np.sqrt(annotation_rank) * w_2[:, None]])
        w_S = np.column_stack([w_S_1, w_S_2])

        w_A_1 = np.column_stack([
            w_1**2 / base**2,
            annotation_rank * w_1[:, None] ** 2 / base[:, None] ** 2,
        ])
        w_A_2 = np.column_stack([
            w_2**2 / base**2,
            annotation_rank * w_2[:, None] ** 2 / base[:, None] ** 2,
        ])
        w_A = np.column_stack([w_A_1, w_A_2])

    return w_B, w_S, w_A


def _pchisq_upper(x: float, df: float) -> float:
    x = float(x)
    df = float(df)
    if x <= 0.0:
        return 1.0
    if df == 1.0:
        return float(special.erfc(math.sqrt(0.5 * x)))
    return float(special.gammaincc(0.5 * df, 0.5 * x))


def _pt_upper(x: float, df: float) -> float:
    return float(special.stdtr(float(df), -float(x)))


def _pnorm_upper(x: float) -> float:
    return float(special.ndtr(-float(x)))


def _safe_chisq_pvalue(numerator: float, denominator: float, eps: float = 1e-12) -> float:
    if abs(denominator) <= eps:
        if abs(numerator) <= eps:
            return 1.0
        return 0.0
    stat = (numerator**2) / denominator
    return _pchisq_upper(stat, 1)


def _safe_ratio_square(numerator: float, denominator: float, eps: float = 1e-12) -> float:
    if abs(denominator) <= eps:
        if abs(numerator) <= eps:
            return 0.0
        return np.inf
    return (numerator**2) / denominator


def _k_binary_spa(x: float, muhat: np.ndarray, g: np.ndarray) -> float:
    xg = np.clip(x * g, -700.0, 700.0)
    return float(np.sum(-x * muhat * g + np.log1p(-muhat + muhat * np.exp(xg))))


def _k_binary_spa_alt(x: float, muhat: np.ndarray, g: np.ndarray) -> float:
    xg = np.clip(-x * g, -700.0, 700.0)
    return float(np.sum(np.log((1.0 - muhat) * np.exp(xg) + muhat)))


def _k1_binary_spa(x: float, muhat: np.ndarray, g: np.ndarray, q: float) -> float:
    xg = np.clip(-x * g, -700.0, 700.0)
    term = muhat * g / (muhat + (1.0 - muhat) * np.exp(xg))
    return float(np.sum(-muhat * g + term) - q)


def _k1_binary_spa_alt(x: float, muhat: np.ndarray, g: np.ndarray, q: float) -> float:
    ex = np.exp(np.clip(x * g, -700.0, 700.0))
    term = muhat * g * ex / (muhat * ex + (1.0 - muhat))
    return float(np.sum(-muhat * g + term) - q)


def _k2_binary_spa(x: float, muhat: np.ndarray, g: np.ndarray) -> float:
    ex = np.exp(np.clip(-x * g, -700.0, 700.0))
    num = muhat * (1.0 - muhat) * (g**2) * ex
    den = (muhat + (1.0 - muhat) * ex) ** 2
    return float(np.sum(num / den))


def _k2_binary_spa_alt(x: float, muhat: np.ndarray, g: np.ndarray) -> float:
    ex = np.exp(np.clip(x * g, -700.0, 700.0))
    num = muhat * (1.0 - muhat) * (g**2)
    den = (muhat * ex + (1.0 - muhat)) ** 2
    return float(np.sum(num / den))


def _is_bad_number(x: float) -> bool:
    return not np.isfinite(x)


def _have_same_sign(a: float, b: float) -> bool:
    return (a >= 0.0) == (b >= 0.0)


def _nr_binary_spa(
    muhat: np.ndarray,
    g: np.ndarray,
    q: float,
    init: float,
    tol: float,
    max_iter: int,
) -> float:
    xi = float(init)
    xi_update = float(init)

    k1 = _k1_binary_spa(xi, muhat, g, q)
    if abs(k1) > tol:
        k2 = _k2_binary_spa(xi, muhat, g)
        xi_update = xi - k1 / k2

    no_iter = 0
    while (
        np.isfinite(xi_update)
        and abs(xi_update - xi) > tol
        and abs(_k1_binary_spa(xi_update, muhat, g, q)) > tol
        and no_iter < max_iter
    ):
        no_iter += 1
        xi = xi_update

        numerator = _k1_binary_spa(xi, muhat, g, q)
        if _is_bad_number(numerator):
            numerator = _k1_binary_spa_alt(xi, muhat, g, q)

        denominator = _k2_binary_spa(xi, muhat, g)
        if _is_bad_number(denominator):
            denominator = _k2_binary_spa_alt(xi, muhat, g)

        xi_update = xi - numerator / denominator

    if _is_bad_number(xi_update):
        xi_update = xi

    if no_iter == max_iter:
        xhat = float(xi_update)
        w = math.sqrt(max(0.0, 2.0 * (xhat * q - _k_binary_spa(xhat, muhat, g))))
        if xhat < 0.0:
            w = -w
        ki = xhat * math.sqrt(max(0.0, _k2_binary_spa(xhat, muhat, g)))
        if _is_bad_number(ki):
            ki = xhat * math.sqrt(max(0.0, _k2_binary_spa_alt(xhat, muhat, g)))
        ratio = ki / w if w != 0.0 else np.nan
        z = w + math.log(ratio) / w if w != 0.0 and ratio > 0.0 else np.inf
        if abs(z) < 38.0:
            xi_update = 0.0

    return float(xi_update)


def _golden_section_search_sign_change(
    a: float,
    b: float,
    muhat: np.ndarray,
    g: np.ndarray,
    q: float,
    tol: float,
    max_iter: int,
) -> Tuple[float, float]:
    x0 = float(a)
    x1 = float(b)
    iter_count = 0
    phi = (1.0 + math.sqrt(5.0)) / 2.0

    k1x0 = _k1_binary_spa(x0, muhat, g, q)
    if _is_bad_number(k1x0):
        k1x0 = _k1_binary_spa_alt(x0, muhat, g, q)
    k1x1 = _k1_binary_spa(x1, muhat, g, q)
    if _is_bad_number(k1x1):
        k1x1 = _k1_binary_spa_alt(x1, muhat, g, q)

    while iter_count < max_iter and abs(x1 - x0) > tol and _have_same_sign(k1x0, k1x1):
        x1 = b - (b - a) / phi
        x0 = a + (b - a) / phi

        k1x0 = _k1_binary_spa(x0, muhat, g, q)
        if _is_bad_number(k1x0):
            k1x0 = _k1_binary_spa_alt(x0, muhat, g, q)

        k1x1 = _k1_binary_spa(x1, muhat, g, q)
        if _is_bad_number(k1x1):
            k1x1 = _k1_binary_spa_alt(x1, muhat, g, q)

        if k1x1 < k1x0:
            b = x0
        else:
            a = x1
        iter_count += 1

    return float(x0), float(x1)


def _bisection_binary_spa(
    muhat: np.ndarray,
    g: np.ndarray,
    q: float,
    xmin: float,
    xmax: float,
    tol: float,
) -> float:
    xupper = float(xmax)
    xlower = float(xmin)
    x0 = 0.0
    k1x0 = 1.0

    k1left = _k1_binary_spa(xlower, muhat, g, q)
    if _is_bad_number(k1left):
        k1left = _k1_binary_spa_alt(xlower, muhat, g, q)

    k1right = _k1_binary_spa(xupper, muhat, g, q)
    if _is_bad_number(k1right):
        k1right = _k1_binary_spa_alt(xupper, muhat, g, q)

    if (not _is_bad_number(k1left)) and (not _is_bad_number(k1right)) and (not _have_same_sign(k1left, k1right)):
        while abs(xupper - xlower) > tol and abs(k1x0) > tol:
            x0 = (xupper + xlower) / 2.0
            k1x0 = _k1_binary_spa(x0, muhat, g, q)
            if _is_bad_number(k1x0):
                k1x0 = _k1_binary_spa_alt(x0, muhat, g, q)

            if k1x0 == 0.0:
                break
            if _have_same_sign(k1left, k1x0):
                xlower = x0
            else:
                xupper = x0

            k1left = _k1_binary_spa(xlower, muhat, g, q)
            if _is_bad_number(k1left):
                k1left = _k1_binary_spa_alt(xlower, muhat, g, q)

            k1right = _k1_binary_spa(xupper, muhat, g, q)
            if _is_bad_number(k1right):
                k1right = _k1_binary_spa_alt(xupper, muhat, g, q)

    return float(x0)


def _saddle_binary_spa(
    q: float,
    muhat: np.ndarray,
    g: np.ndarray,
    tol: float,
    max_iter: int,
    lower: bool,
) -> float:
    xhat = _nr_binary_spa(muhat, g, q, init=0.0, tol=tol, max_iter=max_iter)
    w = math.sqrt(max(0.0, 2.0 * (xhat * q - _k_binary_spa(xhat, muhat, g))))
    if _is_bad_number(w):
        w = math.sqrt(max(0.0, 2.0 * (xhat * q - _k_binary_spa_alt(xhat, muhat, g))))

    if xhat < 0.0:
        w = -w

    ki = xhat * math.sqrt(max(0.0, _k2_binary_spa(xhat, muhat, g)))
    if _is_bad_number(ki):
        ki = xhat * math.sqrt(max(0.0, _k2_binary_spa_alt(xhat, muhat, g)))

    ratio = ki / w if w != 0.0 else np.nan
    if abs(xhat) < 1e-4 or w == 0.0 or ratio <= 0.0 or (not np.isfinite(ratio)):
        return 1.0

    z = w + math.log(ratio) / w
    return float(special.ndtr(z) if lower else special.ndtr(-z))


def _saddle_binary_spa_bisection(
    q: float,
    muhat: np.ndarray,
    g: np.ndarray,
    tol: float,
    max_iter: int,
    xmin: float,
    xmax: float,
    lower: bool,
) -> float:
    x0, x1 = _golden_section_search_sign_change(xmin, xmax, muhat, g, q, tol, max_iter)
    xmin_use, xmax_use = (x0, x1) if x0 < x1 else (x1, x0)
    xhat = _bisection_binary_spa(muhat, g, q, xmin_use, xmax_use, tol)

    w = math.sqrt(max(0.0, 2.0 * (xhat * q - _k_binary_spa(xhat, muhat, g))))
    if _is_bad_number(w):
        w = math.sqrt(max(0.0, 2.0 * (xhat * q - _k_binary_spa_alt(xhat, muhat, g))))
    if xhat < 0.0:
        w = -w

    ki = xhat * math.sqrt(max(0.0, _k2_binary_spa(xhat, muhat, g)))
    if _is_bad_number(ki):
        ki = xhat * math.sqrt(max(0.0, _k2_binary_spa_alt(xhat, muhat, g)))

    ratio = ki / w if w != 0.0 else np.nan
    if abs(xhat) < 1e-4 or w == 0.0 or ratio <= 0.0 or (not np.isfinite(ratio)):
        return 1.0

    z = w + math.log(ratio) / w
    return float(special.ndtr(z) if lower else special.ndtr(-z))


def _staartest_burden_binary_spa(
    G: np.ndarray,
    XW: np.ndarray,
    XXWX_inv: np.ndarray,
    residuals: np.ndarray,
    muhat: np.ndarray,
    weights_B: np.ndarray,
    tol: float,
    max_iter: int,
) -> np.ndarray:
    wn = weights_B.shape[1]
    res = np.ones(wn, dtype=float)

    G_tilde = G - XXWX_inv @ (XW @ G)
    x = residuals @ G
    G_cumu = G_tilde @ weights_B

    xmin = -100.0
    xmax = 100.0

    for i in range(wn):
        sum0 = float(np.sum(x * weights_B[:, i]))
        g_col = G_cumu[:, i]
        res[i] = _burden_spa_two_sided_pvalue(
            score=sum0,
            muhat=muhat,
            g_col=g_col,
            tol=tol,
            max_iter=max_iter,
            xmin=xmin,
            xmax=xmax,
        )

    return res


def _burden_spa_two_sided_pvalue(
    score: float,
    muhat: np.ndarray,
    g_col: np.ndarray,
    tol: float,
    max_iter: int,
    xmin: float,
    xmax: float,
) -> float:
    q_abs = abs(float(score))

    respart1 = _saddle_binary_spa(q_abs, muhat, g_col, tol, max_iter, lower=False)
    if np.isnan(respart1) or respart1 == 1.0:
        respart1 = _saddle_binary_spa_bisection(
            q_abs, muhat, g_col, tol, max_iter, xmin, xmax, lower=False
        )
    if np.isnan(respart1) or respart1 == 1.0:
        return 1.0

    respart2 = _saddle_binary_spa(-q_abs, muhat, g_col, tol, max_iter, lower=True)
    if np.isnan(respart2) or respart2 == 1.0:
        respart2 = _saddle_binary_spa_bisection(
            -q_abs, muhat, g_col, tol, max_iter, xmin, xmax, lower=True
        )
    if np.isnan(respart2) or respart2 == 1.0:
        return 1.0

    return float(min(1.0, respart1 + respart2))


def _binary_spa_filter_covariance(G: np.ndarray, obj_nullmodel) -> np.ndarray | None:
    precomputed_cov_filter = getattr(obj_nullmodel, "precomputed_cov_filter", None)
    if precomputed_cov_filter is not None:
        Cov = np.asarray(precomputed_cov_filter, dtype=float)
        if Cov.shape == (G.shape[1], G.shape[1]):
            return 0.5 * (Cov + Cov.T)

    if obj_nullmodel.relatedness:
        P = getattr(obj_nullmodel, "P", None)
        if P is not None:
            PG = P @ G
            Cov = PG.T @ G
            return 0.5 * (Cov + Cov.T)

        Sigma_i = getattr(obj_nullmodel, "Sigma_i", None)
        Sigma_iX = getattr(obj_nullmodel, "Sigma_iX", None)
        cov = getattr(obj_nullmodel, "cov", None)
        if Sigma_i is not None and Sigma_iX is not None and cov is not None:
            Sigma_iG = Sigma_i @ G
            tSigma_iX_G = Sigma_iX.T @ G
            Cov = Sigma_iG.T @ G - tSigma_iX_G.T @ cov @ tSigma_iX_G
            return 0.5 * (Cov + Cov.T)

        sigma_solver = getattr(obj_nullmodel, "sigma_solver", None)
        Sigma_iX = getattr(obj_nullmodel, "Sigma_iX", None)
        cov = getattr(obj_nullmodel, "cov", None)
        if sigma_solver is not None and Sigma_iX is not None and cov is not None:
            Sigma_iG = sigma_solver.solve(G)
            tSigma_iX_G = Sigma_iX.T @ G
            Cov = Sigma_iG.T @ G - tSigma_iX_G.T @ cov @ tSigma_iX_G
            return 0.5 * (Cov + Cov.T)

        return None

    X = getattr(obj_nullmodel, "X", None)
    working = getattr(obj_nullmodel, "weights", None)
    family = getattr(obj_nullmodel, "family", None)
    if X is None or working is None or family is None:
        return None

    if family == "gaussian":
        tX_G = X.T @ G
        Cov = G.T @ G - tX_G.T @ np.linalg.inv(X.T @ X) @ tX_G
    else:
        WX = X * working[:, None]
        tX_G = X.T @ (working[:, None] * G)
        Cov = (working[:, None] * G).T @ G - tX_G.T @ np.linalg.inv(X.T @ WX) @ tX_G
    return 0.5 * (Cov + Cov.T)


def _staartest_burden_binary_spa_filtered(
    G: np.ndarray,
    XW: np.ndarray,
    XXWX_inv: np.ndarray,
    residuals: np.ndarray,
    muhat: np.ndarray,
    weights_B: np.ndarray,
    Cov: np.ndarray,
    p_filter_cutoff: float,
    tol: float,
    max_iter: int,
) -> np.ndarray:
    wn = weights_B.shape[1]
    res = np.ones(wn, dtype=float)

    G_tilde = G - XXWX_inv @ (XW @ G)
    x = residuals @ G
    G_cumu = G_tilde @ weights_B

    xmin = -100.0
    xmax = 100.0

    for i in range(wn):
        w = weights_B[:, i]
        score = float(np.sum(x * w))
        variance = float(w.T @ Cov @ w)
        p_normal = _safe_chisq_pvalue(score, variance)
        res[i] = p_normal

        if p_normal < p_filter_cutoff:
            res[i] = _burden_spa_two_sided_pvalue(
                score=score,
                muhat=muhat,
                g_col=G_cumu[:, i],
                tol=tol,
                max_iter=max_iter,
                xmin=xmin,
                xmax=xmax,
            )

    return res


def _cct_ignore_one_na(pvalues: np.ndarray) -> float:
    pvalues = np.asarray(pvalues, dtype=float)
    pvalues = pvalues[~np.isnan(pvalues)]
    if pvalues.size == 0:
        return 1.0
    pvalues = pvalues[pvalues < 1.0]
    if pvalues.size == 0:
        return 1.0
    return float(cct(pvalues))


def _staartest_pvalues(
    G: np.ndarray,
    residuals: np.ndarray,
    Cov: np.ndarray,
    weights_B: np.ndarray,
    weights_S: np.ndarray,
    weights_A: np.ndarray,
    mac: np.ndarray,
    sigma: float | None,
    fam: int,
    df_resid: int | None = None,
) -> np.ndarray:
    n = G.shape[0]
    un = G.shape[1]
    wn = weights_B.shape[1]
    res = np.zeros(3 * wn, dtype=float)

    x = residuals @ G

    id_veryrare = np.where(mac <= 10)[0]
    id_common = np.where(mac > 10)[0]

    n0 = len(id_veryrare)
    n1 = len(id_common)

    pseq = np.zeros(un, dtype=float)
    wseq = np.zeros(un, dtype=float)

    if n1 > 0:
        for idx_pos, idx in enumerate(id_common):
            if fam == 0:
                SSR = x[idx] ** 2 / Cov[idx, idx]
                if df_resid is None:
                    raise ValueError("df_resid is required for gaussian ACAT calculations.")
                SST = sigma**2 * df_resid - SSR
                tval = math.sqrt(SSR / SST * (df_resid - 1))
                pseq[idx_pos] = 2.0 * _pt_upper(tval, df_resid - 1)
            else:
                pseq[idx_pos] = _safe_chisq_pvalue(x[idx], Cov[idx, idx])

    for i in range(wn):
        # SKAT
        w = weights_S[:, i]
        Covw = Cov * np.outer(w, w)
        if sigma is not None:
            Covw = Covw * (sigma**2)

        sum0 = np.sum((x**2) * (w**2))
        eigenvals = np.linalg.eigvalsh(Covw)
        eigenvals[eigenvals < 1e-8] = 0.0
        p_skat = saddle(sum0, eigenvals)
        if p_skat == 2:
            c1 = np.trace(Covw)
            Covw2 = Covw @ Covw
            c2 = np.trace(Covw2)
            Covw4 = Covw2 @ Covw2
            c4 = np.trace(Covw4)
            sum0_adj = (sum0 - c1) / math.sqrt(2 * c2)
            lval = (c2**2) / c4
            p_skat = _pchisq_upper(sum0_adj * math.sqrt(2 * lval) + lval, lval)
        res[i] = p_skat

        # Burden
        w = weights_B[:, i]
        Covw = Cov * np.outer(w, w)
        if sigma is not None:
            Covw = Covw * (sigma**2)
        sumw = np.sum(Covw)
        sum0 = np.sum(x * w)
        sumx = _safe_ratio_square(sum0, sumw)
        res[wn + i] = _pchisq_upper(sumx, 1)

        # ACAT-V
        if n1 > 0:
            for idx_pos, idx in enumerate(id_common):
                wseq[idx_pos] = weights_A[idx, i]
        if n0 == 0:
            res[2 * wn + i] = cct_pval(pseq[:n1], wseq[:n1])
        else:
            sum0 = np.sum(x[id_veryrare] * weights_B[id_veryrare, i])
            sumw = np.sum(weights_A[id_veryrare, i])
            sumx = _safe_ratio_square(sum0, np.sum(Covw[np.ix_(id_veryrare, id_veryrare)]))
            pseq[n1] = _pchisq_upper(sumx, 1)
            wseq[n1] = sumw / n0
            res[2 * wn + i] = cct_pval(pseq[: n1 + 1], wseq[: n1 + 1])

    return res


def staar_binary_spa(
    genotype: np.ndarray,
    obj_nullmodel,
    annotation_phred: np.ndarray | None = None,
    rare_maf_cutoff: float = 0.01,
    rv_num_cutoff: int = 2,
    rv_num_cutoff_max: float = 1e9,
    tol: float = float(np.finfo(float).eps ** 0.25),
    max_iter: int = 1000,
    SPA_p_filter: bool = False,
    p_filter_cutoff: float = 0.05,
) -> Dict[str, object]:
    genotype = _as_2d_matrix(genotype, "genotype")
    if genotype.shape[1] == 1:
        raise ValueError("Number of rare variant in the set is less than 2!")

    annotation_phred = _as_annotation_matrix(annotation_phred, genotype.shape[1])
    if annotation_phred.size > 0 and genotype.shape[1] != annotation_phred.shape[0]:
        raise ValueError("Dimensions don't match for genotype and annotation!")

    G_flip, _, maf = matrix_flip(genotype)
    rv_label = (maf < rare_maf_cutoff) & (maf > 0)
    G = G_flip[:, rv_label]
    annotation_phred = annotation_phred[rv_label, :]

    if np.sum(rv_label) >= rv_num_cutoff_max:
        raise ValueError("Number of rare variant in the set is more than rv_num_cutoff_max!")
    if np.sum(rv_label) < rv_num_cutoff:
        raise ValueError("Number of rare variant in the set is less than rv_num_cutoff!")

    maf = maf[rv_label]
    w_B, _, _ = _compute_weights(maf, annotation_phred)

    if obj_nullmodel.XW is None or obj_nullmodel.XXWX_inv is None:
        raise ValueError("Binary SPA null model requires XW and XXWX_inv.")
    if obj_nullmodel.relatedness:
        if getattr(obj_nullmodel, "scaled_residuals", None) is None:
            raise ValueError("Related binary SPA null model requires scaled_residuals.")
        residuals = obj_nullmodel.scaled_residuals
    else:
        residuals = obj_nullmodel.y - obj_nullmodel.fitted
    muhat = obj_nullmodel.fitted
    if SPA_p_filter:
        if not (0.0 < p_filter_cutoff <= 1.0):
            raise ValueError("p_filter_cutoff must be in (0, 1].")
        Cov_filter = _binary_spa_filter_covariance(G, obj_nullmodel)
        if Cov_filter is not None:
            pvalues = _staartest_burden_binary_spa_filtered(
                G=G,
                XW=obj_nullmodel.XW,
                XXWX_inv=obj_nullmodel.XXWX_inv,
                residuals=residuals,
                muhat=muhat,
                weights_B=w_B,
                Cov=Cov_filter,
                p_filter_cutoff=p_filter_cutoff,
                tol=tol,
                max_iter=max_iter,
            )
        else:
            pvalues = _staartest_burden_binary_spa(
                G=G,
                XW=obj_nullmodel.XW,
                XXWX_inv=obj_nullmodel.XXWX_inv,
                residuals=residuals,
                muhat=muhat,
                weights_B=w_B,
                tol=tol,
                max_iter=max_iter,
            )
    else:
        pvalues = _staartest_burden_binary_spa(
            G=G,
            XW=obj_nullmodel.XW,
            XXWX_inv=obj_nullmodel.XXWX_inv,
            residuals=residuals,
            muhat=muhat,
            weights_B=w_B,
            tol=tol,
            max_iter=max_iter,
        )

    num_variant = int(np.sum(rv_label))
    num_annotation = annotation_phred.shape[1] + 1

    results = {}
    results["num_variant"] = float(num_variant)
    results["cMAC"] = float(np.sum(G))
    results["RV_label"] = rv_label
    results["results_STAAR_B"] = _cct_ignore_one_na(pvalues)

    def _assemble(start_idx: int, prefix: str, label: str) -> Dict[str, float]:
        segment = pvalues[start_idx : start_idx + num_annotation]
        agg = _cct_ignore_one_na(segment)
        values = list(segment) + [agg]
        cols = [
            f"{prefix}",
            *[f"{prefix}-{c}" for c in _annotation_columns(annotation_phred)],
            label,
        ]
        return dict(zip(cols, map(float, values)))

    results["results_STAAR_B_1_25"] = _assemble(0, "Burden(1,25)", "STAAR-B(1,25)")
    results["results_STAAR_B_1_1"] = _assemble(num_annotation, "Burden(1,1)", "STAAR-B(1,1)")

    return results


def _staartest_pvalues_smmat(
    G: np.ndarray,
    residuals: np.ndarray,
    Cov: np.ndarray,
    weights_B: np.ndarray,
    weights_S: np.ndarray,
    weights_A: np.ndarray,
    mac: np.ndarray,
) -> np.ndarray:
    """SMMAT-style variant-set p-values used by STAAR_cond."""
    un = G.shape[1]
    wn = weights_B.shape[1]
    res = np.zeros(3 * wn, dtype=float)

    x = residuals @ G

    id_veryrare = np.where(mac <= 10)[0]
    id_common = np.where(mac > 10)[0]

    n0 = len(id_veryrare)
    n1 = len(id_common)

    pseq = np.zeros(un, dtype=float)
    wseq = np.zeros(un, dtype=float)

    if n1 > 0:
        for idx_pos, idx in enumerate(id_common):
            pseq[idx_pos] = _safe_chisq_pvalue(x[idx], Cov[idx, idx])

    for i in range(wn):
        # SKAT
        w = weights_S[:, i]
        Covw = Cov * np.outer(w, w)

        sum0 = np.sum((x**2) * (w**2))
        eigenvals = np.linalg.eigvalsh(Covw)
        eigenvals[eigenvals < 1e-8] = 0.0
        p_skat = saddle(sum0, eigenvals)
        if p_skat == 2:
            c1 = np.trace(Covw)
            Covw2 = Covw @ Covw
            c2 = np.trace(Covw2)
            Covw4 = Covw2 @ Covw2
            c4 = np.trace(Covw4)
            sum0_adj = (sum0 - c1) / math.sqrt(2 * c2)
            lval = (c2**2) / c4
            p_skat = _pchisq_upper(sum0_adj * math.sqrt(2 * lval) + lval, lval)
        res[i] = p_skat

        # Burden
        w = weights_B[:, i]
        Covw = Cov * np.outer(w, w)
        sumw = np.sum(Covw)
        sum0 = np.sum(x * w)
        sumx = _safe_ratio_square(sum0, sumw)
        res[wn + i] = _pchisq_upper(sumx, 1)

        # ACAT-V
        if n1 > 0:
            for idx_pos, idx in enumerate(id_common):
                wseq[idx_pos] = weights_A[idx, i]
        if n0 == 0:
            res[2 * wn + i] = cct_pval(pseq[:n1], wseq[:n1])
        else:
            sum0 = np.sum(x[id_veryrare] * weights_B[id_veryrare, i])
            sumw = np.sum(weights_A[id_veryrare, i])
            sumx = _safe_ratio_square(sum0, np.sum(Covw[np.ix_(id_veryrare, id_veryrare)]))
            pseq[n1] = _pchisq_upper(sumx, 1)
            wseq[n1] = sumw / n0
            res[2 * wn + i] = cct_pval(pseq[: n1 + 1], wseq[: n1 + 1])

    return res


def _as_annotation_matrix(annotation_phred, num_variants: int) -> np.ndarray:
    if annotation_phred is None:
        return np.empty((num_variants, 0), dtype=float)
    arr = np.asarray(annotation_phred, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return arr


def _as_2d_matrix(arr: np.ndarray | list | tuple, name: str) -> np.ndarray:
    out = np.asarray(arr, dtype=float)
    if out.ndim == 1:
        out = out.reshape(-1, 1)
    if out.ndim != 2:
        raise ValueError(f"{name} is not a matrix!")
    return out


def _residualize(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def _build_conditional_design(
    genotype_adj: np.ndarray, X_base: np.ndarray, method_cond: str
) -> np.ndarray:
    if method_cond == "optimal":
        return np.column_stack([genotype_adj, X_base])
    if method_cond == "naive":
        return np.column_stack([np.ones(genotype_adj.shape[0]), genotype_adj])
    raise ValueError("method_cond must be one of {'optimal', 'naive'}.")


def _compute_related_covariance(G: np.ndarray, obj_nullmodel) -> np.ndarray:
    P = getattr(obj_nullmodel, "P", None)
    if P is not None:
        Cov = (P @ G).T @ G
        return 0.5 * (Cov + Cov.T)

    sigma_solver = getattr(obj_nullmodel, "sigma_solver", None)
    Sigma_iX = getattr(obj_nullmodel, "Sigma_iX", None)
    cov = getattr(obj_nullmodel, "cov", None)
    if sigma_solver is not None and Sigma_iX is not None and cov is not None:
        Sigma_iG = sigma_solver.solve(G)
        tSigma_iX_G = Sigma_iX.T @ G
        Cov = Sigma_iG.T @ G - tSigma_iX_G.T @ cov @ tSigma_iX_G
        return 0.5 * (Cov + Cov.T)

    Sigma_i = getattr(obj_nullmodel, "Sigma_i", None)
    if Sigma_i is not None and Sigma_iX is not None and cov is not None:
        Sigma_iG = Sigma_i @ G
        tSigma_iX_G = Sigma_iX.T @ G
        Cov = Sigma_iG.T @ G - tSigma_iX_G.T @ cov @ tSigma_iX_G
        return 0.5 * (Cov + Cov.T)

    raise ValueError("Related null model is missing covariance components for score tests.")


def _score_test_from_covariance(score: np.ndarray, Cov: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    var = np.diag(Cov).astype(float)
    var_safe = np.where(var > 0.0, var, np.nan)

    se = np.sqrt(np.where(var > 0.0, var, 0.0))
    stat = np.divide(score**2, var_safe, out=np.zeros_like(score), where=~np.isnan(var_safe))
    pvalue = np.where(np.isnan(var_safe), 1.0, special.erfc(np.sqrt(np.maximum(stat, 0.0) * 0.5)))
    return se, pvalue


def indiv_score_test_region(
    genotype: np.ndarray,
    obj_nullmodel,
    rare_maf_cutoff: float = 0.01,
    rv_num_cutoff: int = 2,
    rv_num_cutoff_max: float = 1e9,
) -> Dict[str, object]:
    genotype = _as_2d_matrix(genotype, "genotype")
    if genotype.shape[1] == 1:
        raise ValueError("Number of rare variant in the set is less than 2!")

    G_flip, _, maf = matrix_flip(genotype)
    rv_label = (maf < rare_maf_cutoff) & (maf > 0)
    G = G_flip[:, rv_label]

    if np.sum(rv_label) >= rv_num_cutoff_max:
        raise ValueError("Number of rare variant in the set is more than rv_num_cutoff_max!")
    if np.sum(rv_label) < rv_num_cutoff:
        raise ValueError("Number of rare variant in the set is less than rv_num_cutoff!")

    if obj_nullmodel.relatedness:
        residuals = obj_nullmodel.scaled_residuals
        if getattr(obj_nullmodel, "precomputed_cov", None) is not None:
            Cov = np.asarray(obj_nullmodel.precomputed_cov, dtype=float)
            if Cov.shape != (G.shape[1], G.shape[1]):
                raise ValueError("precomputed_cov has incompatible dimensions for score tests.")
            Cov = 0.5 * (Cov + Cov.T)
        else:
            Cov = _compute_related_covariance(G, obj_nullmodel)
    else:
        residuals = obj_nullmodel.y - obj_nullmodel.fitted
        X = obj_nullmodel.X
        working = obj_nullmodel.weights
        if obj_nullmodel.family == "gaussian":
            tX_G = X.T @ G
            Cov = G.T @ G - tX_G.T @ np.linalg.inv(X.T @ X) @ tX_G
            Cov = Cov * float(obj_nullmodel.dispersion)
        elif obj_nullmodel.family == "binomial":
            WX = X * working[:, None]
            tX_G = X.T @ (working[:, None] * G)
            Cov = (working[:, None] * G).T @ G - tX_G.T @ np.linalg.inv(X.T @ WX) @ tX_G
        else:
            raise NotImplementedError(
                f"Unsupported family for Indiv_Score_Test_Region: {obj_nullmodel.family}"
            )
        Cov = 0.5 * (Cov + Cov.T)

    score = np.asarray(residuals @ G, dtype=float).reshape(-1)
    se, pvalue = _score_test_from_covariance(score, Cov)

    score_full = np.full(genotype.shape[1], np.nan, dtype=float)
    se_full = np.full(genotype.shape[1], np.nan, dtype=float)
    pvalue_full = np.full(genotype.shape[1], np.nan, dtype=float)

    score_full[rv_label] = score
    se_full[rv_label] = se
    pvalue_full[rv_label] = pvalue

    return {
        "RV_label": rv_label,
        "Score": score_full,
        "SE": se_full,
        "pvalue": pvalue_full,
    }


def indiv_score_test_region_cond(
    genotype: np.ndarray,
    genotype_adj: np.ndarray,
    obj_nullmodel,
    rare_maf_cutoff: float = 0.01,
    rv_num_cutoff: int = 2,
    rv_num_cutoff_max: float = 1e9,
    method_cond: str = "optimal",
    precomputed_cov_cond: np.ndarray | None = None,
) -> Dict[str, object]:
    genotype = _as_2d_matrix(genotype, "genotype")
    if genotype.shape[1] == 1:
        raise ValueError("Number of rare variant in the set is less than 2!")

    genotype_adj = _as_2d_matrix(genotype_adj, "genotype_adj")
    if genotype.shape[0] != genotype_adj.shape[0]:
        raise ValueError("Dimensions don't match for genotype and genotype_adj!")

    G_flip, _, maf = matrix_flip(genotype)
    rv_label = (maf < rare_maf_cutoff) & (maf > 0)
    G = G_flip[:, rv_label]

    if np.sum(rv_label) >= rv_num_cutoff_max:
        raise ValueError("Number of rare variant in the set is more than rv_num_cutoff_max!")
    if np.sum(rv_label) < rv_num_cutoff:
        raise ValueError("Number of rare variant in the set is less than rv_num_cutoff!")

    X_cond = _build_conditional_design(genotype_adj, obj_nullmodel.X, method_cond)

    if obj_nullmodel.relatedness:
        residuals = _residualize(obj_nullmodel.scaled_residuals, X_cond)
        if precomputed_cov_cond is not None:
            Cov = np.asarray(precomputed_cov_cond, dtype=float)
            if Cov.shape != (G.shape[1], G.shape[1]):
                raise ValueError("precomputed_cov_cond has incompatible dimensions.")
            Cov = 0.5 * (Cov + Cov.T)
        else:
            sigma_solver = getattr(obj_nullmodel, "sigma_solver", None)
            Sigma_iX = getattr(obj_nullmodel, "Sigma_iX", None)
            cov = getattr(obj_nullmodel, "cov", None)
            if sigma_solver is not None and Sigma_iX is not None and cov is not None:
                Sigma_iG = sigma_solver.solve(G)
                tSigma_iX_G = Sigma_iX.T @ G
                Cov = Sigma_iG.T @ G - tSigma_iX_G.T @ cov @ tSigma_iX_G

                Sigma_iX_adj = sigma_solver.solve(X_cond)
                PX_adj = Sigma_iX_adj - Sigma_iX @ cov @ (Sigma_iX.T @ X_cond)

                A = np.linalg.inv(X_cond.T @ X_cond)
                XtG = X_cond.T @ G
                PXtG = PX_adj.T @ G
                Cov = (
                    Cov
                    - XtG.T @ A @ PXtG
                    - PXtG.T @ A @ XtG
                    + XtG.T @ A @ (PX_adj.T @ X_cond) @ A @ XtG
                )
                Cov = 0.5 * (Cov + Cov.T)
            else:
                P = getattr(obj_nullmodel, "P", None)
                if P is None:
                    raise ValueError(
                        "Related null model is missing components for conditional score-test covariance."
                    )
                A = np.linalg.inv(X_cond.T @ X_cond)
                XtG = X_cond.T @ G
                PX_adj = P @ X_cond
                PXtG = PX_adj.T @ G
                PG = P @ G
                Cov = (
                    PG.T @ G
                    - XtG.T @ A @ PXtG
                    - PXtG.T @ A @ XtG
                    + XtG.T @ A @ (PX_adj.T @ X_cond) @ A @ XtG
                )
                Cov = 0.5 * (Cov + Cov.T)
    else:
        residuals = _residualize(obj_nullmodel.y - obj_nullmodel.fitted, X_cond)
        PG = _apply_projection_unrelated(
            M=G,
            X=obj_nullmodel.X,
            working=obj_nullmodel.weights,
            family=obj_nullmodel.family,
        )
        PX_adj = _apply_projection_unrelated(
            M=X_cond,
            X=obj_nullmodel.X,
            working=obj_nullmodel.weights,
            family=obj_nullmodel.family,
        )
        A = np.linalg.inv(X_cond.T @ X_cond)
        XtG = X_cond.T @ G
        PXtG = PX_adj.T @ G
        Cov = (
            PG.T @ G
            - XtG.T @ A @ PXtG
            - PXtG.T @ A @ XtG
            + XtG.T @ A @ (PX_adj.T @ X_cond) @ A @ XtG
        )
        Cov = 0.5 * (Cov + Cov.T)

    score = np.asarray(residuals @ G, dtype=float).reshape(-1)
    se, pvalue = _score_test_from_covariance(score, Cov)

    score_full = np.full(genotype.shape[1], np.nan, dtype=float)
    se_full = np.full(genotype.shape[1], np.nan, dtype=float)
    pvalue_full = np.full(genotype.shape[1], np.nan, dtype=float)

    score_full[rv_label] = score
    se_full[rv_label] = se
    pvalue_full[rv_label] = pvalue

    return {
        "RV_label": rv_label,
        "Score_cond": score_full,
        "SE_cond": se_full,
        "pvalue_cond": pvalue_full,
    }


def _apply_projection_unrelated(
    M: np.ndarray,
    X: np.ndarray,
    working: np.ndarray,
    family: str,
) -> np.ndarray:
    if family == "gaussian":
        return M - X @ np.linalg.inv(X.T @ X) @ (X.T @ M)
    if family == "binomial":
        WX = X * working[:, None]
        WM = working[:, None] * M
        return WM - WX @ np.linalg.inv(X.T @ WX) @ (X.T @ WM)
    raise NotImplementedError(f"Unsupported family for STAAR_cond: {family}")


def staar_cond(
    genotype: np.ndarray,
    genotype_adj: np.ndarray,
    obj_nullmodel,
    annotation_phred: np.ndarray | None = None,
    rare_maf_cutoff: float = 0.01,
    rv_num_cutoff: int = 2,
    rv_num_cutoff_max: float = 1e9,
    method_cond: str = "optimal",
    precomputed_cov_cond: np.ndarray | None = None,
) -> Dict[str, object]:
    genotype = _as_2d_matrix(genotype, "genotype")
    if genotype.shape[1] == 1:
        raise ValueError("Number of rare variant in the set is less than 2!")

    genotype_adj = _as_2d_matrix(genotype_adj, "genotype_adj")
    if genotype.shape[0] != genotype_adj.shape[0]:
        raise ValueError("Dimensions don't match for genotype and genotype_adj!")

    annotation_phred = _as_annotation_matrix(annotation_phred, genotype.shape[1])
    if annotation_phred.size > 0 and genotype.shape[1] != annotation_phred.shape[0]:
        raise ValueError("Dimensions don't match for genotype and annotation!")

    G_flip, _, maf = matrix_flip(genotype)
    rv_label = (maf < rare_maf_cutoff) & (maf > 0)
    G = G_flip[:, rv_label]
    annotation_phred = annotation_phred[rv_label, :]

    if np.sum(rv_label) >= rv_num_cutoff_max:
        raise ValueError("Number of rare variant in the set is more than rv_num_cutoff_max!")
    if np.sum(rv_label) < rv_num_cutoff:
        raise ValueError("Number of rare variant in the set is less than rv_num_cutoff!")

    maf = maf[rv_label]
    w_B, w_S, w_A = _compute_weights(maf, annotation_phred)

    X_cond = _build_conditional_design(genotype_adj, obj_nullmodel.X, method_cond)

    if obj_nullmodel.relatedness:
        residuals = _residualize(obj_nullmodel.scaled_residuals, X_cond)
        if precomputed_cov_cond is not None:
            Cov = np.asarray(precomputed_cov_cond, dtype=float)
            if Cov.shape != (G.shape[1], G.shape[1]):
                raise ValueError("precomputed_cov_cond has incompatible dimensions.")
        else:
            sigma_solver = obj_nullmodel.sigma_solver
            Sigma_iX = obj_nullmodel.Sigma_iX
            cov = obj_nullmodel.cov

            Sigma_iG = sigma_solver.solve(G)
            tSigma_iX_G = Sigma_iX.T @ G
            Cov = Sigma_iG.T @ G - tSigma_iX_G.T @ cov @ tSigma_iX_G

            Sigma_iX_adj = sigma_solver.solve(X_cond)
            PX_adj = Sigma_iX_adj - Sigma_iX @ cov @ (Sigma_iX.T @ X_cond)

            A = np.linalg.inv(X_cond.T @ X_cond)
            XtG = X_cond.T @ G
            PXtG = PX_adj.T @ G
            Cov = (
                Cov
                - XtG.T @ A @ PXtG
                - PXtG.T @ A @ XtG
                + XtG.T @ A @ (PX_adj.T @ X_cond) @ A @ XtG
            )
    else:
        residuals = _residualize(obj_nullmodel.y - obj_nullmodel.fitted, X_cond)
        PG = _apply_projection_unrelated(
            M=G,
            X=obj_nullmodel.X,
            working=obj_nullmodel.weights,
            family=obj_nullmodel.family,
        )
        PX_adj = _apply_projection_unrelated(
            M=X_cond,
            X=obj_nullmodel.X,
            working=obj_nullmodel.weights,
            family=obj_nullmodel.family,
        )
        A = np.linalg.inv(X_cond.T @ X_cond)
        XtG = X_cond.T @ G
        PXtG = PX_adj.T @ G
        Cov = (
            PG.T @ G
            - XtG.T @ A @ PXtG
            - PXtG.T @ A @ XtG
            + XtG.T @ A @ (PX_adj.T @ X_cond) @ A @ XtG
        )

    Cov = 0.5 * (Cov + Cov.T)
    pvalues = _staartest_pvalues_smmat(
        G=G,
        residuals=residuals,
        Cov=Cov,
        weights_B=w_B,
        weights_S=w_S,
        weights_A=w_A,
        mac=np.round(maf * 2 * G.shape[0]).astype(int),
    )

    num_variant = int(np.sum(rv_label))
    num_annotation = annotation_phred.shape[1] + 1

    results = {}
    results["num_variant"] = float(num_variant)
    results["cMAC"] = float(np.sum(G))
    results["RV_label"] = rv_label

    results["results_STAAR_O_cond"] = float(cct(pvalues))
    idx_acat_o = [
        0,
        num_annotation,
        2 * num_annotation,
        3 * num_annotation,
        4 * num_annotation,
        5 * num_annotation,
    ]
    results["results_ACAT_O_cond"] = float(cct(pvalues[idx_acat_o]))

    def _assemble(start_idx: int, prefix: str, label: str) -> Dict[str, float]:
        segment = pvalues[start_idx : start_idx + num_annotation]
        agg = cct(segment)
        values = list(segment) + [agg]
        cols = [
            f"{prefix}",
            *[f"{prefix}-{c}" for c in _annotation_columns(annotation_phred)],
            label,
        ]
        return dict(zip(cols, map(float, values)))

    results["results_STAAR_S_1_25_cond"] = _assemble(0, "SKAT(1,25)", "STAAR-S(1,25)")
    results["results_STAAR_S_1_1_cond"] = _assemble(num_annotation, "SKAT(1,1)", "STAAR-S(1,1)")
    results["results_STAAR_B_1_25_cond"] = _assemble(
        2 * num_annotation, "Burden(1,25)", "STAAR-B(1,25)"
    )
    results["results_STAAR_B_1_1_cond"] = _assemble(
        3 * num_annotation, "Burden(1,1)", "STAAR-B(1,1)"
    )
    results["results_STAAR_A_1_25_cond"] = _assemble(
        4 * num_annotation, "ACAT-V(1,25)", "STAAR-A(1,25)"
    )
    results["results_STAAR_A_1_1_cond"] = _assemble(
        5 * num_annotation, "ACAT-V(1,1)", "STAAR-A(1,1)"
    )

    return results


def staar(
    genotype: np.ndarray,
    obj_nullmodel,
    annotation_phred: np.ndarray,
    rare_maf_cutoff: float = 0.01,
    rv_num_cutoff: int = 2,
    rv_num_cutoff_max: float = 1e9,
) -> Dict[str, object]:
    if genotype.ndim != 2:
        raise ValueError("genotype is not a matrix!")
    if genotype.shape[1] == 1:
        raise ValueError("Number of rare variant in the set is less than 2!")

    annotation_phred = np.asarray(annotation_phred, dtype=float)
    if annotation_phred.size > 0 and genotype.shape[1] != annotation_phred.shape[0]:
        raise ValueError("Dimensions don't match for genotype and annotation!")

    G_flip, _, maf = matrix_flip(genotype)
    rv_label = (maf < rare_maf_cutoff) & (maf > 0)
    Geno_rare = G_flip[:, rv_label]

    annotation_phred = annotation_phred[rv_label, :]

    if np.sum(rv_label) >= rv_num_cutoff_max:
        raise ValueError("Number of rare variant in the set is more than rv_num_cutoff_max!")

    if np.sum(rv_label) < rv_num_cutoff:
        raise ValueError("Number of rare variant in the set is less than rv_num_cutoff!")

    maf = maf[rv_label]
    G = Geno_rare

    w_B, w_S, w_A = _compute_weights(maf, annotation_phred)

    if obj_nullmodel.relatedness:
        # Use relatedness path with Sigma_i solver
        sigma_solver = obj_nullmodel.sigma_solver
        Sigma_iX = obj_nullmodel.Sigma_iX
        cov = obj_nullmodel.cov
        residuals = obj_nullmodel.scaled_residuals

        if obj_nullmodel.precomputed_cov is None:
            Sigma_iG = sigma_solver.solve(G)
            tSigma_iX_G = Sigma_iX.T @ G
            Cov = Sigma_iG.T @ G - tSigma_iX_G.T @ cov @ tSigma_iX_G
        else:
            Cov = obj_nullmodel.precomputed_cov
        pvalues = _staartest_pvalues(
            G=G,
            residuals=residuals,
            Cov=Cov,
            weights_B=w_B,
            weights_S=w_S,
            weights_A=w_A,
            mac=np.round(maf * 2 * G.shape[0]).astype(int),
            sigma=None,
            fam=1,
        )
    else:
        X = obj_nullmodel.X
        working = obj_nullmodel.weights
        sigma = math.sqrt(obj_nullmodel.dispersion)
        fam = 0 if obj_nullmodel.family == "gaussian" else 1
        residuals = obj_nullmodel.y - obj_nullmodel.fitted
        df_resid = X.shape[0] - X.shape[1]

        if fam == 0:
            tX_G = X.T @ G
            Cov = G.T @ G - tX_G.T @ np.linalg.inv(X.T @ X) @ tX_G
        else:
            WX = X * working[:, None]
            tX_G = X.T @ (working[:, None] * G)
            Cov = (working[:, None] * G).T @ G - tX_G.T @ np.linalg.inv(X.T @ WX) @ tX_G

        pvalues = _staartest_pvalues(
            G=G,
            residuals=residuals,
            Cov=Cov,
            weights_B=w_B,
            weights_S=w_S,
            weights_A=w_A,
            mac=np.round(maf * 2 * G.shape[0]).astype(int),
            sigma=sigma,
            fam=fam,
            df_resid=df_resid,
        )

    num_variant = int(np.sum(rv_label))
    num_annotation = annotation_phred.shape[1] + 1

    results = {}
    results["num_variant"] = float(num_variant)

    results["results_STAAR_O"] = float(cct(pvalues))
    idx_acat_o = [
        0,
        num_annotation,
        2 * num_annotation,
        3 * num_annotation,
        4 * num_annotation,
        5 * num_annotation,
    ]
    results["results_ACAT_O"] = float(cct(pvalues[idx_acat_o]))

    def _assemble(start_idx: int, prefix: str, label: str) -> Dict[str, float]:
        segment = pvalues[start_idx : start_idx + num_annotation]
        agg = cct(segment)
        values = list(segment) + [agg]
        cols = [
            f"{prefix}",
            *[f"{prefix}-{c}" for c in _annotation_columns(annotation_phred)],
            label,
        ]
        return dict(zip(cols, map(float, values)))

    results["results_STAAR_S_1_25"] = _assemble(0, "SKAT(1,25)", "STAAR-S(1,25)")
    results["results_STAAR_S_1_1"] = _assemble(num_annotation, "SKAT(1,1)", "STAAR-S(1,1)")

    results["results_STAAR_B_1_25"] = _assemble(
        2 * num_annotation, "Burden(1,25)", "STAAR-B(1,25)"
    )
    results["results_STAAR_B_1_1"] = _assemble(
        3 * num_annotation, "Burden(1,1)", "STAAR-B(1,1)"
    )

    results["results_STAAR_A_1_25"] = _assemble(
        4 * num_annotation, "ACAT-V(1,25)", "STAAR-A(1,25)"
    )
    results["results_STAAR_A_1_1"] = _assemble(
        5 * num_annotation, "ACAT-V(1,1)", "STAAR-A(1,1)"
    )

    return results


def _ordered_unique(values: np.ndarray) -> List[object]:
    seen = set()
    out: List[object] = []
    for val in values.tolist():
        if val not in seen:
            seen.add(val)
            out.append(val)
    return out


def _variant_set_pvalues(
    G: np.ndarray,
    maf: np.ndarray,
    annotation_phred: np.ndarray,
    obj_nullmodel,
    precomputed_cov: np.ndarray | None = None,
    precomputed_weights: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
    precomputed_mac: np.ndarray | None = None,
) -> np.ndarray:
    if precomputed_weights is None:
        w_B, w_S, w_A = _compute_weights(maf, annotation_phred)
    else:
        w_B, w_S, w_A = precomputed_weights
    if precomputed_mac is None:
        mac = np.round(maf * 2 * G.shape[0]).astype(int)
    else:
        mac = precomputed_mac

    if obj_nullmodel.relatedness:
        residuals = obj_nullmodel.scaled_residuals
        if precomputed_cov is not None:
            Cov = np.asarray(precomputed_cov, dtype=float)
            if Cov.shape != (G.shape[1], G.shape[1]):
                raise ValueError(
                    "precomputed_cov shape does not match the rare-variant dimension."
                )
        else:
            sigma_solver = getattr(obj_nullmodel, "sigma_solver", None)
            Sigma_iX = getattr(obj_nullmodel, "Sigma_iX", None)
            cov = getattr(obj_nullmodel, "cov", None)
            if sigma_solver is not None and Sigma_iX is not None and cov is not None:
                Sigma_iG = sigma_solver.solve(G)
                tSigma_iX_G = Sigma_iX.T @ G
                Cov = Sigma_iG.T @ G - tSigma_iX_G.T @ cov @ tSigma_iX_G
            else:
                Sigma_i = getattr(obj_nullmodel, "Sigma_i", None)
                if Sigma_i is not None and Sigma_iX is not None and cov is not None:
                    Sigma_iG = Sigma_i @ G
                    tSigma_iX_G = Sigma_iX.T @ G
                    Cov = Sigma_iG.T @ G - tSigma_iX_G.T @ cov @ tSigma_iX_G
                else:
                    P = getattr(obj_nullmodel, "P", None)
                    if P is None:
                        raise ValueError(
                            "Related null model is missing covariance components for AI-STAAR."
                        )
                    Cov = (P @ G).T @ G
        Cov = 0.5 * (Cov + Cov.T)
        return _staartest_pvalues(
            G=G,
            residuals=residuals,
            Cov=Cov,
            weights_B=w_B,
            weights_S=w_S,
            weights_A=w_A,
            mac=mac,
            sigma=None,
            fam=1,
        )

    X = obj_nullmodel.X
    working = obj_nullmodel.weights
    sigma = math.sqrt(obj_nullmodel.dispersion)
    fam = 0 if obj_nullmodel.family == "gaussian" else 1
    residuals = obj_nullmodel.y - obj_nullmodel.fitted
    df_resid = X.shape[0] - X.shape[1]
    if precomputed_cov is not None:
        Cov = np.asarray(precomputed_cov, dtype=float)
        if Cov.shape != (G.shape[1], G.shape[1]):
            raise ValueError("precomputed_cov shape does not match the rare-variant dimension.")
        Cov = 0.5 * (Cov + Cov.T)
    else:
        if fam == 0:
            tX_G = X.T @ G
            Cov = G.T @ G - tX_G.T @ np.linalg.inv(X.T @ X) @ tX_G
        else:
            WX = X * working[:, None]
            tX_G = X.T @ (working[:, None] * G)
            Cov = (working[:, None] * G).T @ G - tX_G.T @ np.linalg.inv(X.T @ WX) @ tX_G
        Cov = 0.5 * (Cov + Cov.T)
    return _staartest_pvalues(
        G=G,
        residuals=residuals,
        Cov=Cov,
        weights_B=w_B,
        weights_S=w_S,
        weights_A=w_A,
        mac=mac,
        sigma=sigma,
        fam=fam,
        df_resid=df_resid,
    )


def _ai_unrelated_gaussian_cov_from_group_stats(
    pop_weights: np.ndarray,
    group_gram: np.ndarray,
    group_xtg: np.ndarray,
    inv_xtx: np.ndarray,
) -> np.ndarray:
    pop_weights = np.asarray(pop_weights, dtype=float)
    cov = np.tensordot(pop_weights * pop_weights, group_gram, axes=(0, 0))
    xtg = np.tensordot(pop_weights, group_xtg, axes=(0, 0))
    cov = cov - xtg.T @ inv_xtx @ xtg
    return 0.5 * (cov + cov.T)


def _assemble_ai_outputs_from_pvalues(
    pvalues: np.ndarray,
    annotation_phred: np.ndarray,
) -> Dict[str, object]:
    num_annotation = annotation_phred.shape[1] + 1

    results: Dict[str, object] = {}
    results["results_STAAR_O"] = float(cct(pvalues))
    idx_acat_o = [
        0,
        num_annotation,
        2 * num_annotation,
        3 * num_annotation,
        4 * num_annotation,
        5 * num_annotation,
    ]
    results["results_ACAT_O"] = float(cct(pvalues[idx_acat_o]))

    def _assemble(start_idx: int, prefix: str, label: str) -> Dict[str, float]:
        segment = pvalues[start_idx : start_idx + num_annotation]
        agg = cct(segment)
        values = list(segment) + [agg]
        cols = [
            f"{prefix}",
            *[f"{prefix}-{c}" for c in _annotation_columns(annotation_phred)],
            label,
        ]
        return dict(zip(cols, map(float, values)))

    results["results_STAAR_S_1_25"] = _assemble(0, "SKAT(1,25)", "STAAR-S(1,25)")
    results["results_STAAR_S_1_1"] = _assemble(num_annotation, "SKAT(1,1)", "STAAR-S(1,1)")
    results["results_STAAR_B_1_25"] = _assemble(
        2 * num_annotation, "Burden(1,25)", "STAAR-B(1,25)"
    )
    results["results_STAAR_B_1_1"] = _assemble(
        3 * num_annotation, "Burden(1,1)", "STAAR-B(1,1)"
    )
    results["results_STAAR_A_1_25"] = _assemble(
        4 * num_annotation, "ACAT-V(1,25)", "STAAR-A(1,25)"
    )
    results["results_STAAR_A_1_1"] = _assemble(
        5 * num_annotation, "ACAT-V(1,1)", "STAAR-A(1,1)"
    )
    return results


def ai_staar(
    genotype: np.ndarray,
    obj_nullmodel,
    annotation_phred: np.ndarray | None = None,
    rare_maf_cutoff: float = 0.01,
    rv_num_cutoff: int = 2,
    rv_num_cutoff_max: float = 1e9,
    find_weight: bool = False,
) -> Dict[str, object]:
    genotype = _as_2d_matrix(genotype, "genotype")
    if genotype.shape[1] == 1:
        raise ValueError("Number of rare variant in the set is less than 2!")

    annotation_phred = _as_annotation_matrix(annotation_phred, genotype.shape[1])
    if annotation_phred.size > 0 and genotype.shape[1] != annotation_phred.shape[0]:
        raise ValueError("Dimensions don't match for genotype and annotation!")

    pop_groups = getattr(obj_nullmodel, "pop_groups", None)
    pop_weights_1_1 = getattr(obj_nullmodel, "pop_weights_1_1", None)
    pop_weights_1_25 = getattr(obj_nullmodel, "pop_weights_1_25", None)
    if pop_groups is None or pop_weights_1_1 is None or pop_weights_1_25 is None:
        raise ValueError(
            "AI-STAAR requires pop_groups, pop_weights_1_1, and pop_weights_1_25 on obj_nullmodel."
        )

    pop_groups = np.asarray(pop_groups)
    pop_weights_1_1 = np.asarray(pop_weights_1_1, dtype=float)
    pop_weights_1_25 = np.asarray(pop_weights_1_25, dtype=float)

    G_flip, _, maf_full = matrix_flip(genotype)
    rv_label = (maf_full < rare_maf_cutoff) & (maf_full > 0)
    G_base = G_flip[:, rv_label]
    annotation_phred = annotation_phred[rv_label, :]

    if np.sum(rv_label) >= rv_num_cutoff_max:
        raise ValueError("Number of rare variant in the set is more than rv_num_cutoff_max!")
    if np.sum(rv_label) < rv_num_cutoff:
        raise ValueError("Number of rare variant in the set is less than rv_num_cutoff!")

    if pop_groups.shape[0] != G_base.shape[0]:
        raise ValueError("pop_groups length does not match sample size.")

    pop_levels = getattr(obj_nullmodel, "pop_levels", None)
    if pop_levels is None:
        pop_levels = _ordered_unique(pop_groups)
    else:
        pop_levels = [lvl for lvl in list(pop_levels)]
    n_pop = len(pop_levels)
    if pop_weights_1_1.shape[0] != n_pop or pop_weights_1_25.shape[0] != n_pop:
        raise ValueError("AI-STAAR pop-weight rows must match the number of population levels.")
    if pop_weights_1_1.shape[1] != pop_weights_1_25.shape[1]:
        raise ValueError("AI-STAAR pop-weight matrices must have the same number of base tests.")
    B = pop_weights_1_1.shape[1]

    level_to_idx = {level: idx for idx, level in enumerate(pop_levels)}
    pop_idx = np.fromiter(
        (level_to_idx.get(level, -1) for level in pop_groups.tolist()),
        dtype=np.int64,
        count=pop_groups.shape[0],
    )
    if np.any(pop_idx < 0):
        raise ValueError("pop_groups contains levels that are not in pop_levels.")

    indices = [np.where(pop_idx == i)[0] for i in range(n_pop)]
    a_p = np.zeros(n_pop, dtype=float)
    for i, idx in enumerate(indices):
        if idx.size == 0:
            a_p[i] = 0.0
            continue
        g_sub = G_base[idx, :]
        maf_sub = np.minimum(np.mean(g_sub, axis=0) / 2.0, 1.0 - np.mean(g_sub, axis=0) / 2.0)
        a_mean = float(np.mean(maf_sub))
        a_p[i] = stats.beta.pdf(a_mean, 1, 25) if a_mean > 0 else 0.0

    maf = maf_full[rv_label]
    precomputed_weights = _compute_weights(maf, annotation_phred)
    precomputed_mac = np.round(maf * 2 * G_base.shape[0]).astype(int)
    precomputed_cov_s1 = getattr(obj_nullmodel, "precomputed_ai_cov_s1", None)
    precomputed_cov_s2 = getattr(obj_nullmodel, "precomputed_ai_cov_s2", None)

    group_gram = None
    group_xtg = None
    inv_xtx = None
    if (not obj_nullmodel.relatedness) and getattr(obj_nullmodel, "family", None) == "gaussian":
        X = obj_nullmodel.X
        inv_xtx = np.linalg.inv(X.T @ X)
        group_gram = np.stack(
            [G_base[idx, :].T @ G_base[idx, :] for idx in indices],
            axis=0,
        )
        group_xtg = np.stack(
            [X[idx, :].T @ G_base[idx, :] for idx in indices],
            axis=0,
        )

    pvalues_1_tot: list[np.ndarray] = []
    pvalues_2_tot: list[np.ndarray] = []
    weight_all_1: list[np.ndarray] = []
    weight_all_2: list[np.ndarray] = []

    for b in range(B):
        w_b_1 = pop_weights_1_1[:, b]
        w_b_2 = a_p * pop_weights_1_25[:, b]
        if find_weight:
            weight_all_1.append(w_b_1.copy())
            weight_all_2.append(w_b_2.copy())

        G1 = G_base * w_b_1[pop_idx][:, None]
        G2 = G_base * w_b_2[pop_idx][:, None]

        precomputed_cov_1 = None
        precomputed_cov_2 = None
        if isinstance(precomputed_cov_s1, (list, tuple)) and b < len(precomputed_cov_s1):
            candidate = np.asarray(precomputed_cov_s1[b], dtype=float)
            if candidate.shape == (maf.size, maf.size):
                precomputed_cov_1 = candidate
        if isinstance(precomputed_cov_s2, (list, tuple)) and b < len(precomputed_cov_s2):
            candidate = np.asarray(precomputed_cov_s2[b], dtype=float)
            if candidate.shape == (maf.size, maf.size):
                precomputed_cov_2 = candidate
        if precomputed_cov_1 is None and group_gram is not None and group_xtg is not None and inv_xtx is not None:
            precomputed_cov_1 = _ai_unrelated_gaussian_cov_from_group_stats(
                pop_weights=w_b_1,
                group_gram=group_gram,
                group_xtg=group_xtg,
                inv_xtx=inv_xtx,
            )
        if precomputed_cov_2 is None and group_gram is not None and group_xtg is not None and inv_xtx is not None:
            precomputed_cov_2 = _ai_unrelated_gaussian_cov_from_group_stats(
                pop_weights=w_b_2,
                group_gram=group_gram,
                group_xtg=group_xtg,
                inv_xtx=inv_xtx,
            )

        pvalues_1 = _variant_set_pvalues(
            G=G1,
            maf=maf,
            annotation_phred=annotation_phred,
            obj_nullmodel=obj_nullmodel,
            precomputed_cov=precomputed_cov_1,
            precomputed_weights=precomputed_weights,
            precomputed_mac=precomputed_mac,
        )
        pvalues_2 = _variant_set_pvalues(
            G=G2,
            maf=maf,
            annotation_phred=annotation_phred,
            obj_nullmodel=obj_nullmodel,
            precomputed_cov=precomputed_cov_2,
            precomputed_weights=precomputed_weights,
            precomputed_mac=precomputed_mac,
        )
        pvalues_1_tot.append(pvalues_1)
        pvalues_2_tot.append(pvalues_2)

    pvalues_1_tot_arr = np.column_stack(pvalues_1_tot)
    pvalues_2_tot_arr = np.column_stack(pvalues_2_tot)
    pvalues_tot = np.column_stack([pvalues_1_tot_arr, pvalues_2_tot_arr])
    pvalues_aggregate = np.array([cct(row) for row in pvalues_tot], dtype=float)

    num_variant = int(np.sum(rv_label))
    cmac = float(np.sum(G_base))

    results: Dict[str, object] = {}
    results["num_variant"] = float(num_variant)
    results["cMAC"] = cmac
    results["RV_label"] = rv_label
    results.update(_assemble_ai_outputs_from_pvalues(pvalues_aggregate, annotation_phred))

    if find_weight:
        if weight_all_1:
            results["weight_all_1"] = np.column_stack(weight_all_1)
            results["weight_all_2"] = np.column_stack(weight_all_2)
        else:
            results["weight_all_1"] = np.zeros((n_pop, 0), dtype=float)
            results["weight_all_2"] = np.zeros((n_pop, 0), dtype=float)

        num_tests = pvalues_1_tot_arr.shape[0]
        results_weight = []
        results_weight1 = []
        results_weight2 = []
        for b in range(B):
            pvalues_aggregate_weight_b = np.array(
                [
                    cct(np.array([pvalues_1_tot_arr[row_idx, b], pvalues_2_tot_arr[row_idx, b]]))
                    for row_idx in range(num_tests)
                ],
                dtype=float,
            )
            payload_prefix = {"num_variant": float(num_variant), "cMAC": cmac}
            results_weight.append(
                {
                    **payload_prefix,
                    **_assemble_ai_outputs_from_pvalues(
                        pvalues_aggregate_weight_b,
                        annotation_phred,
                    ),
                }
            )
            results_weight1.append(
                {
                    **payload_prefix,
                    **_assemble_ai_outputs_from_pvalues(
                        pvalues_1_tot_arr[:, b],
                        annotation_phred,
                    ),
                }
            )
            results_weight2.append(
                {
                    **payload_prefix,
                    **_assemble_ai_outputs_from_pvalues(
                        pvalues_2_tot_arr[:, b],
                        annotation_phred,
                    ),
                }
            )
        results["results_weight"] = results_weight
        results["results_weight1"] = results_weight1
        results["results_weight2"] = results_weight2

    return results


def _annotation_columns(annotation_phred: np.ndarray) -> List[str]:
    if annotation_phred.shape[1] == 0:
        return []
    # Default column names for example data
    return [f"Z{i}" for i in range(1, annotation_phred.shape[1] + 1)]
