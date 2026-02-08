"""Python migration of the STAAR R package.

Public API includes:
- R-compatible names (e.g., ``STAAR`` / ``AI_STAAR`` / ``CCT``).
- Python workflow entry points used by parity scenarios.
"""

from .models import (
    fit_null_glm,
    fit_null_glm_binary_spa,
    fit_null_glmmkin,
    fit_null_glmmkin_binary_spa,
)
from .staar_core import (
    _get_eigensolver_runtime_info,
    ai_staar,
    indiv_score_test_region,
    indiv_score_test_region_cond,
    matrix_flip,
    staar,
    staar_binary_spa,
    staar_cond,
)
from .staar_stats import cct
from .workflows import (
    ai_staar_related_dense_glmmkin,
    ai_staar_related_dense_glmmkin_find_weight,
    ai_staar_related_sparse_glmmkin,
    ai_staar_related_sparse_glmmkin_find_weight,
    ai_staar_unrelated_glm,
    ai_staar_unrelated_glm_find_weight,
    clear_runtime_caches,
    get_runtime_cache_info,
    indiv_score_related_dense_glmmkin,
    indiv_score_related_dense_glmmkin_cond,
    indiv_score_related_sparse_glmmkin,
    indiv_score_related_sparse_glmmkin_cond,
    indiv_score_unrelated_glm,
    indiv_score_unrelated_glm_cond,
    staar_related_dense_binary_spa,
    staar_related_dense_glmmkin,
    staar_related_dense_glmmkin_cond,
    staar_related_sparse_binary_spa,
    staar_related_sparse_glmmkin,
    staar_related_sparse_glmmkin_cond,
    staar_unrelated_binary_spa,
    staar_unrelated_glm,
    staar_unrelated_glm_cond,
)

# R-compatible aliases from STAAR NAMESPACE exports.
CCT = cct
fit_null_glm_Binary_SPA = fit_null_glm_binary_spa
fit_null_glmmkin_Binary_SPA = fit_null_glmmkin_binary_spa
STAAR = staar
STAAR_cond = staar_cond
Indiv_Score_Test_Region = indiv_score_test_region
Indiv_Score_Test_Region_cond = indiv_score_test_region_cond
STAAR_Binary_SPA = staar_binary_spa
AI_STAAR = ai_staar

__all__ = [
    "CCT",
    "fit_null_glm",
    "fit_null_glmmkin",
    "fit_null_glm_Binary_SPA",
    "fit_null_glmmkin_Binary_SPA",
    "matrix_flip",
    "STAAR",
    "STAAR_cond",
    "STAAR_Binary_SPA",
    "Indiv_Score_Test_Region",
    "Indiv_Score_Test_Region_cond",
    "AI_STAAR",
    "staar_unrelated_glm",
    "staar_related_sparse_glmmkin",
    "staar_related_dense_glmmkin",
    "staar_unrelated_glm_cond",
    "staar_related_sparse_glmmkin_cond",
    "staar_related_dense_glmmkin_cond",
    "staar_unrelated_binary_spa",
    "staar_related_sparse_binary_spa",
    "staar_related_dense_binary_spa",
    "indiv_score_unrelated_glm",
    "indiv_score_related_sparse_glmmkin",
    "indiv_score_related_dense_glmmkin",
    "indiv_score_unrelated_glm_cond",
    "indiv_score_related_sparse_glmmkin_cond",
    "indiv_score_related_dense_glmmkin_cond",
    "ai_staar_unrelated_glm",
    "ai_staar_related_sparse_glmmkin",
    "ai_staar_related_dense_glmmkin",
    "ai_staar_unrelated_glm_find_weight",
    "ai_staar_related_sparse_glmmkin_find_weight",
    "ai_staar_related_dense_glmmkin_find_weight",
    "get_runtime_cache_info",
    "clear_runtime_caches",
    "get_eigensolver_runtime_info",
]


def get_eigensolver_runtime_info():
    """Return backend-aware eigensolver runtime selection metadata."""
    return _get_eigensolver_runtime_info()
