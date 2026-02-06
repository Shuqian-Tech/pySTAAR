import pystaar
from pystaar import models, staar_core, staar_stats


def test_r_namespace_aliases_are_exposed():
    assert pystaar.CCT is staar_stats.cct
    assert pystaar.fit_null_glm is models.fit_null_glm
    assert pystaar.fit_null_glmmkin is models.fit_null_glmmkin
    assert pystaar.fit_null_glm_Binary_SPA is models.fit_null_glm_binary_spa
    assert pystaar.fit_null_glmmkin_Binary_SPA is models.fit_null_glmmkin_binary_spa
    assert pystaar.STAAR is staar_core.staar
    assert pystaar.STAAR_cond is staar_core.staar_cond
    assert pystaar.Indiv_Score_Test_Region is staar_core.indiv_score_test_region
    assert pystaar.Indiv_Score_Test_Region_cond is staar_core.indiv_score_test_region_cond
    assert pystaar.STAAR_Binary_SPA is staar_core.staar_binary_spa
    assert pystaar.AI_STAAR is staar_core.ai_staar


def test_r_namespace_aliases_are_listed_in_all():
    required = {
        "CCT",
        "fit_null_glm",
        "fit_null_glmmkin",
        "fit_null_glm_Binary_SPA",
        "fit_null_glmmkin_Binary_SPA",
        "STAAR",
        "STAAR_cond",
        "Indiv_Score_Test_Region",
        "Indiv_Score_Test_Region_cond",
        "STAAR_Binary_SPA",
        "AI_STAAR",
    }
    assert required.issubset(set(pystaar.__all__))


def test_cct_alias_smoke():
    pvalue = pystaar.CCT([0.1, 0.2, 0.8])
    assert 0.0 <= pvalue <= 1.0
