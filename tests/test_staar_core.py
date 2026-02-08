import numpy as np

from pystaar import staar_core


def test_matrix_flip_all_missing_columns_match_legacy_behavior():
    genotype = np.full((4, 3), -1.0, dtype=float)

    flipped, af, maf = staar_core.matrix_flip(genotype)

    np.testing.assert_allclose(flipped, np.zeros_like(flipped))
    assert np.isnan(af).all()
    assert np.isnan(maf).all()


def _clear_eigensolver_caches():
    for fn_name in (
        "_numpy_config_text",
        "_blas_backend_hint",
        "_eigensolver_size_threshold",
        "_selected_eigensolver",
    ):
        fn = getattr(staar_core, fn_name, None)
        if fn is not None and hasattr(fn, "cache_clear"):
            fn.cache_clear()


def test_eigensolver_policy_env_override(monkeypatch):
    monkeypatch.setenv("PYSTAAR_EIGENSOLVER", "scipy")
    _clear_eigensolver_caches()
    assert staar_core._selected_eigensolver() == "scipy"

    monkeypatch.setenv("PYSTAAR_EIGENSOLVER", "numpy")
    _clear_eigensolver_caches()
    assert staar_core._selected_eigensolver() == "numpy"


def test_eigensolver_backend_autoselect(monkeypatch):
    monkeypatch.delenv("PYSTAAR_EIGENSOLVER", raising=False)
    monkeypatch.setattr(staar_core, "_blas_backend_hint", lambda: "openblas")
    _clear_eigensolver_caches()
    assert staar_core._selected_eigensolver() == "scipy"

    monkeypatch.setattr(staar_core, "_blas_backend_hint", lambda: "accelerate")
    _clear_eigensolver_caches()
    assert staar_core._selected_eigensolver() == "numpy"


def test_eigensolver_small_matrix_prefers_numpy(monkeypatch):
    monkeypatch.delenv("PYSTAAR_EIGENSOLVER", raising=False)
    monkeypatch.setenv("PYSTAAR_EIGENSOLVER_SIZE_THRESHOLD", "256")
    monkeypatch.setattr(staar_core, "_selected_eigensolver", lambda: "scipy")
    _clear_eigensolver_caches()

    calls = {"numpy": 0, "scipy": 0}
    original_numpy_eigvalsh = staar_core.np.linalg.eigvalsh

    def fake_numpy(A):
        calls["numpy"] += 1
        return original_numpy_eigvalsh(np.asarray(A, dtype=float))

    def fake_scipy(A):
        calls["scipy"] += 1
        return original_numpy_eigvalsh(np.asarray(A, dtype=float))

    monkeypatch.setattr(staar_core.np.linalg, "eigvalsh", fake_numpy)
    monkeypatch.setattr(staar_core, "_eigvalsh_scipy", fake_scipy)

    mat = np.eye(8, dtype=float)
    vals = staar_core._eigvalsh_symmetric(mat)
    assert np.all(np.isfinite(vals))
    assert calls["numpy"] == 1
    assert calls["scipy"] == 0
