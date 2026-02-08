# Release Readiness Smoke Check

- Generated: 2026-02-08T08:46:02Z
- Platform: Darwin 24.6.0 (arm64)
- Python: 3.13.5
- Overall status: PASS

## Steps

- `PASS` `build_wheel` in `9.53s`: /opt/anaconda3/bin/python -m pip wheel . --no-deps -w /var/folders/4g/l2tqzrz960vgmb3nq4t9xxlh0000gn/T/pystaar_release_smoke_g_s8zo96/wheelhouse
- `PASS` `upgrade_pip` in `5.35s`: /var/folders/4g/l2tqzrz960vgmb3nq4t9xxlh0000gn/T/pystaar_release_smoke_g_s8zo96/venv/bin/python -m pip install --upgrade pip
- `PASS` `install_wheel` in `26.14s`: /var/folders/4g/l2tqzrz960vgmb3nq4t9xxlh0000gn/T/pystaar_release_smoke_g_s8zo96/venv/bin/python -m pip install /var/folders/4g/l2tqzrz960vgmb3nq4t9xxlh0000gn/T/pystaar_release_smoke_g_s8zo96/wheelhouse/pystaar-1.0.0-py3-none-any.whl
- `PASS` `import_and_workflow_smoke` in `29.67s`: /var/folders/4g/l2tqzrz960vgmb3nq4t9xxlh0000gn/T/pystaar_release_smoke_g_s8zo96/venv/bin/python -c import pystaar
info = pystaar.get_runtime_cache_info()
assert "workflows" in info and "data" in info
result = pystaar.staar_unrelated_glm(dataset="example", seed=600, rare_maf_cutoff=0.05)
assert "results_STAAR_O" in result
assert 0.0 <= float(result["results_STAAR_O"]) <= 1.0
after = pystaar.clear_runtime_caches(include_dataset_cache=True)
assert after["after"]["workflows"]["related_nullmodel"]["currsize"] == 0
print("SMOKE_OK", float(result["results_STAAR_O"]))

## Artifacts

- JSON details: `reports/release_smoke.json`
