# pySTAAR Quickstart (English)

`pySTAAR` is the Python migration of the STAAR R package.

For Chinese users, start here first: [`README_CN.md`](README_CN.md).

## Install

```bash
pip install -e .
```

## Minimal example

```python
from pystaar import staar_unrelated_glm

res = staar_unrelated_glm(dataset="example", seed=600, rare_maf_cutoff=0.05)
print(res["results_STAAR_O"])
```

## Main workflow APIs

- `staar_unrelated_glm`, `staar_related_sparse_glmmkin`, `staar_related_dense_glmmkin`
- `staar_unrelated_glm_cond`, `staar_related_sparse_glmmkin_cond`, `staar_related_dense_glmmkin_cond`
- `staar_unrelated_binary_spa`, `staar_related_sparse_binary_spa`, `staar_related_dense_binary_spa`
- `ai_staar_unrelated_glm`, `ai_staar_related_sparse_glmmkin`, `ai_staar_related_dense_glmmkin`

## More docs

- Installation: [`installation.md`](installation.md)
- Performance comparison (Python vs R): [`performance_comparison.md`](performance_comparison.md)
- Tutorials: [`tutorials/`](tutorials/)
- API output fields: [`api/output_fields.md`](api/output_fields.md)
- API stability policy: [`api/stability.md`](api/stability.md)
- R migration guide: [`migration_from_r.md`](migration_from_r.md)
- R quick migration checklist (CN): [`migration_r_quickstart_cn.md`](migration_r_quickstart_cn.md)
- Dataset directory template (CN): [`data_directory_template_cn.md`](data_directory_template_cn.md)
- Changelog: [`../CHANGELOG.md`](../CHANGELOG.md)

## Parity note

Functional migration scope is complete through `STAAR-56`. Related-workflow parity runs on pure Python paths; current parity and deviation history are tracked in `reports/summary.md` and `reports/deviations.md`.
