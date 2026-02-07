# pySTAAR

Python migration of the STAAR R package.

Primary audience: Chinese genomics/statistical genetics users migrating from R STAAR workflows.

## Documentation

- Chinese quickstart: [`docs/README_CN.md`](docs/README_CN.md)
- English quickstart: [`docs/README.md`](docs/README.md)
- Installation: [`docs/installation.md`](docs/installation.md)
- Tutorials:
  - [`docs/tutorials/01_basic_staar.md`](docs/tutorials/01_basic_staar.md)
  - [`docs/tutorials/02_binary_spa.md`](docs/tutorials/02_binary_spa.md)
  - [`docs/tutorials/03_related_samples.md`](docs/tutorials/03_related_samples.md)
  - [`docs/tutorials/04_conditional.md`](docs/tutorials/04_conditional.md)
  - [`docs/tutorials/05_ai_staar.md`](docs/tutorials/05_ai_staar.md)
- API references:
  - [`docs/api/null_models.md`](docs/api/null_models.md)
  - [`docs/api/staar_functions.md`](docs/api/staar_functions.md)
  - [`docs/api/output_fields.md`](docs/api/output_fields.md)
  - [`docs/api/utilities.md`](docs/api/utilities.md)
- R migration guide: [`docs/migration_from_r.md`](docs/migration_from_r.md)
- Changelog: [`CHANGELOG.md`](CHANGELOG.md)

## Quick Install

```bash
pip install -e .
```

## Quick Run

```python
from pystaar import staar_unrelated_glm

res = staar_unrelated_glm(dataset="example", seed=600, rare_maf_cutoff=0.05)
print(res["results_STAAR_O"])
```

## Current Parity Policy

- Functional migration scope is complete (`STAAR-1` through `STAAR-41`).
- Related-workflow parity runs on pure Python paths.
- Approved parity tolerance relaxation exists for affected related sentinels (up to `rtol=5e-4`), tracked in `reports/deviations.md` (`DEV-001`).
