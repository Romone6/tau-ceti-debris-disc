# Data provenance

Every headline result is linked to a release path:

| Result | Configuration | Inputs | Output |
|---|---|---|---|
| Accepted population headline | archived baseline manifest | processed FGK/DEBRIS tables | `tables/reproduction_comparison.csv` |
| Collisional ratios | `config/collisional_reference.yaml` | SED posterior and stellar priors | `tables/collisional_geometry_summary.csv` |
| Grain scenarios | `config/grain_scenarios.yaml` | stellar inputs and composition grid | `results/grains/` |
| Bayesian sensitivity | `config/population_mcmc.yaml` | strict merged catalogue | `posterior_samples/` |

Overlapping survey identities and exclusions are retained in `data/processed/overlap_resolution.csv` and `data/processed/exclusion_log.csv`.
