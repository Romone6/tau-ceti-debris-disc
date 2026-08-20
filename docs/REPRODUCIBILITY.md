# Reproducibility protocol

1. Create the clean Conda environment from `environment.yml`.
2. Run `make verify`; this checks the accepted values, posterior archives and grain samples without requiring telescope archive access.
3. Run `make reproduce`; this reruns the seeded 200,000-draw collisional geometry analysis.
4. Compare the regenerated `tables/collisional_geometry_summary.csv` with the committed release table.

The exact collisional priors are in `config/collisional_reference.yaml`. The grain-grid construction and continuous-draw weights are in `config/grain_scenarios.yaml`. Population MCMC settings are in `config/population_mcmc.yaml`; the NetCDF archives are validation sensitivities.

The published Zenodo archive DOI for the v1.0.2 release is [10.5281/zenodo.22021388](https://doi.org/10.5281/zenodo.22021388). The GitHub release and Zenodo record are immutable references for the frozen scientific product.
