# Reproducibility protocol

1. Create the clean Conda environment from `environment.yml`.
2. Run `make verify`; this checks the accepted values, posterior archives and grain samples without requiring telescope archive access.
3. Run `make reproduce`; this reruns the seeded 200,000-draw collisional geometry analysis.
4. Compare the regenerated `tables/collisional_geometry_summary.csv` with the committed release table.

The exact collisional priors are in `config/collisional_reference.yaml`. The grain-grid construction and continuous-draw weights are in `config/grain_scenarios.yaml`. Population MCMC settings are in `config/population_mcmc.yaml`; the NetCDF archives are validation sensitivities.

The public DOI is intentionally blank until a Zenodo record is actually published. Do not replace it with a placeholder DOI.
