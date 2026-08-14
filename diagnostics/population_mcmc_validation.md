# Population MCMC validation

The frozen strict merged input contains 233 stars and 44 detections, with Tau Ceti excluded from the fit. A PyMC NUTS validation model was run with four independent chains, 1,500 tuning draws and 1,500 retained draws per chain. A separate DUNES-only sensitivity fit contains 126 stars and is archived alongside the merged fit.

| Check | Result |
|---|---:|
| Split rank-normalized R-hat | 1.000 for all reported parameters |
| Minimum bulk ESS | 2,300 |
| Minimum tail ESS | 2,800 |
| Divergences | 0 |
| Maximum tree depth | 6 |
| BFMI by chain | 0.968, 1.019, 0.991, 1.007 |

The ArviZ-compatible posterior is archived at `posterior_samples/population_strict_merged_idata.nc`; trace and rank plots are in `diagnostics/`. The resulting Bayesian validation sensitivity gives (q_\tau\) median approximately 0.933 and (z_\tau\) median approximately 1.50. This is close to, but not numerically identical to, the accepted MLE/bootstrap baseline (0.927 and 1.46), because it introduces explicit Bayesian priors and standardized predictors. It is archived as a validation sensitivity, not silently promoted as the headline result.

The DUNES-only ArviZ posterior is archived at `posterior_samples/population_dunes_only_idata.nc`, with diagnostics in `tables/population_dunes_only_sampling_diagnostics.csv`. It also has R-hat=1.000, minimum bulk ESS=2,600, minimum tail ESS=2,700, zero divergences, maximum tree depth 5 and BFMI 0.926--1.027 by chain. Both fits used the explicit three-limit CDF marginalisation documented in `config/population_mcmc.yaml`; no posterior is claimed to supersede the accepted baseline.
