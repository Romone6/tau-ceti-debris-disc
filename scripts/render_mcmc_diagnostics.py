from pathlib import Path
import arviz as az
import matplotlib.pyplot as plt

F = Path(__file__).resolve().parent
idata = az.from_netcdf(F / "posterior_samples/population_strict_merged_idata.nc")
vars_ = ["alpha", "beta_age", "beta_mass", "survey_intercept", "sigma_dunes", "sigma_debris"]
az.plot_trace(idata, var_names=vars_)
plt.savefig(F / "diagnostics/population_traceplots.pdf", bbox_inches="tight")
plt.close("all")
az.plot_rank(idata, var_names=vars_)
plt.savefig(F / "diagnostics/population_rank_plots.pdf", bbox_inches="tight")
plt.close("all")
print("diagnostic plots written")
