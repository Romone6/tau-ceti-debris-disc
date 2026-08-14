"""Four-chain PyMC validation of the frozen strict merged population model."""
from __future__ import annotations

import ast
import csv
import json
from pathlib import Path

import arviz as az
import numpy as np
import pymc as pm
from scipy.special import ndtr

ROOT = Path(__file__).resolve().parents[1]
F = ROOT / "finalisation"
SEEDS = [20260811, 20260812, 20260813, 20260814]


def read_rows():
    rows = []
    dunes_ids = set()
    dunes_path = ROOT / "data/processed/fgk_fractional_luminosities.csv"
    with dunes_path.open(encoding="utf-8", newline="") as f:
        dunes_ids = {f"HIP {r['hip_id']}" for r in csv.DictReader(f)}
    for path in [ROOT / "data/processed/fgk_fractional_luminosities.csv", ROOT / "data/processed/debris_fractional_luminosities.csv"]:
        survey_default = "DUNES" if "fgk_" in path.name else "DEBRIS"
        with path.open(encoding="utf-8", newline="") as f:
            for r in csv.DictReader(f):
                detected = (r.get("catalogue_excess", "") == "1") if survey_default == "DUNES" else (r.get("detected", "False").lower() == "true")
                if survey_default == "DEBRIS" and r.get("canonical_star_id", "") in dunes_ids:
                    continue
                age = float(r["age_gyr"]); mass = float(r.get("mass_solar", "nan"))
                if not np.isfinite(mass): mass = float(r.get("luminosity_solar", "nan")) ** 0.25
                if detected:
                    logf = float(r["logf"]); sigma = float(r["logf_sigma"])
                    limits = [logf - 10, logf - 10, logf - 10]
                else:
                    logf = np.nan; sigma = np.nan
                    if survey_default == "DUNES": limits = [float(r["logf_limit_p16"]), float(r["logf_limit_p50"]), float(r["logf_limit_p84"])]
                    else: limits = [float(x) for x in ast.literal_eval(r["limits"].replace("nan", "float('nan')"))]
                rows.append({"survey": survey_default, "age": age, "mass": mass, "detected": detected, "logf": logf, "sigma": sigma, "limits": limits})
    # Exact HIP overlap resolution: DUNES is primary. The DEBRIS table was
    # already filtered by run_debris_comparison into the strict non-overlap set.
    return rows


def fit(rows, out_path: Path, dunes_only: bool = False):
    detected = np.array([r["detected"] for r in rows], dtype=bool)
    age = np.log10(np.array([r["age"] for r in rows], dtype=float))
    mass = np.log10(np.array([r["mass"] for r in rows], dtype=float))
    survey = np.array([r["survey"] == "DEBRIS" for r in rows], dtype=bool)
    d_age = (age - np.mean(age)) / np.std(age)
    d_mass = (mass - np.mean(mass)) / np.std(mass)
    yi = np.flatnonzero(detected); ni = np.flatnonzero(~detected)
    y = np.array([r["logf"] for r in rows], dtype=float)
    ysig = np.array([r["sigma"] for r in rows], dtype=float)
    lim = np.array([r["limits"] for r in rows], dtype=float)
    with pm.Model() as model:
        alpha = pm.Normal("alpha", -5.0, 2.0)
        beta_age = pm.Normal("beta_age", -1.0, 1.0)
        beta_mass = pm.Normal("beta_mass", 0.0, 1.0)
        if dunes_only:
            sigma_dunes = pm.HalfNormal("sigma_dunes", 1.0)
            parameter_names = ["alpha", "beta_age", "beta_mass", "sigma_dunes"]
            mu = alpha + beta_age * d_age + beta_mass * d_mass
            sig = pm.math.ones_like(mu) * sigma_dunes
        else:
            survey_intercept = pm.Normal("survey_intercept", 0.0, 1.0)
            sigma_dunes = pm.HalfNormal("sigma_dunes", 1.0)
            sigma_debris = pm.HalfNormal("sigma_debris", 1.0)
            parameter_names = ["alpha", "beta_age", "beta_mass", "survey_intercept", "sigma_dunes", "sigma_debris"]
            mu = alpha + beta_age * d_age + beta_mass * d_mass + survey_intercept * survey
            sig = pm.math.switch(survey, sigma_debris, sigma_dunes)
        pm.Normal("detected", mu=mu[yi], sigma=pm.math.sqrt(sig[yi] ** 2 + ysig[yi] ** 2), observed=y[yi])
        lim_dist = pm.Normal.dist(mu=mu[ni, None], sigma=sig[ni, None])
        pm.Potential("censored", pm.math.logsumexp(pm.logcdf(lim_dist, lim[ni]), axis=1).sum() - len(ni) * np.log(3.0))
        idata = pm.sample(draws=1500, tune=1500, chains=4, cores=4, target_accept=0.92, random_seed=SEEDS, return_inferencedata=True, progressbar=False)
    idata.to_netcdf(out_path, engine="h5netcdf")
    return idata, d_age, d_mass, parameter_names


def main():
    rows = read_rows()
    out = F / "posterior_samples/population_strict_merged_idata.nc"
    idata, d_age, d_mass, parameter_names = fit(rows, out)
    summary = az.summary(idata, var_names=parameter_names, round_to=None)
    summary.to_csv(F / "tables/population_posterior_summary.csv")
    diag = az.summary(idata, var_names=parameter_names, round_to=None)
    rows_diag = []
    for name, row in diag.iterrows():
        rows_diag.append({"parameter": name, "r_hat": float(row["r_hat"]), "ess_bulk": float(row["ess_bulk"]), "ess_tail": float(row["ess_tail"]), "mean": float(row["mean"]), "sd": float(row["sd"])})
    write = lambda p, rs: (p.parent.mkdir(parents=True, exist_ok=True), p.write_text("parameter,r_hat,ess_bulk,ess_tail,mean,sd\n" + "\n".join(",".join(str(x[k]) for k in ["parameter", "r_hat", "ess_bulk", "ess_tail", "mean", "sd"]) for x in rs) + "\n"))
    write(F / "tables/population_sampling_diagnostics.csv", rows_diag)
    posterior = idata.posterior
    # Predictive draws conditioned on the accepted SED posterior and Tau age/mass.
    rng = np.random.default_rng(20260815)
    flat = {name: np.asarray(posterior[name]).reshape(-1) for name in ["alpha", "beta_age", "beta_mass", "sigma_dunes"]}
    n = len(flat["alpha"])
    tau_f = np.load(ROOT / "results/sed/single_mbb_posterior.npz")["fractional_luminosity"]
    take = rng.integers(0, len(tau_f), n)
    tau_age = np.maximum(rng.normal(7.63, 0.87, n), 4.0)
    tau_mass = np.full(n, 0.78)
    mu = flat["alpha"] + flat["beta_age"] * ((np.log10(tau_age) - np.mean(np.array([r["age"] for r in rows])) * 0) / 1.0) + flat["beta_mass"] * ((np.log10(tau_mass) - np.mean(np.log10([r["mass"] for r in rows]))) / np.std(np.log10([r["mass"] for r in rows])))
    # Recompute the same centring used in fit().
    ages = np.log10(np.array([r["age"] for r in rows])); masses = np.log10(np.array([r["mass"] for r in rows]))
    mu = flat["alpha"] + flat["beta_age"] * ((np.log10(tau_age) - ages.mean()) / ages.std()) + flat["beta_mass"] * ((np.log10(tau_mass) - masses.mean()) / masses.std())
    z = (np.log10(tau_f[take]) - mu) / flat["sigma_dunes"]
    q = ndtr(z)
    np.savez(F / "posterior_samples/tau_ceti_reference_prediction.npz", q_tau=q, z_tau=z, tau_f=tau_f[take], tau_age=tau_age)
    np.savez(F / "posterior_samples/tau_ceti_archival_prediction.npz", q_tau=q, z_tau=z, tau_f=np.full(n, 2.8513e-6), tau_age=tau_age)
    (F / "posterior_samples/population_model_metadata.json").write_text(json.dumps({"rows": len(rows), "detections": int(detected_count(rows)), "chains": 4, "draws": 1500, "tune": 1500, "target_accept": 0.92, "tau_holdout": True, "likelihood": "detections Normal with group scatter; non-detections averaged Normal CDF over three temperature-marginal limits"}, indent=2) + "\n")
    print("MCMC complete", len(rows), "rows", "q median", float(np.median(q)), "z median", float(np.median(z)))

    # A separately archived DUNES-only fit is required for survey sensitivity.
    dunes_rows = [r for r in rows if r["survey"] == "DUNES"]
    dunes_out = F / "posterior_samples/population_dunes_only_idata.nc"
    dunes_idata, _, _, dunes_names = fit(dunes_rows, dunes_out, dunes_only=True)
    dunes_summary = az.summary(dunes_idata, var_names=dunes_names, round_to=None)
    dunes_summary.to_csv(F / "tables/population_dunes_only_posterior_summary.csv")
    dunes_diag = []
    for name, row in dunes_summary.iterrows():
        dunes_diag.append({"parameter": name, "r_hat": float(row["r_hat"]), "ess_bulk": float(row["ess_bulk"]), "ess_tail": float(row["ess_tail"]), "mean": float(row["mean"]), "sd": float(row["sd"])})
    write(F / "tables/population_dunes_only_sampling_diagnostics.csv", dunes_diag)
    (F / "posterior_samples/population_dunes_only_metadata.json").write_text(json.dumps({"rows": len(dunes_rows), "detections": int(detected_count(dunes_rows)), "chains": 4, "draws": 1500, "tune": 1500, "target_accept": 0.92, "tau_holdout": True, "likelihood": "detections Normal with one intrinsic scatter; non-detections averaged Normal CDF over three temperature-marginal limits"}, indent=2) + "\n")
    print("DUNES-only MCMC complete", len(dunes_rows), "rows")


def detected_count(rows):
    return sum(r["detected"] for r in rows)


if __name__ == "__main__":
    main()
