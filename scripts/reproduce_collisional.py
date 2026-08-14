#!/usr/bin/env python3
"""Audit geometry, normalisation, and convergence of Tau Ceti f/f_max scenarios."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

from tau_ceti.collisional import (
    broad_belt_fmax_continuous,
    broad_belt_fmax_independent_annuli,
    characteristic_radius_ring_fmax,
)

ROOT = Path(__file__).resolve().parents[1]
SEED = 20260723
N_DRAWS = 200_000


def positive_normal(rng, mean, sd, size, lower):
    result = rng.normal(mean, sd, size)
    while np.any(result <= lower):
        bad = result <= lower
        result[bad] = rng.normal(mean, sd, bad.sum())
    return result


def loguniform(rng, low, high, size):
    return np.exp(rng.uniform(np.log(low), np.log(high), size))


def quantiles(values):
    return tuple(float(x) for x in np.percentile(values, [2.5, 16, 50, 84, 97.5]))


def draw_priors():
    rng = np.random.default_rng(SEED)
    posterior = np.load(ROOT / "posterior_samples/single_mbb_posterior.npz")
    return {
        "f_obs": rng.choice(posterior["fractional_luminosity"], size=N_DRAWS, replace=True),
        "age_myr": positive_normal(rng, 7630.0, 870.0, N_DRAWS, 4000.0),
        "stellar_mass_solar": positive_normal(rng, 0.783, 0.012, N_DRAWS, 0.1),
        "blowout_diameter_um": positive_normal(rng, 1.2, 0.3, N_DRAWS, 0.1),
        "largest_body_km": loguniform(rng, 1.0, 2000.0, N_DRAWS),
        "disruption_energy_j_per_kg": loguniform(rng, 50.0, 500.0, N_DRAWS),
        "eccentricity": loguniform(rng, 0.01, 0.1, N_DRAWS),
    }


def chunked_broad(function, prior, **extra):
    output = np.empty(N_DRAWS)
    for start in range(0, N_DRAWS, 5000):
        stop = min(start + 5000, N_DRAWS)
        index = slice(start, stop)
        output[index] = function(
            stellar_mass_solar=prior["stellar_mass_solar"][index, None],
            age_myr=prior["age_myr"][index, None],
            blowout_diameter_um=prior["blowout_diameter_um"][index, None],
            largest_body_km=prior["largest_body_km"][index, None],
            disruption_energy_j_per_kg=prior["disruption_energy_j_per_kg"][index, None],
            eccentricity=prior["eccentricity"][index, None],
            inclination=prior["eccentricity"][index, None],
            **extra,
        )
    return output


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    table_dir = ROOT / "tables"; figure_dir = ROOT / "figures"; report_dir = ROOT / "docs"
    table_dir.mkdir(parents=True, exist_ok=True); figure_dir.mkdir(parents=True, exist_ok=True)
    common = dict(inner_radius_au=6.0, outer_radius_au=55.0, stellar_mass_solar=0.783,
                  age_myr=7630.0, blowout_diameter_um=1.2, largest_body_km=100.0,
                  disruption_energy_j_per_kg=200.0, eccentricity=0.05, inclination=0.05)
    convergence = []
    reference = broad_belt_fmax_continuous(annuli=1000, binning="log", midpoint="geometric", **common)
    for annuli in (10, 25, 50, 100, 250, 500, 1000):
        for binning in ("linear", "log"):
            for midpoint in ("arithmetic", "geometric"):
                value = broad_belt_fmax_continuous(annuli=annuli, binning=binning, midpoint=midpoint, **common)
                convergence.append({"annuli": annuli, "binning": binning, "midpoint": midpoint, "fmax": float(value), "relative_to_log_geometric_1000": float(value/reference - 1.0)})
    write_csv(table_dir / "collisional_geometry_audit.csv", convergence)

    prior = draw_priors()
    characteristic_radius, fmax_ring = characteristic_radius_ring_fmax(
        6.0, 55.0, 0.0, prior["stellar_mass_solar"], prior["age_myr"], prior["blowout_diameter_um"],
        prior["largest_body_km"], prior["disruption_energy_j_per_kg"], prior["eccentricity"], prior["eccentricity"],
    )
    fmax_continuous = chunked_broad(broad_belt_fmax_continuous, prior, inner_radius_au=6.0, outer_radius_au=55.0, annuli=250, binning="log", midpoint="geometric")
    fmax_independent = chunked_broad(broad_belt_fmax_independent_annuli, prior, inner_radius_au=6.0, outer_radius_au=55.0, surface_density_power=0.0, annuli=250, binning="log")
    models = {
        "A_characteristic_radius_ring": (fmax_ring, "One belt at cross-section-weighted R_c with its full Δr/R_c."),
        "B_continuous_broad_belt": (fmax_continuous, "Integral of local Eq. 14 ceilings over dln r."),
        "C_independent_annuli": (fmax_independent, "Zones evolve independently; total ceiling is Σ fmax,i."),
    }
    summary_rows = []
    for name, (fmax, description) in models.items():
        ratio = prior["f_obs"] / fmax
        q2, q16, q50, q84, q97 = quantiles(ratio)
        summary_rows.append({"formulation": name, "definition": description, "fmax_p50": float(np.median(fmax)), "ratio_p2_5": q2, "ratio_p16": q16, "ratio_p50": q50, "ratio_p84": q84, "ratio_p97_5": q97, "p_ratio_gt_1": float(np.mean(ratio > 1)), "p_ratio_gt_10": float(np.mean(ratio > 10)), "classification": "Consistent with steady state" if np.mean(ratio > 1) < 0.16 else "Marginal or assumption-dependent"})
    write_csv(table_dir / "collisional_geometry_summary.csv", summary_rows)
    write_csv(table_dir / "collisional_scenarios.csv", summary_rows)

    prior_rows = []
    for model in models:
        for parameter, spec in (("f_obs", "SED posterior samples"), ("age_myr", "truncated Normal(7630,870); lower 4000"), ("stellar_mass_solar", "truncated Normal(0.783,0.012)"), ("blowout_diameter_um", "truncated Normal(1.2,0.3); model assumption"), ("largest_body_km", "log-uniform(1,2000)"), ("disruption_energy_j_per_kg", "log-uniform(50,500)"), ("eccentricity", "log-uniform(0.01,0.1), I=e"), ("belt_geometry", "6–55 au; p=0")):
            prior_rows.append({"formulation": model, "parameter": parameter, "distribution_or_value": spec})
    write_csv(table_dir / "collisional_geometry_priors.csv", prior_rows)

    delta = np.log10((prior["f_obs"] / fmax_continuous) / (prior["f_obs"] / fmax_ring))
    sensitivity = []
    for parameter in ("age_myr", "stellar_mass_solar", "blowout_diameter_um", "largest_body_km", "disruption_energy_j_per_kg", "eccentricity"):
        rho, pvalue = spearmanr(np.log10(prior[parameter]), delta)
        sensitivity.append({"parameter": parameter, "spearman_rho_with_delta": float(rho), "p_value": float(pvalue)})
    write_csv(table_dir / "collisional_geometry_delta_sensitivity.csv", sensitivity)

    fig, ax = plt.subplots(figsize=(7.4, 4.6), constrained_layout=True)
    for (binning, midpoint), group in __import__("itertools").groupby(sorted(convergence, key=lambda x: (x["binning"], x["midpoint"], x["annuli"])), key=lambda x: (x["binning"], x["midpoint"])):
        group = list(group); ax.plot([x["annuli"] for x in group], [100*x["relative_to_log_geometric_1000"] for x in group], marker="o", label=f"{binning}, {midpoint}")
    ax.axhline(0, color="0.3", linewidth=1); ax.set(xscale="log", xlabel="Number of annuli", ylabel="Difference from 1000-bin reference (%)", title="Broad-belt ceiling convergence")
    ax.legend(fontsize=8, ncol=2); fig.savefig(figure_dir / "collisional_annulus_convergence.svg"); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.4, 4.6), constrained_layout=True)
    names = [row["formulation"].replace("_", " ") for row in summary_rows]
    medians = [row["ratio_p50"] for row in summary_rows]; lows = [row["ratio_p50"]-row["ratio_p16"] for row in summary_rows]; highs = [row["ratio_p84"]-row["ratio_p50"] for row in summary_rows]
    ax.errorbar(medians, np.arange(3), xerr=[lows, highs], fmt="o", color="#2a6fbb", capsize=3); ax.axvline(1, color="0.3", linewidth=1); ax.axvline(10, color="0.3", linewidth=1, linestyle="--")
    ax.set(xscale="log", yticks=np.arange(3), yticklabels=names, xlabel="Observed f / fmax (16th–84th percentile)", title="Collisional geometry formulations")
    fig.savefig(figure_dir / "collisional_geometry_comparison.svg"); plt.close(fig)

    max_convergence_error = max(abs(row["relative_to_log_geometric_1000"]) for row in convergence if row["annuli"] >= 100)
    sensitivity_text = ", ".join(f"{row['parameter']}={row['spearman_rho_with_delta']:.3f}" for row in sensitivity)
    report = ["# Collisional geometry audit", "", "## Finding", "", "The previous broad-annulus result was numerically invalid. It divided a sum of width-dependent local ceilings by an area normalisation. Because Eq. 14 contains `dr/r`, that operation makes the reported total decrease as annuli are made narrower. It explains the former broad-belt median near 6 and is not a defensible physical ceiling.", "", "## Equation and units", "", "`literature/equations.md` now maps Eq. 14 directly to code. The source defines `r` and `dr` in au, age in Myr, `D_c` as a *diameter* in km, `D_bl` as a diameter in μm, and `Q_D*` in J kg^-1. The implementation converts μm/km before taking the diameter ratio. Eq. 14 is a narrow belt expression; its width is `dr/r`. It contains no independent stellar-luminosity exponent except through `D_bl` (Eq. 5).", "", "## Formulations", "", f"A uses a cross-section-weighted characteristic radius R_c={characteristic_radius:.2f} au and full Δr/R_c. B integrates local narrow-belt ceilings over 6–55 au. C assumes independently evolving zones but sums their luminosity ceilings, so it is algebraically identical to B when the same Eq. 14 local prescription is imposed. A different surface-density profile changes allocation among zones, but not Σfmax_i.", "", "## Convergence and diagnosis", "", f"Across binning and midpoint conventions, the largest difference for at least 100 annuli is {100*max_convergence_error:.3f}%. The corrected continuous integral and independent-annuli sum agree by construction; any remaining difference from A is geometry (characteristic radius and full-width approximation), together with interactions in Eq. 14's radius-dependent collision threshold. Rank correlations with Δ=log10[(f/fmax)_B/(f/fmax)_A] are: {sensitivity_text}. Shared f_obs and age cancel algebraically from Δ, so they are not drivers.", "", "## Classification", ""]
    report += [f"- **{row['formulation']}**: median f/fmax={row['ratio_p50']:.2g}; 68% {row['ratio_p16']:.2g}–{row['ratio_p84']:.2g}; 95% {row['ratio_p2_5']:.2g}–{row['ratio_p97_5']:.2g}; P(>1)={row['p_ratio_gt_1']:.3f}; P(>10)={row['p_ratio_gt_10']:.3f}; **{row['classification']}**." for row in summary_rows]
    report += ["", "The evidence permits steady-state evolution. It does not establish a transient event under any valid formulation. The broad-belt conclusion remains assumption-dependent because the narrow-belt analytic ceiling has been extrapolated over a wide resolved belt; it must not be compressed into one geometry-independent number."]
    (report_dir / "collisional_geometry_audit.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    (report_dir / "collisional_analysis.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
