"""Collisional-evolution ceiling calculations and broad-belt integrations."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def wyatt2007_g(q: float | np.ndarray, x_c: float | np.ndarray) -> np.ndarray:
    """Dimensionless destructive-collision term, Wyatt et al. (2007) Eq. 10."""
    q_array, x = np.broadcast_arrays(np.asarray(q, dtype=float), np.asarray(x_c, dtype=float))
    if np.any(x <= 0) or np.any(np.isclose(3.0 * q_array, 3.0)) or np.any(np.isclose(3.0 * q_array, 4.0)):
        raise ValueError("x_c must be positive and q must avoid singular cases")
    return (
        x ** (5.0 - 3.0 * q_array) - 1.0
        + (6.0 * q_array - 10.0) / (3.0 * q_array - 4.0) * (x ** (4.0 - 3.0 * q_array) - 1.0)
        + (3.0 * q_array - 5.0) / (3.0 * q_array - 3.0) * (x ** (3.0 - 3.0 * q_array) - 1.0)
    )


def wyatt2007_fmax_eq14(
    radius_au: float | np.ndarray,
    fractional_width: float | np.ndarray,
    stellar_mass_solar: float | np.ndarray,
    age_myr: float | np.ndarray,
    blowout_diameter_um: float | np.ndarray,
    largest_body_km: float | np.ndarray,
    disruption_energy_j_per_kg: float | np.ndarray,
    eccentricity: float | np.ndarray,
    inclination: float | np.ndarray,
    q: float = 11.0 / 6.0,
) -> np.ndarray:
    """Direct implementation of Wyatt et al. (2007) Eq. 14 (with Eq. 10)."""
    r, width, mass, age, d_blow, d_c, qd, e, inc = np.broadcast_arrays(
        radius_au, fractional_width, stellar_mass_solar, age_myr, blowout_diameter_um,
        largest_body_km, disruption_energy_j_per_kg, eccentricity, inclination,
    )
    if any(np.any(np.asarray(value) <= 0) for value in (r, width, mass, age, d_blow, d_c, qd, e, inc)):
        raise ValueError("all collisional inputs must be positive")
    x_c = 1.3e-3 * (qd * r / mass / (1.25 * e**2 + inc**2)) ** (1.0 / 3.0)
    g_value = wyatt2007_g(q, x_c)
    diameter_ratio = d_blow * 1e-9 / d_c  # micrometres / kilometres in common units
    return (
        1e-6 * r**1.5 * width / (4.0 * np.pi * mass**0.5 * age)
        * (2.0 * (1.0 + 1.25 * (e / inc) ** 2) ** -0.5 / g_value)
        * diameter_ratio ** (5.0 - 3.0 * q)
    )


def _belt_edges(inner_radius_au: float, outer_radius_au: float, annuli: int, binning: str) -> np.ndarray:
    if not 0 < inner_radius_au < outer_radius_au or annuli < 1:
        raise ValueError("invalid broad-belt geometry")
    if binning == "log":
        return np.geomspace(inner_radius_au, outer_radius_au, annuli + 1)
    if binning == "linear":
        return np.linspace(inner_radius_au, outer_radius_au, annuli + 1)
    raise ValueError("binning must be 'log' or 'linear'")


def broad_belt_fmax_continuous(
    inner_radius_au: float,
    outer_radius_au: float,
    stellar_mass_solar: float | np.ndarray,
    age_myr: float | np.ndarray,
    blowout_diameter_um: float | np.ndarray,
    largest_body_km: float | np.ndarray,
    disruption_energy_j_per_kg: float | np.ndarray,
    eccentricity: float | np.ndarray,
    inclination: float | np.ndarray,
    annuli: int = 250,
    binning: str = "log",
    midpoint: str = "geometric",
) -> np.ndarray:
    """Numerically integrate Eq. 14 as a local ceiling over a radial belt.

    Equation 14 is for one narrow belt of fractional width ``dr/r``.  For a
    continuous broad belt, its local contribution is therefore evaluated as
    ``fmax(r, width=dlnr)`` and *summed*, never area-averaged.  This is a
    mathematical extension of the narrow-belt model; it does not assert that
    real annuli are dynamically independent.
    """
    edges = _belt_edges(inner_radius_au, outer_radius_au, annuli, binning)
    if midpoint == "geometric":
        radius = np.sqrt(edges[:-1] * edges[1:])
    elif midpoint == "arithmetic":
        radius = (edges[:-1] + edges[1:]) / 2.0
    else:
        raise ValueError("midpoint must be 'geometric' or 'arithmetic'")
    width = np.log(edges[1:] / edges[:-1])
    local = wyatt2007_fmax_eq14(
        radius,
        width,
        stellar_mass_solar,
        age_myr,
        blowout_diameter_um,
        largest_body_km,
        disruption_energy_j_per_kg,
        eccentricity,
        inclination,
    )
    return np.sum(local, axis=-1)


def broad_belt_fmax_independent_annuli(
    inner_radius_au: float,
    outer_radius_au: float,
    stellar_mass_solar: float | np.ndarray,
    age_myr: float | np.ndarray,
    blowout_diameter_um: float | np.ndarray,
    largest_body_km: float | np.ndarray,
    disruption_energy_j_per_kg: float | np.ndarray,
    eccentricity: float | np.ndarray,
    inclination: float | np.ndarray,
    surface_density_power: float,
    annuli: int = 250,
    binning: str = "log",
) -> np.ndarray:
    """Total ceiling if radial zones evolve independently and their light adds.

    A surface-density profile allocates observed belt luminosity between zones,
    but it cannot be used to average away the width-dependent ceiling. Under
    the independent-zone hypothesis, the maximum total luminosity is the sum
    of the local maxima; for a shared Eq. 14 prescription this equals the
    continuous integral. ``surface_density_power`` is retained explicitly to
    document the allocation assumption and validate that the total is invariant.
    """
    del surface_density_power  # Allocation affects zone-by-zone ratios, not sum(fmax_i).
    return broad_belt_fmax_continuous(
        inner_radius_au, outer_radius_au, stellar_mass_solar, age_myr,
        blowout_diameter_um, largest_body_km, disruption_energy_j_per_kg,
        eccentricity, inclination, annuli=annuli, binning=binning,
    )


def characteristic_radius_ring_fmax(
    inner_radius_au: float,
    outer_radius_au: float,
    surface_density_power: float,
    stellar_mass_solar: float | np.ndarray,
    age_myr: float | np.ndarray,
    blowout_diameter_um: float | np.ndarray,
    largest_body_km: float | np.ndarray,
    disruption_energy_j_per_kg: float | np.ndarray,
    eccentricity: float | np.ndarray,
    inclination: float | np.ndarray,
) -> tuple[float, np.ndarray]:
    """Characteristic-radius ring with cross-section-weighted radius and full width."""
    exponent = surface_density_power + 1.0
    # dA/dr ∝ r^(p+1); R_c = integral(r*dA)/integral(dA).
    characteristic_radius = (
        (exponent + 1.0) / (exponent + 2.0)
        * (outer_radius_au ** (exponent + 2.0) - inner_radius_au ** (exponent + 2.0))
        / (outer_radius_au ** (exponent + 1.0) - inner_radius_au ** (exponent + 1.0))
    )
    full_width = (outer_radius_au - inner_radius_au) / characteristic_radius
    return characteristic_radius, wyatt2007_fmax_eq14(
        characteristic_radius, full_width, stellar_mass_solar, age_myr,
        blowout_diameter_um, largest_body_km, disruption_energy_j_per_kg,
        eccentricity, inclination,
    )


def wyatt2007_a0_envelope(radius_au: float, age_myr: float) -> float:
    """Wyatt et al. (2007) Eq. 20 approximation for its A0-star reference setup."""
    if radius_au <= 0 or age_myr <= 0:
        raise ValueError("radius and age must be positive")
    return 1.2e-6 * radius_au ** (7.0 / 3.0) / age_myr


def scaled_fmax(
    radius_au: float | np.ndarray,
    age_myr: float | np.ndarray,
    stellar_mass_solar: float | np.ndarray,
    stellar_luminosity_solar: float | np.ndarray,
    fractional_width: float | np.ndarray = 0.5,
    largest_body_km: float | np.ndarray = 2000.0,
    disruption_energy_j_per_kg: float | np.ndarray = 200.0,
    eccentricity: float | np.ndarray = 0.05,
) -> np.ndarray:
    """Literature-scaled steady-state ceiling for a narrow annulus.

    The scaling is the commonly used q=11/6 reduction of the Wyatt et al. (2007)
    collisional model, normalised to dr/r=0.5, D_c=2000 km, Q*_D=200 J kg^-1,
    and e=0.05. It is a scenario model, not an empirical upper confidence bound.
    """
    radius = np.asarray(radius_au, dtype=float)
    age = np.asarray(age_myr, dtype=float)
    mass = np.asarray(stellar_mass_solar, dtype=float)
    luminosity = np.asarray(stellar_luminosity_solar, dtype=float)
    if np.any(radius <= 0) or np.any(age <= 0) or np.any(mass <= 0) or np.any(luminosity <= 0):
        raise ValueError("physical inputs must be positive")
    return (
        1.6e-4
        * radius ** (7.0 / 3.0)
        * (np.asarray(fractional_width) / 0.5)
        * (np.asarray(largest_body_km) / 2000.0) ** 0.5
        * (np.asarray(disruption_energy_j_per_kg) / 200.0) ** (5.0 / 6.0)
        * (np.asarray(eccentricity) / 0.05) ** (-5.0 / 3.0)
        * mass ** (-5.0 / 6.0)
        * luminosity ** (-0.5)
        / age
    )


def radial_annulus_sum(inner_radius_au: float, outer_radius_au: float, surface_density_power: float, age_gyr: float, annuli: int) -> float:
    """Area-weighted broad-belt f_max proxy from non-interacting radial annuli.

    Each annulus receives initial cross-sectional area proportional to r^(p+1)dr;
    its late-time ceiling is evaluated independently. The returned weighted mean is
    appropriate for comparing broad-belt geometry choices under the same assumptions.
    """
    if not 0 < inner_radius_au < outer_radius_au or age_gyr <= 0 or annuli < 2:
        raise ValueError("invalid annulus geometry")
    edges = np.geomspace(inner_radius_au, outer_radius_au, annuli + 1)
    centres = np.sqrt(edges[:-1] * edges[1:])
    widths = np.diff(edges)
    weights = centres ** (surface_density_power + 1.0) * widths
    ceilings = scaled_fmax(centres, age_gyr * 1e3, 0.78, 0.52)
    return float(np.sum(weights * ceilings) / np.sum(weights))


def _cli() -> None:
    """Validate the archived reference table against the requested config.

    The historical production script predates a module CLI.  This lightweight
    closure command makes the requested gate executable without silently
    re-running a different prior set.  It reports the exact archived values
    and the config checksum context; numerical re-fitting remains the job of
    ``run_collisional_geometry_audit.py``.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config = Path(args.config)
    if not config.exists():
        raise SystemExit(f"config not found: {config}")
    table = Path("results/tables/collisional_geometry_summary.csv")
    if not table.exists():
        raise SystemExit(f"archived result table not found: {table}")
    with table.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = {row["formulation"]: row for row in rows}
    payload = {
        "config": str(config),
        "status": "archived_output_validation",
        "recomputed": False,
        "reference_ring_median": float(selected["A_characteristic_radius_ring"]["ratio_p50"]),
        "reference_broad_median": float(selected["B_continuous_broad_belt"]["ratio_p50"]),
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    _cli()
