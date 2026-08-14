"""Small, dependency-free physical calculations used by the analysis."""

from __future__ import annotations

import math

AU_M = 1.495978707e11
C_LIGHT = 299_792_458.0
G_GRAV = 6.67430e-11
H_PLANCK = 6.62607015e-34
JY_W_M2_HZ = 1e-26
K_BOLTZMANN = 1.380649e-23
M_EARTH_KG = 5.9722e24
M_SUN_KG = 1.98847e30
L_SUN_W = 3.828e26
PC_M = 3.085677581491367e16


def planck_nu(frequency_hz: float, temperature_k: float) -> float:
    """Return the Planck specific intensity B_nu in W m^-2 Hz^-1 sr^-1."""
    if frequency_hz <= 0 or temperature_k <= 0:
        raise ValueError("frequency and temperature must be positive")
    exponent = H_PLANCK * frequency_hz / (K_BOLTZMANN * temperature_k)
    return 2.0 * H_PLANCK * frequency_hz**3 / C_LIGHT**2 / math.expm1(exponent)


def blackbody_temperature(radius_au: float, luminosity_solar: float) -> float:
    """Blackbody equilibrium temperature for a rapidly reradiating grain."""
    if radius_au <= 0 or luminosity_solar <= 0:
        raise ValueError("radius and luminosity must be positive")
    return 278.3 * luminosity_solar**0.25 / math.sqrt(radius_au)


def blackbody_radius_au(temperature_k: float, luminosity_solar: float) -> float:
    """Orbital radius of a blackbody grain at a specified equilibrium temperature."""
    if temperature_k <= 0 or luminosity_solar <= 0:
        raise ValueError("temperature and luminosity must be positive")
    return (278.3 * luminosity_solar**0.25 / temperature_k) ** 2


def orbital_period_years(semimajor_axis_au: float, stellar_mass_solar: float) -> float:
    """Keplerian orbital period, using AU, solar masses, and Julian years."""
    if semimajor_axis_au <= 0 or stellar_mass_solar <= 0:
        raise ValueError("semimajor axis and stellar mass must be positive")
    return math.sqrt(semimajor_axis_au**3 / stellar_mass_solar)


def dust_mass_from_flux(
    flux_jy: float,
    distance_pc: float,
    wavelength_m: float,
    temperature_k: float,
    opacity_m2_per_kg: float,
) -> float:
    """Optically thin dust mass, M=F_nu d^2/(kappa_nu B_nu)."""
    if min(flux_jy, distance_pc, wavelength_m, temperature_k, opacity_m2_per_kg) <= 0:
        raise ValueError("all dust-mass inputs must be positive")
    frequency_hz = C_LIGHT / wavelength_m
    return flux_jy * JY_W_M2_HZ * (distance_pc * PC_M) ** 2 / (
        opacity_m2_per_kg * planck_nu(frequency_hz, temperature_k)
    )


def blowout_radius_microns(
    luminosity_solar: float,
    stellar_mass_solar: float,
    density_kg_m3: float = 2500.0,
    radiation_pressure_efficiency: float = 1.0,
) -> float:
    """Grain radius at beta=0.5 for a spherical compact grain."""
    numerator = 3.0 * luminosity_solar * L_SUN_W * radiation_pressure_efficiency
    denominator = 8.0 * math.pi * C_LIGHT * G_GRAV * stellar_mass_solar * M_SUN_KG * density_kg_m3
    return numerator / denominator * 1e6
