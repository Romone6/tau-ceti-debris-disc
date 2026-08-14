"""SED components and likelihood functions for Tau Ceti's cold belt."""

from __future__ import annotations

import math

import numpy as np
from scipy.special import log_ndtr

from .physics import C_LIGHT, H_PLANCK, JY_W_M2_HZ, K_BOLTZMANN, L_SUN_W, PC_M


def modified_blackbody_emissivity(wavelength_um: float, lambda0_um: float, beta: float) -> float:
    """Return the unit-normalized long-wavelength modified-blackbody emissivity."""
    if min(wavelength_um, lambda0_um) <= 0 or beta < 0:
        raise ValueError("wavelengths must be positive and beta non-negative")
    return 1.0 if wavelength_um <= lambda0_um else (wavelength_um / lambda0_um) ** (-beta)


def greybody_shape(wavelength_um: np.ndarray | float, temperature_k: float, lambda0_um: float = 210.0, beta: float = 1.0) -> np.ndarray:
    """Unnormalised F_nu shape for a modified blackbody."""
    wavelengths = np.asarray(wavelength_um, dtype=float)
    if temperature_k <= 0:
        raise ValueError("temperature must be positive")
    frequency = C_LIGHT / (wavelengths * 1e-6)
    exponent = H_PLANCK * frequency / (K_BOLTZMANN * temperature_k)
    with np.errstate(over="ignore"):
        planck = 2.0 * H_PLANCK * frequency**3 / C_LIGHT**2 / np.expm1(exponent)
    emissivity = np.where(wavelengths <= lambda0_um, 1.0, (wavelengths / lambda0_um) ** (-beta))
    return planck * emissivity


def normalized_greybody(wavelength_um: np.ndarray | float, amplitude_mjy_at_160: float, temperature_k: float, lambda0_um: float = 210.0, beta: float = 1.0) -> np.ndarray:
    """Modified blackbody normalised to a dust-only 160-um flux density."""
    return amplitude_mjy_at_160 * greybody_shape(wavelength_um, temperature_k, lambda0_um, beta) / greybody_shape(160.0, temperature_k, lambda0_um, beta)


def gaussian_log_likelihood(observed_mjy: float, sigma_mjy: float, model_mjy: float) -> float:
    """Independent Gaussian log likelihood for a flux-density detection."""
    if sigma_mjy <= 0:
        raise ValueError("sigma must be positive")
    residual = (observed_mjy - model_mjy) / sigma_mjy
    return -0.5 * residual**2 - math.log(sigma_mjy * math.sqrt(2.0 * math.pi))


def censored_log_likelihood(limit_mjy: float, sigma_mjy: float, model_mjy: float) -> float:
    """Log P(F < limit | model, sigma), evaluated stably for Gaussian censoring."""
    if sigma_mjy <= 0:
        raise ValueError("sigma must be positive")
    return float(log_ndtr((limit_mjy - model_mjy) / sigma_mjy))


def dust_fractional_luminosity(
    wavelength_um: np.ndarray,
    dust_flux_mjy: np.ndarray,
    stellar_luminosity_solar: float,
    distance_pc: float,
) -> float:
    """Integrate an F_nu dust model and divide by stellar bolometric flux."""
    wavelengths = np.asarray(wavelength_um, dtype=float)
    fluxes = np.asarray(dust_flux_mjy, dtype=float)
    if np.any(wavelengths <= 0) or stellar_luminosity_solar <= 0 or distance_pc <= 0:
        raise ValueError("invalid SED integration inputs")
    frequency = C_LIGHT / (wavelengths * 1e-6)
    order = np.argsort(frequency)
    # Input flux densities are mJy, whereas JY_W_M2_HZ is defined per Jy.
    dust_bolometric_flux = np.trapezoid(fluxes[order] * 1e-3 * JY_W_M2_HZ, frequency[order])
    stellar_bolometric_flux = stellar_luminosity_solar * L_SUN_W / (4.0 * math.pi * (distance_pc * PC_M) ** 2)
    return float(dust_bolometric_flux / stellar_bolometric_flux)
