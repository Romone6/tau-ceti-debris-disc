"""Restricted, reproducible censored fractional-luminosity population utilities."""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.special import log_ndtr


def censored_normal_negative_log_likelihood(theta, x, detected, y, y_sigma, limits):
    """Normal latent log-f population with star-specific, temperature-grid limits."""
    alpha, beta, log_sigma = theta
    sigma = np.exp(log_sigma)
    mu = alpha + beta * x
    value = 0.0
    for index in range(len(mu)):
        if detected[index]:
            total = np.sqrt(sigma**2 + y_sigma[index]**2)
            value += 0.5 * ((y[index] - mu[index]) / total) ** 2 + np.log(total)
        else:
            # Temperature-marginalised censoring probability p(y < limit(T)).
            log_terms = log_ndtr((limits[index] - mu[index]) / sigma)
            value -= float(np.logaddexp.reduce(log_terms) - np.log(len(log_terms)))
    return float(value)


def fit_censored_age_model(x, detected, y, y_sigma, limits):
    """Maximum-likelihood censored age relation; a compact baseline, not MCMC."""
    result = minimize(
        censored_normal_negative_log_likelihood,
        x0=np.array([-5.5, -0.3, np.log(0.8)]),
        args=(np.asarray(x), np.asarray(detected), np.asarray(y), np.asarray(y_sigma), np.asarray(limits)),
        method="L-BFGS-B",
        bounds=[(-10, -2), (-5, 5), (np.log(0.05), np.log(3.0))],
    )
    if not result.success:
        raise RuntimeError(f"censored fit failed: {result.message}")
    return result.x


def censored_regression_negative_log_likelihood(theta, design, detected, y, y_sigma, limits, groups=None):
    """Censored Gaussian regression with optional survey-specific scatter.

    ``design`` includes the intercept.  When ``groups`` is supplied, the final
    two parameters are log intrinsic scatters for groups 0 and 1; otherwise the
    final parameter is one shared log scatter.  This is a likelihood fit used
    for sensitivity analysis, not a Bayesian posterior sampler.
    """
    design = np.asarray(design, dtype=float)
    detected = np.asarray(detected, dtype=bool)
    y = np.asarray(y, dtype=float)
    y_sigma = np.asarray(y_sigma, dtype=float)
    limits = np.asarray(limits, dtype=float)
    coefficient_count = design.shape[1]
    mu = design @ np.asarray(theta[:coefficient_count])
    if groups is None:
        sigma = np.full(len(mu), np.exp(theta[coefficient_count]))
    else:
        groups = np.asarray(groups, dtype=int)
        sigma = np.where(groups == 0, np.exp(theta[coefficient_count]), np.exp(theta[coefficient_count + 1]))
    value = 0.0
    for index in range(len(mu)):
        if detected[index]:
            total = np.sqrt(sigma[index] ** 2 + y_sigma[index] ** 2)
            value += 0.5 * ((y[index] - mu[index]) / total) ** 2 + np.log(total)
        else:
            log_terms = log_ndtr((limits[index] - mu[index]) / sigma[index])
            value -= float(np.logaddexp.reduce(log_terms) - np.log(len(log_terms)))
    return float(value)


def fit_censored_regression(design, detected, y, y_sigma, limits, groups=None):
    """Fit a reproducible censored regression and return its MLE parameters."""
    design = np.asarray(design, dtype=float)
    coefficient_count = design.shape[1]
    initial = np.r_[np.full(coefficient_count, 0.0), np.log(0.8)]
    initial[0] = -5.5
    if groups is not None:
        initial = np.r_[initial, np.log(0.8)]
    bounds = [(-10, -2)] + [(-5, 5)] * (coefficient_count - 1) + [(np.log(0.05), np.log(3.0))] * (2 if groups is not None else 1)
    result = minimize(
        censored_regression_negative_log_likelihood,
        x0=initial,
        args=(design, detected, y, y_sigma, limits, groups),
        method="L-BFGS-B",
        bounds=bounds,
    )
    if not result.success:
        raise RuntimeError(f"censored regression failed: {result.message}")
    return result.x
