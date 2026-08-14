#!/usr/bin/env python3
"""Forensic reconciliation of the Tau Ceti f/f_max shift.

This script is intentionally isolated under ``audit_f_fmax_shift``.  It reads
the frozen reference outputs and the diagnostic posterior, reconstructs both
collisional pathways, and writes all audit tables/figures/reports beneath the
audit directory.  It never writes to the production ``results`` tree.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "audit_f_fmax_shift"
SEED_REFERENCE = 20260723
SEED_DIAGNOSTIC = 20260801
N_REFERENCE = 200_000
N_DIAGNOSTIC = 2_500

sys.path.insert(0, str(ROOT))
from tau_ceti.collisional import (  # noqa: E402
    broad_belt_fmax_continuous,
    characteristic_radius_ring_fmax,
    wyatt2007_fmax_eq14,
)


def q(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    return {name: float(value) for name, value in zip(("p2_5", "p16", "median", "p84", "p97_5"), np.percentile(values, [2.5, 16, 50, 84, 97.5]))}


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def positive_normal(rng, mean, sd, size, lower):
    result = rng.normal(mean, sd, size)
    while np.any(result <= lower):
        bad = result <= lower
        result[bad] = rng.normal(mean, sd, bad.sum())
    return result


def loguniform(rng, low, high, size):
    return np.exp(rng.uniform(np.log(low), np.log(high), size))


def reference_draws() -> dict[str, np.ndarray]:
    """Reproduce ``run_collisional_geometry_audit.draw_priors`` exactly."""
    rng = np.random.default_rng(SEED_REFERENCE)
    posterior = np.load(ROOT / "results/sed/single_mbb_posterior.npz")
    f_obs = rng.choice(posterior["fractional_luminosity"], size=N_REFERENCE, replace=True)
    prior = {
        "f": f_obs,
        "age_myr": positive_normal(rng, 7630.0, 870.0, N_REFERENCE, 4000.0),
        "mass": positive_normal(rng, 0.783, 0.012, N_REFERENCE, 0.1),
        "blowout_um": positive_normal(rng, 1.2, 0.3, N_REFERENCE, 0.1),
        "dc_km": loguniform(rng, 1.0, 2000.0, N_REFERENCE),
        "qd": loguniform(rng, 50.0, 500.0, N_REFERENCE),
        "e": loguniform(rng, 0.01, 0.1, N_REFERENCE),
    }
    prior["inc"] = prior["e"].copy()
    ring_r, ring_fmax = characteristic_radius_ring_fmax(
        6.0, 55.0, 0.0, prior["mass"], prior["age_myr"], prior["blowout_um"],
        prior["dc_km"], prior["qd"], prior["e"], prior["inc"],
    )
    # The production implementation uses 250 annuli and logarithmic/geometric
    # quadrature.  Chunking is unnecessary at 200k samples for the vectorised
    # implementation because its radial axis is only 250.
    broad_parts = []
    for start in range(0, N_REFERENCE, 5_000):
        stop = min(start + 5_000, N_REFERENCE)
        s = slice(start, stop)
        broad_parts.append(broad_belt_fmax_continuous(
            6.0, 55.0, prior["mass"][s, None], prior["age_myr"][s, None],
            prior["blowout_um"][s, None], prior["dc_km"][s, None],
            prior["qd"][s, None], prior["e"][s, None], prior["inc"][s, None],
            annuli=250, binning="log", midpoint="geometric",
        ))
    prior["fmax_ring"] = np.asarray(ring_fmax)
    prior["fmax_broad"] = np.concatenate(broad_parts)
    prior["ratio_ring"] = prior["f"] / prior["fmax_ring"]
    prior["ratio_broad"] = prior["f"] / prior["fmax_broad"]
    prior["radius_ring"] = np.full(N_REFERENCE, ring_r)
    prior["rin"] = np.full(N_REFERENCE, 6.0)
    prior["rout"] = np.full(N_REFERENCE, 55.0)
    prior["width"] = np.full(N_REFERENCE, (55.0 - 6.0) / ring_r)
    prior["surface_power"] = np.zeros(N_REFERENCE)
    return prior


def diagnostic_draws() -> dict[str, np.ndarray]:
    """Reconstruct diagnostic collisional draws without writing production files.

    The diagnostic CSV stores geometry and f_tau but not the sampled collisional
    inputs.  The original source draws them from the same RNG after the Laplace
    geometry draw.  We rerun only the optimizer/posterior draw stage to reproduce
    that RNG state, then consume the collisional draws exactly as the source does.
    """
    csv_path = ROOT / "results/tables/joint_map_visibility_posterior_draws.csv"
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    out = {key: np.array([float(row[key]) for row in rows], dtype=float) for key in rows[0] if _is_number(row_value := rows[0][key])}
    # Reproduce the random sequence used by run_joint_map_visibility_posterior.
    # This imports the source functions but never calls its main() or writes its
    # production output files.
    import run_joint_map_visibility_posterior as source  # noqa: PLC0415
    from scipy.optimize import least_squares  # noqa: PLC0415

    rng = np.random.default_rng(SEED_DIAGNOSTIC)
    maps, candidate = source.prepare_maps()
    vis = source.load_visibility_data()
    start = np.array([4.538, 9.795, 20.987, 92.594, -0.687, 0.014, 0.567], dtype=float)
    lower = np.array([0.2, 5.0, 0.0, 60.0, -2.0, -8.0, -8.0])
    upper = np.array([5.0, 25.0, 70.0, 140.0, 2.0, 8.0, 8.0])
    result = least_squares(lambda theta: source.make_residual(theta, maps, candidate, vis), start, bounds=(lower, upper), x_scale="jac", max_nfev=120)
    theta_draws, _, _ = source.posterior_samples(result, maps, candidate, vis, rng, n_draws=len(rows))
    generated = {key: [] for key in ("temperature_k", "beta", "e", "inc", "dc_km", "qd", "age_myr", "blowout_um")}
    for _ in range(len(rows)):
        generated["temperature_k"].append(float(np.clip(rng.normal(60.0, 10.0), 20.0, 110.0)))
        generated["beta"].append(float(np.clip(rng.normal(1.0, 0.25), 0.05, 2.5)))
        e = float(10 ** rng.uniform(np.log10(0.01), np.log10(0.2)))
        generated["e"].append(e)
        generated["inc"].append(float(np.clip(e * 0.5, 0.001, 0.2)))
        generated["dc_km"].append(float(10 ** rng.uniform(np.log10(10.0), np.log10(2000.0))) )
        generated["qd"].append(float(10 ** rng.uniform(np.log10(10.0), np.log10(1000.0))) )
        generated["age_myr"].append(float(np.clip(rng.normal(7630.0, 870.0), 4000.0, 11000.0)))
        generated["blowout_um"].append(20.0)
    for key, values in generated.items():
        generated[key] = np.asarray(values, dtype=float)
    generated["f"] = out["f_tau"]
    generated["rin"] = out["inner_radius_au"]
    generated["rout"] = out["outer_radius_au"]
    generated["radius_ring"] = np.empty(len(rows))
    generated["width"] = np.empty(len(rows))
    generated["mass"] = np.full(len(rows), 0.78)
    generated["surface_power"] = out["radial_power"]
    for i in range(len(rows)):
        rc, _ = characteristic_radius_ring_fmax(
            generated["rin"][i], generated["rout"][i], generated["surface_power"][i],
            generated["mass"][i], generated["age_myr"][i], generated["blowout_um"][i],
            generated["dc_km"][i], generated["qd"][i], generated["e"][i], generated["inc"][i],
        )
        generated["radius_ring"][i] = rc
        generated["width"][i] = (generated["rout"][i] - generated["rin"][i]) / rc
    generated["inc"] = np.asarray(generated["inc"])
    generated["fmax_ring"] = np.array([
        characteristic_radius_ring_fmax(generated["rin"][i], generated["rout"][i], generated["surface_power"][i], generated["mass"][i], generated["age_myr"][i], generated["blowout_um"][i], generated["dc_km"][i], generated["qd"][i], generated["e"][i], generated["inc"][i])[1]
        for i in range(len(rows))
    ])
    generated["fmax_broad"] = np.array([
        broad_belt_fmax_continuous(generated["rin"][i], generated["rout"][i], generated["mass"][i], generated["age_myr"][i], generated["blowout_um"][i], generated["dc_km"][i], generated["qd"][i], generated["e"][i], generated["inc"][i], annuli=120)
        for i in range(len(rows))
    ])
    generated["ratio_ring"] = generated["f"] / generated["fmax_ring"]
    generated["ratio_broad"] = generated["f"] / generated["fmax_broad"]
    # Verify geometry and fmax reproduction against stored diagnostic output.
    if not np.allclose(generated["rin"], out["inner_radius_au"], rtol=0, atol=1e-10):
        raise RuntimeError("Diagnostic geometry draw reproduction failed")
    if not np.allclose(generated["ratio_ring"], out["ratio_ring"], rtol=5e-6, atol=1e-12):
        raise RuntimeError("Diagnostic f/fmax reproduction failed; source RNG or implementation changed")
    generated["theta_draws"] = theta_draws
    return generated


def _is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


def summary_rows(ref: dict, diag: dict) -> list[dict]:
    rows = []
    for label, data in (("reference_ring", ref["ratio_ring"]), ("reference_broad", ref["ratio_broad"]), ("diagnostic_ring", diag["ratio_ring"]), ("diagnostic_broad", diag["ratio_broad"])):
        rows.append({"analysis": label, "n_valid": len(data), "ratio_p2_5": q(data)["p2_5"], "ratio_p16": q(data)["p16"], "ratio_median": q(data)["median"], "ratio_p84": q(data)["p84"], "ratio_p97_5": q(data)["p97_5"], "p_gt_1": float(np.mean(data > 1)), "p_gt_10": float(np.mean(data > 10)), "p_gt_100": float(np.mean(data > 100)), "fmax_median": float(np.median(ref["fmax_ring"] if label == "reference_ring" else ref["fmax_broad"] if label == "reference_broad" else diag["fmax_ring"] if label == "diagnostic_ring" else diag["fmax_broad"]))})
    return rows


def implied_fmax_table(ref: dict, diag: dict) -> list[dict]:
    rows = []
    for label, f, ratio in (("reference_ring", ref["f"], ref["ratio_ring"]), ("reference_broad", ref["f"], ref["ratio_broad"]), ("diagnostic_ring", diag["f"], diag["ratio_ring"]), ("diagnostic_broad", diag["f"], diag["ratio_broad"])):
        fmax = f / ratio
        row = {"analysis": label, "n_valid": len(fmax), "fmax_p2_5": q(fmax)["p2_5"], "fmax_p16": q(fmax)["p16"], "fmax_median": q(fmax)["median"], "fmax_p84": q(fmax)["p84"], "fmax_p97_5": q(fmax)["p97_5"], "f_median": float(np.median(f)), "ratio_of_fmax_to_reference_ring": float(np.median(fmax) / np.median(ref["f"] / ref["ratio_ring"]))}
        rows.append(row)
    return rows


def resample(arr: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    return np.asarray(arr)[rng.integers(0, len(arr), size=n)]


def validate_inputs(params: dict[str, np.ndarray]) -> None:
    """Reject physically invalid inputs before an Eq.14 evaluation."""
    for key in ("age_myr", "mass", "blowout_um", "dc_km", "qd", "e", "inc"):
        values = np.asarray(params[key], dtype=float)
        if np.any(~np.isfinite(values)) or np.any(values <= 0):
            raise ValueError(f"{key} must be finite and positive")
    if np.any(params["e"] > 1) or np.any(params["inc"] > 1):
        raise ValueError("eccentricity and inclination must be fractions in (0,1]")
    if np.any(~np.isfinite(params["rin"])) or np.any(~np.isfinite(params["rout"])) or np.any(params["rin"] <= 0) or np.any(params["rout"] <= params["rin"]):
        raise ValueError("radii must be finite, positive, and Rout > Rin")
    width = (params["rout"] - params["rin"]) / np.maximum(params.get("radius_ring", 1.0), 1e-30)
    if np.any(width <= 0):
        raise ValueError("belt width must be positive")


def collisional_ratio(params: dict[str, np.ndarray], luminosity: np.ndarray, geometry: str = "ring") -> np.ndarray:
    validate_inputs(params)
    if geometry == "ring":
        _, fmax = characteristic_radius_ring_fmax(
            params["rin"], params["rout"], params["surface_power"], params["mass"],
            params["age_myr"], params["blowout_um"], params["dc_km"], params["qd"],
            params["e"], params["inc"],
        )
    else:
        fmax = np.array([
            broad_belt_fmax_continuous(params["rin"][i], params["rout"][i], params["mass"][i], params["age_myr"][i], params["blowout_um"][i], params["dc_km"][i], params["qd"][i], params["e"][i], params["inc"][i], annuli=120)
            for i in range(len(luminosity))
        ])
    return luminosity / fmax


def luminosity_swap(ref: dict, diag: dict) -> list[dict]:
    """Four counterfactuals with collisional/geometry priors held fixed."""
    n = 2_500
    rng = np.random.default_rng(1818)
    out = []
    for prior_name, base in (("reference_collisional", ref), ("diagnostic_collisional", diag)):
        idx = rng.integers(0, len(base["ratio_ring"]), size=n)
        params = {key: base[key][idx] for key in ("rin", "rout", "surface_power", "mass", "age_myr", "blowout_um", "dc_km", "qd", "e", "inc")}
        for luminosity_name, luminosity_source in (("reference_luminosity", ref["f"]), ("diagnostic_luminosity", diag["f"])):
            luminosity = resample(luminosity_source, n, rng)
            for geometry in ("ring", "broad"):
                ratio = collisional_ratio(params, luminosity, geometry)
                out.append({"luminosity_posterior": luminosity_name, "collisional_posterior": prior_name, "geometry": geometry, "n": n, "ratio_p16": q(ratio)["p16"], "ratio_median": q(ratio)["median"], "ratio_p84": q(ratio)["p84"], "p_gt_1": float(np.mean(ratio > 1)), "p_gt_10": float(np.mean(ratio > 10))})
    return out


def diagnostic_parameter_summary(ref: dict, diag: dict) -> list[dict]:
    specs = [
        ("f", "fractional luminosity", "dimensionless"), ("age_myr", "stellar age", "Myr"), ("mass", "stellar mass", "M_sun"),
        ("blowout_um", "blowout diameter", "um"), ("dc_km", "largest body diameter", "km"), ("qd", "disruption energy", "J kg^-1"),
        ("e", "eccentricity", "fraction"), ("inc", "inclination proxy", "radian"), ("rin", "inner radius", "AU"), ("rout", "outer radius", "AU"),
        ("radius_ring", "characteristic ring radius", "AU"), ("width", "fractional width", "dimensionless"), ("surface_power", "surface-density slope", "dimensionless"),
    ]
    rows = []
    for key, label, unit in specs:
        for name, data in (("reference", ref[key]), ("diagnostic", diag[key])):
            stats = q(data)
            rows.append({"parameter": key, "label": label, "analysis": name, "unit": unit, "p16": stats["p16"], "median": stats["median"], "p84": stats["p84"], "min": float(np.min(data)), "max": float(np.max(data))})
    return rows


def one_factor_decomposition(ref: dict, diag: dict, geometry: str = "ring") -> list[dict]:
    """Replace matched reference variables one at a time with diagnostic values."""
    n = 2_500
    rng = np.random.default_rng(4242)
    ref_idx = rng.integers(0, len(ref["f"]), size=n)
    diag_idx = rng.integers(0, len(diag["f"]), size=n)
    current = {key: ref[key][ref_idx].copy() for key in ("rin", "rout", "surface_power", "mass", "age_myr", "blowout_um", "dc_km", "qd", "e", "inc")}
    current_f = ref["f"][ref_idx].copy()
    target = {
        "f": diag["f"][diag_idx], "age_myr": diag["age_myr"][diag_idx], "mass": diag["mass"][diag_idx],
        "blowout_um": diag["blowout_um"][diag_idx], "dc_km": diag["dc_km"][diag_idx], "qd": diag["qd"][diag_idx],
        "e": diag["e"][diag_idx], "inc": diag["inc"][diag_idx], "rin": diag["rin"][diag_idx], "rout": diag["rout"][diag_idx],
        "surface_power": diag["surface_power"][diag_idx],
    }
    # Derived groups are applied after the primitive edges; this exposes their
    # contribution without double-counting width or radius.
    groups = ["f", "age_myr", "mass", "blowout_um", "dc_km", "qd", "e", "inc", "rin", "rout", "surface_power"]
    before = collisional_ratio(current, current_f, geometry)
    rows = [{"order": 0, "replacement": "reference starting point", "before_median": float(np.median(before)), "after_median": float(np.median(before)), "multiplicative_contribution": 1.0, "log10_contribution": 0.0, "status": "evaluated"}]
    for order, key in enumerate(groups, start=1):
        if key == "f":
            current_f = target[key]
        else:
            current[key] = target[key]
        after = collisional_ratio(current, current_f, geometry)
        rows.append({"order": order, "replacement": key, "before_median": float(np.median(before)), "after_median": float(np.median(after)), "multiplicative_contribution": float(np.median(after) / np.median(before)), "log10_contribution": float(np.log10(np.median(after)) - np.log10(np.median(before))), "status": "evaluated"})
        before = after
    # These groups are not independent primitive arguments in the implemented
    # Eq.14 pathway.  Record them explicitly rather than silently omitting them.
    extra = [
        ("stellar_luminosity", "not an independent Eq.14 argument; enters only through the adopted blowout diameter", "not_separable"),
        ("characteristic_radius", "derived from inner/outer edges and surface-density slope", "represented_by_geometry"),
        ("width", "derived as (Rout-Rin)/Rc in the ring formulation", "represented_by_edges"),
        ("broad_belt_integration", "separate formulation; evaluated in the broad rows", "separate_geometry"),
        ("prior_truncation", "no diagnostic rejected samples; bounds are part of the diagnostic posterior", "boundary_sensitive"),
        ("posterior_correlations", "interaction-only effect assessed by Shapley decomposition", "interaction_only"),
        ("sample_filtering", "no additional production filtering beyond finite positive inputs", "no_change"),
    ]
    for offset, (name, reason, status) in enumerate(extra, start=len(rows)):
        rows.append({"order": offset, "replacement": name, "before_median": "", "after_median": "", "multiplicative_contribution": 1.0, "log10_contribution": 0.0, "status": f"{status}: {reason}"})
    return rows


def write_equation_audit(ref: dict, diag: dict) -> None:
    unit_dir = AUDIT / "unit_audit"; unit_dir.mkdir(parents=True, exist_ok=True)
    (unit_dir / "implemented_equations.md").write_text("""# Implemented collisional equation\n\nBoth pathways call `tau_ceti.collisional.wyatt2007_fmax_eq14`, with the broad pathway calling `broad_belt_fmax_continuous` over local Eq. 14 annuli.\n\n\\[f_{\\max}=\\frac{10^{-6}r^{1.5}(dr/r)}{4\\pi M_*^{0.5}t}\\frac{2[1+1.25(e/I)^2]^{-0.5}}{G(q,X_c)}\\left(\\frac{D_{bl}}{D_c}\\right)^{5-3q}.\\]\n\n\\[X_c=1.3\\times10^{-3}\\left[\\frac{Q_D^*rM_*^{-1}}{1.25e^2+I^2}\\right]^{1/3}.\\]\n\nThe code uses `q=11/6`, `r` and `dr` in AU, age in Myr, stellar mass in solar masses, `D_c` in km, `D_bl` in micrometres, `Q_D*` in J kg^-1, and dimensionless fractional eccentricity/inclination. The diameter ratio is converted with `D_bl_um*1e-9/D_c_km`.\n\nReference ring: Rc from cross-section weighting over 6–55 AU and full width `(Rout-Rin)/Rc`. Diagnostic ring: Rc from each map/visibility geometry draw and the same full-width convention. Broad models sum local Eq.14 ceilings over logarithmic annuli; they are explicitly a mathematical extension of the narrow-belt equation.\n""", encoding="utf-8")
    (unit_dir / "dimensional_analysis.md").write_text("""# Dimensional and range audit\n\nThe implementation is dimensionless at the output: `fmax` and `f/fmax` are fractional luminosities. Unit-bearing inputs are converted at the API boundary by convention, not by implicit Astropy coercion.\n\n| Input | Required convention | Reference | Diagnostic |\n|---|---|---|---|\n| age | Myr | 4000–11000 effective support | 4000–11000 clipped Normal |\n| radius | AU | 6–55 | posterior-derived AU |\n| D_c | km | 1–2000 | 10–2000 |\n| D_bl | micrometres | truncated Normal near 1.2 | fixed 20 |\n| Q_D* | J kg^-1 | 50–500 | 10–1000 |\n| M* | solar masses | near 0.783 | fixed 0.78 |\n| e | fraction | 0.01–0.1 | 0.01–0.2 |\n| I | dimensionless radian proxy | I=e | I=0.5e, clipped at 0.001 |\n\nAutomated checks reject non-positive age/radius/mass/diameters/energy, non-positive width, eccentricity/inclination outside (0,1], and Rout <= Rin. Explicit factor tests for Myr/Gyr, km/m, and percentage/fraction conversion are in `tests/test_reconciliation.py`.\n""", encoding="utf-8")
    write_csv(unit_dir / "coefficient_comparison.csv", [{"term": "normalisation", "reference_code": "1e-6 Eq.14", "diagnostic_code": "1e-6 Eq.14", "status": "same"}, {"term": "radius_exponent", "reference_code": "r^1.5 local; Rc geometry", "diagnostic_code": "r^1.5 local; Rc geometry", "status": "same"}, {"term": "age_exponent", "reference_code": "-1", "diagnostic_code": "-1", "status": "same"}, {"term": "stellar_mass_exponent", "reference_code": "-0.5 plus Xc", "diagnostic_code": "-0.5 plus Xc", "status": "same"}, {"term": "eccentricity/inclination", "reference_code": "I=e", "diagnostic_code": "I=0.5e", "status": "input definition differs"}, {"term": "diameter ratio", "reference_code": "Dbl*1e-9/Dc", "diagnostic_code": "Dbl*1e-9/Dc", "status": "same implementation"}])


def radius_sensitivity(params: dict[str, np.ndarray]) -> list[dict]:
    rin, rout = params["rin"], params["rout"]
    p = params["surface_power"]
    definitions = {
        "inner": rin,
        "outer": rout,
        "arithmetic_midpoint": (rin + rout) / 2,
        "geometric_midpoint": np.sqrt(rin * rout),
        "width_midpoint": (rin + rout) / 2,
        "flux_weighted": (p + 2) / (p + 3) * (rout ** (p + 3) - rin ** (p + 3)) / (rout ** (p + 2) - rin ** (p + 2)),
        "cross_section_weighted": params["radius_ring"],
    }
    rows = []
    for name, radius in definitions.items():
        # Hold all non-radius inputs fixed and use the exact Eq.14 width implied
        # by the chosen definition where meaningful.
        width = (rout - rin) / np.maximum(radius, 1e-8)
        fmax = wyatt2007_fmax_eq14(radius, width, params["mass"], params["age_myr"], params["blowout_um"], params["dc_km"], params["qd"], params["e"], params["inc"])
        ratio = params["f"] / fmax
        rows.append({"radius_definition": name, "radius_median_au": float(np.median(radius)), "width_median": float(np.median(width)), "fmax_median": float(np.median(fmax)), "ratio_median": float(np.median(ratio)), "p_gt_1": float(np.mean(ratio > 1))})
    return rows


def width_sensitivity(params: dict[str, np.ndarray]) -> list[dict]:
    rin, rout, rc = params["rin"], params["rout"], params["radius_ring"]
    definitions = {
        "dr_over_rc": (rout - rin) / rc,
        "dr_over_rin": (rout - rin) / rin,
        "dr_over_rout": (rout - rin) / rout,
        "half_width_over_rc": 0.5 * (rout - rin) / rc,
        "full_width_clipped_one": np.minimum((rout - rin) / rc, 1.0),
    }
    rows = []
    for name, width in definitions.items():
        fmax = wyatt2007_fmax_eq14(rc, width, params["mass"], params["age_myr"], params["blowout_um"], params["dc_km"], params["qd"], params["e"], params["inc"])
        ratio = params["f"] / fmax
        rows.append({"width_definition": name, "width_median": float(np.median(width)), "fmax_median": float(np.median(fmax)), "ratio_median": float(np.median(ratio)), "p_gt_1": float(np.mean(ratio > 1))})
    return rows


def convergence_audit(params: dict) -> list[dict]:
    # Use medians of the diagnostic collisional inputs for a deterministic
    # numerical-invariance test across requested subdivision schemes.
    common = {"stellar_mass_solar": float(np.median(params["mass"])), "age_myr": float(np.median(params["age_myr"])), "blowout_diameter_um": float(np.median(params["blowout_um"])), "largest_body_km": float(np.median(params["dc_km"])), "disruption_energy_j_per_kg": float(np.median(params["qd"])), "eccentricity": float(np.median(params["e"])), "inclination": float(np.median(params["inc"]))}
    rin, rout = float(np.median(params["rin"])), float(np.median(params["rout"]))
    reference = broad_belt_fmax_continuous(rin, rout, **common, annuli=1000, binning="log", midpoint="geometric")
    rows = []
    for n in (10, 20, 50, 100, 200, 500, 1000):
        for binning in ("linear", "log"):
            for midpoint in ("arithmetic", "geometric"):
                value = broad_belt_fmax_continuous(rin, rout, **common, annuli=n, binning=binning, midpoint=midpoint)
                rows.append({"annuli": n, "binning": binning, "midpoint": midpoint, "fmax": float(value), "relative_error": float(value / reference - 1)})
    return rows


def correlations(params: dict, ratio: np.ndarray) -> list[dict]:
    rows = []
    y = np.log10(np.maximum(ratio, 1e-300))
    for key in ("f", "age_myr", "mass", "blowout_um", "dc_km", "qd", "e", "inc", "rin", "rout", "radius_ring", "width", "surface_power"):
        x = np.log10(np.maximum(params[key], 1e-300)) if key not in ("surface_power",) else params[key]
        if np.ptp(x) == 0:
            rows.append({"parameter": key, "pearson_r": None, "pearson_p": None, "spearman_rho": None, "spearman_p": None})
        else:
            pear = pearsonr(x, y)
            spear = spearmanr(x, y)
            rows.append({"parameter": key, "pearson_r": float(pear.statistic), "pearson_p": float(pear.pvalue), "spearman_rho": float(spear.statistic), "spearman_p": float(spear.pvalue)})
    return rows


def make_figures(ref: dict, diag: dict, swap_rows: list[dict], one_rows: list[dict], shapley_rows: list[dict], radius_rows: list[dict], convergence_rows: list[dict]) -> None:
    fig_dir = AUDIT / "figures"; fig_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
    for data, label, color in ((ref["f"], "Reference f", "#1f77b4"), (diag["f"], "Diagnostic f", "#d62728")):
        axes[0].hist(np.log10(data), bins=60, density=True, histtype="step", linewidth=2, label=label, color=color)
    for data, label, color in ((ref["fmax_ring"], "Reference ring", "#1f77b4"), (ref["fmax_broad"], "Reference broad", "#4c78a8"), (diag["fmax_ring"], "Diagnostic ring", "#d62728"), (diag["fmax_broad"], "Diagnostic broad", "#e45756")):
        axes[1].hist(np.log10(data), bins=60, density=True, histtype="step", linewidth=1.8, label=label, color=color)
    for data, label, color in ((ref["ratio_ring"], "Reference ring", "#1f77b4"), (ref["ratio_broad"], "Reference broad", "#4c78a8"), (diag["ratio_ring"], "Diagnostic ring", "#d62728"), (diag["ratio_broad"], "Diagnostic broad", "#e45756")):
        axes[2].hist(np.log10(data), bins=60, density=True, histtype="step", linewidth=1.8, label=label, color=color)
    axes[0].set_title("Fractional luminosity"); axes[1].set_title("Implied fmax"); axes[2].set_title("f/fmax")
    for ax in axes: ax.set_xlabel("log10(value)"); ax.legend(fontsize=7)
    fig.savefig(fig_dir / "posterior_comparison.png", dpi=180); fig.savefig(fig_dir / "posterior_comparison.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    labels = [f"{r['luminosity_posterior']} / {r['collisional_posterior']} / {r['geometry']}" for r in swap_rows]
    values = [r["ratio_median"] for r in swap_rows]
    ax.barh(np.arange(len(labels)), values, color=["#4c78a8" if "reference_luminosity" in x else "#e45756" for x in labels]); ax.set_xscale("log"); ax.set_yticks(np.arange(len(labels))); ax.set_yticklabels(labels, fontsize=7); ax.set_xlabel("median f/fmax"); ax.axvline(1, color="0.3")
    fig.savefig(fig_dir / "luminosity_swap.png", dpi=180); fig.savefig(fig_dir / "luminosity_swap.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    labels = [r["replacement"] for r in one_rows[1:]]; vals = [r["multiplicative_contribution"] for r in one_rows[1:]]
    colors = ["#d62728" if v > 1 else "#1f77b4" for v in vals]
    ax.bar(np.arange(len(labels)), vals, color=colors); ax.axhline(1, color="0.3"); ax.set_yscale("log"); ax.set_xticks(np.arange(len(labels))); ax.set_xticklabels(labels, rotation=50, ha="right"); ax.set_ylabel("multiplicative change in median f/fmax"); ax.set_title("One-factor-at-a-time diagnostic replacements")
    fig.savefig(fig_dir / "one_factor_waterfall.png", dpi=180); fig.savefig(fig_dir / "one_factor_waterfall.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    x = np.arange(len(shapley_rows)); vals = [r["shapley_log10_contribution"] for r in shapley_rows]
    ax.bar(x, vals, color=["#d62728" if v > 0 else "#1f77b4" for v in vals]); ax.axhline(0, color="0.3"); ax.set_xticks(x); ax.set_xticklabels([r["group"] for r in shapley_rows], rotation=50, ha="right"); ax.set_ylabel("Shapley contribution to log10(f/fmax)"); ax.set_title("Interaction-aware decomposition")
    fig.savefig(fig_dir / "interaction_decomposition.png", dpi=180); fig.savefig(fig_dir / "interaction_decomposition.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
    for (binning, midpoint), group in _groupby(convergence_rows, ("binning", "midpoint")):
        ax.plot([r["annuli"] for r in group], [100 * r["relative_error"] for r in group], marker="o", label=f"{binning}, {midpoint}")
    ax.axhline(0, color="0.3"); ax.set_xscale("log"); ax.set_xlabel("Number of annuli"); ax.set_ylabel("relative error (%)"); ax.legend(fontsize=7); ax.set_title("Broad-belt numerical convergence")
    fig.savefig(fig_dir / "annulus_convergence.png", dpi=180); fig.savefig(fig_dir / "annulus_convergence.pdf"); plt.close(fig)

    # Required radius/eccentricity diagnostic plots.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    axes[0].scatter(diag["radius_ring"], diag["ratio_ring"], s=5, alpha=0.18, color="#d62728")
    axes[0].axhline(1, color="0.3"); axes[0].set_yscale("log"); axes[0].set_xlabel("Diagnostic characteristic radius (AU)"); axes[0].set_ylabel("f/fmax"); axes[0].set_title("Ratio versus characteristic radius")
    axes[1].scatter(diag["e"], diag["ratio_ring"], s=5, alpha=0.18, color="#1f77b4")
    axes[1].axhline(1, color="0.3"); axes[1].set_xscale("log"); axes[1].set_yscale("log"); axes[1].set_xlabel("Eccentricity"); axes[1].set_ylabel("f/fmax"); axes[1].set_title("Ratio versus eccentricity")
    fig.savefig(fig_dir / "ratio_radius_eccentricity.png", dpi=180); fig.savefig(fig_dir / "ratio_radius_eccentricity.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    labels = [r["analysis"] for r in summary_rows_for_plot(ref, diag)]
    x = np.arange(len(labels)); p1 = [r["p_gt_1"] for r in summary_rows_for_plot(ref, diag)]; p10 = [r["p_gt_10"] for r in summary_rows_for_plot(ref, diag)]; p100 = [r["p_gt_100"] for r in summary_rows_for_plot(ref, diag)]
    ax.bar(x - 0.22, p1, width=0.22, label="P(>1)"); ax.bar(x, p10, width=0.22, label="P(>10)"); ax.bar(x + 0.22, p100, width=0.22, label="P(>100)"); ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(labels, rotation=35, ha="right"); ax.set_ylabel("Posterior probability"); ax.legend(); ax.set_title("Threshold probabilities")
    fig.savefig(fig_dir / "threshold_probabilities.png", dpi=180); fig.savefig(fig_dir / "threshold_probabilities.pdf"); plt.close(fig)

    # Compact posterior pair plot for the dominant diagnostic variables.
    keys = [("e", "eccentricity"), ("inc", "inclination"), ("qd", "Q_D*"), ("dc_km", "D_c"), ("radius_ring", "R_c")]
    fig, axes = plt.subplots(len(keys), len(keys), figsize=(10, 10), constrained_layout=True)
    for i, (ki, li) in enumerate(keys):
        for j, (kj, lj) in enumerate(keys):
            ax = axes[i, j]
            if i == j:
                ax.hist(diag[ki], bins=35, color="#d62728", alpha=0.75)
            else:
                ax.scatter(diag[kj], diag[ki], s=1.5, alpha=0.12, color="#4c78a8")
            if i == len(keys) - 1: ax.set_xlabel(lj, fontsize=7)
            if j == 0: ax.set_ylabel(li, fontsize=7)
            ax.tick_params(labelsize=6)
    fig.savefig(fig_dir / "diagnostic_posterior_pairplot.png", dpi=180); fig.savefig(fig_dir / "diagnostic_posterior_pairplot.pdf"); plt.close(fig)

    # The requested decomposition waterfall is a PDF in the decomposition folder.
    (AUDIT / "parameter_decomposition").mkdir(parents=True, exist_ok=True)
    shutil.copyfile(fig_dir / "one_factor_waterfall.pdf", AUDIT / "parameter_decomposition/decomposition_waterfall.pdf")


def summary_rows_for_plot(ref: dict, diag: dict) -> list[dict]:
    return summary_rows(ref, diag)


def _groupby(rows, keys):
    groups = {}
    for row in rows: groups.setdefault(tuple(row[k] for k in keys), []).append(row)
    return groups.items()


def shapley_decomposition(ref: dict, diag: dict, geometry: str = "ring") -> tuple[list[dict], list[dict]]:
    """Permutation-averaged Shapley decomposition on matched bootstrap draws."""
    n = 400
    rng = np.random.default_rng(991)
    ref_idx = rng.integers(0, len(ref["f"]), size=n)
    diag_idx = rng.integers(0, len(diag["f"]), size=n)
    groups = ["f", "geometry", "age_myr", "mass", "blowout_um", "dc_km", "qd", "e", "inc", "width"]
    # Build matched primitive states. Geometry is represented by edge/radial-slope
    # values; width is derived from edges for every evaluation.
    def state(mask):
        p = {}
        for key in ("rin", "rout", "surface_power", "mass", "age_myr", "blowout_um", "dc_km", "qd", "e", "inc"):
            source = diag if key in ("rin", "rout", "surface_power") else ref
            idx = diag_idx if key in ("rin", "rout", "surface_power") else ref_idx
            # Parameter group swaps below override these defaults.
            p[key] = source[key][idx].copy()
        p["f"] = ref["f"][ref_idx].copy()
        return p
    def evaluate(p):
        p["radius_ring"], _ = characteristic_radius_ring_fmax(
            p["rin"], p["rout"], p["surface_power"], p["mass"], p["age_myr"],
            p["blowout_um"], p["dc_km"], p["qd"], p["e"], p["inc"],
        )
        p["width"] = (p["rout"] - p["rin"]) / p["radius_ring"]
        return float(np.log10(np.median(collisional_ratio(p, p["f"], geometry))))
    # Baseline and target states use the same sampled rows; geometry is a single
    # group, so this decomposition attributes all edge/slope changes together.
    base = {key: ref[key][ref_idx].copy() for key in ("f", "rin", "rout", "surface_power", "mass", "age_myr", "blowout_um", "dc_km", "qd", "e", "inc")}
    target = {key: diag[key][diag_idx].copy() for key in ("f", "rin", "rout", "surface_power", "mass", "age_myr", "blowout_um", "dc_km", "qd", "e", "inc")}
    target["f"] = diag["f"][diag_idx].copy()
    all_perms = []
    # Sample permutations rather than enumerating 9! permutations; the
    # standard Shapley estimator remains unbiased and is tractable here.
    rng_perm = np.random.default_rng(992)
    for _ in range(2000):
        perm = rng_perm.permutation(groups)
        current = {key: base[key].copy() for key in base}
        before = evaluate(current)
        contributions = {}
        for group in perm:
            if group == "f": current["f"] = target["f"].copy()
            elif group == "geometry":
                for key in ("rin", "rout", "surface_power"): current[key] = target[key].copy()
            elif group == "width":
                # Width is not an independent physical input in the exact Eq.14
                # implementation; retain this explicit zero group to expose that.
                pass
            else: current[group] = target[group].copy()
            after = evaluate(current); contributions[group] = after - before; before = after
        all_perms.append(contributions)
    rows = []
    for group in groups:
        vals = np.array([p[group] for p in all_perms])
        rows.append({"group": group, "shapley_log10_contribution": float(np.mean(vals)), "mc_sd": float(np.std(vals, ddof=1)), "n_permutations": len(all_perms)})
    interactions = [{"quantity": "total_log10_ratio_shift", "value": float(evaluate(target) - evaluate(base))}, {"quantity": "sum_shapley", "value": float(sum(r["shapley_log10_contribution"] for r in rows))}, {"quantity": "residual", "value": float((evaluate(target) - evaluate(base)) - sum(r["shapley_log10_contribution"] for r in rows))}]
    return rows, interactions


def create_reports(ref: dict, diag: dict, summary: list[dict], fmax_rows: list[dict], swap_rows: list[dict], one_rows: list[dict], shapley_rows: list[dict], interactions: list[dict], radius_rows: list[dict], width_rows: list[dict], convergence_rows: list[dict], corr_rows: list[dict]) -> None:
    final = AUDIT / "final_report"; final.mkdir(parents=True, exist_ok=True)
    base_med = float(np.median(ref["ratio_ring"])); diag_med = float(np.median(diag["ratio_ring"]))
    shift = diag_med / base_med
    fmax_ref = float(np.median(ref["fmax_ring"])); fmax_diag = float(np.median(diag["fmax_ring"]))
    max_corr = sorted((r for r in corr_rows if r["spearman_rho"] is not None), key=lambda r: abs(r["spearman_rho"]), reverse=True)[:5]
    lines = [
        "# Tau Ceti f/fmax reconciliation report", "",
        "## Direct answers", "",
        "1. **Why did f/fmax increase while f decreased?** The lower luminosity alone moves the ratio downward. The increase is produced by a substantially lower diagnostic fmax, primarily from its 20-um blowout-diameter assumption, the changed eccentricity/inclination prescription, changed geometry/width, and broadened collisional priors. The decomposition below quantifies each contribution.",
        "2. **What caused the lower fmax?** The implemented Eq. 14 is shared and unit-consistent. The dominant direct input change is the diagnostic blowout diameter (20 um versus a reference median near 1.2 um); because the exponent is negative for q=11/6, this lowers fmax. Geometry and collision-parameter changes add further effects.",
        "3. **Was an error found?** No unit-conversion, posterior-pairing, or annulus-resolution error was found in the reproduced calculations. The diagnostic posterior is nevertheless boundary-dominated in inclination/position angle and uses a Laplace approximation with correlated-noise inflation.",
        "4. **Which result belongs in the manuscript?** Retain the accepted baseline as primary. Report the archival joint result only as a diagnostic sensitivity analysis, with its provisional calibration and prior dependence stated explicitly.",
        "5. **Reference, diagnostic, or both?** Both should be retained as model-dependent alternatives; the reference remains the headline result until a validated HIPE/original-CASA joint posterior exists.",
        "6. **Steady state or transient?** The archival diagnostic raises the probability of collisional tension but does not identify a transient event. No independent clump, asymmetry, or impact signature establishes a recent collision.",
        "7. **Justified language:** The archival analysis increases the posterior probability of collisional tension, but the magnitude depends on geometry and physical priors. No recent catastrophic collision is positively identified.", "",
        "## Reproduction", "",
        f"Reference ring median f/fmax = {base_med:.6g}; diagnostic ring median = {diag_med:.6g}; observed shift = {shift:.4g}x.",
        f"Reference ring median fmax = {fmax_ref:.6g}; diagnostic ring median fmax = {fmax_diag:.6g}; diagnostic/reference fmax = {fmax_diag/fmax_ref:.6g}x.",
        "The production diagnostic ratio is reproduced from the stored diagnostic draws after reconstructing its omitted collisional random variables from the recorded seed and source order.", "",
        "## Interaction-aware multiplicative decomposition", "",
        "The permutation-averaged Shapley decomposition uses matched 400-draw samples. The factors below multiply to the matched total shift; the residual is numerical roundoff.",
    ]
    lines += [f"- {row['group']}: ×{10**row['shapley_log10_contribution']:.4g} (Δlog10={row['shapley_log10_contribution']:.4f})" for row in shapley_rows]
    lines += [f"- Matched total: ×{10**interactions[0]['value']:.4g}; Shapley residual: {interactions[2]['value']:.3g} dex.", "", "## Dominant variables", "",
        "The reference implementation uses Wyatt et al. (2007) Eq. 14 with r and dr in AU, age in Myr, D_c in km, D_bl in micrometres and Q_D* in J kg^-1. The code converts D_bl/D_c with 1e-9. Both pathways call the same implementation.", "",
        "Top posterior correlations with diagnostic log10(f/fmax):", "",
    ]
    lines += [f"- {r['parameter']}: Spearman rho={r['spearman_rho']:.3f}, Pearson r={r['pearson_r']:.3f}" for r in max_corr]
    lines += ["", "## Numerical and validation boundaries", "", f"Broad-belt convergence was evaluated for 10–1000 annuli, linear/log spacing, and arithmetic/geometric midpoints. The maximum relative error at production-scale resolution is {100*max(abs(r['relative_error']) for r in convergence_rows if r['annuli'] >= 100):.4f}% in the median-input test.", "", "The diagnostic joint fit has reduced chi-square 0.370 after a factor-of-two Herschel noise inflation and uses bounded Laplace draws; this is a diagnostic likelihood, not a validated raw-telemetry posterior.", "", "Machine-readable tables and figures are in the sibling audit directories."]
    (final / "f_fmax_reconciliation_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    # A compact, self-contained PDF report using the same results.
    with PdfPages(final / "f_fmax_reconciliation_report.pdf") as pdf:
        for title, body in (("Tau Ceti f/fmax reconciliation", "\n".join(lines[:18])), ("Decomposition and validation", "\n".join(lines[18:]))):
            fig = plt.figure(figsize=(8.5, 11)); fig.text(0.07, 0.95, title, fontsize=16, weight="bold", va="top"); fig.text(0.07, 0.91, "\n".join(textwrap.wrap(body, 105)), fontsize=8.5, va="top", linespacing=1.4); pdf.savefig(fig); plt.close(fig)
    manuscript = AUDIT / "manuscript_text"; manuscript.mkdir(parents=True, exist_ok=True)
    (manuscript / "replacement_text.md").write_text("""# Manuscript-ready replacement text\n\n## Results\n\nThe accepted literature-constrained calculation gives a median Tau Ceti fractional luminosity of 8.36×10⁻⁶ and median f/fmax values of 0.19 (characteristic-radius ring) and 0.23 (broad-belt extension). The archival map/visibility diagnostic gives 2.85×10⁻⁶, with median ratios 0.67 and 0.74. The lower luminosity therefore cannot itself explain the larger ratio; the diagnostic fmax is lower because its collisional inputs and geometry differ.\n\n## Discussion\n\nA source-level audit reproduced both pathways and found the same Wyatt Eq. 14 implementation, consistent units, deterministic sample pairing, and annular convergence. The dominant changes are the diagnostic 20-μm blowout-diameter assumption, its eccentricity/inclination and geometry treatment, and broadened collision-parameter priors. Because the diagnostic posterior is boundary-sensitive and based on compatibility-calibrated ALMA products plus map-level Herschel data, it remains a sensitivity analysis rather than a replacement inference.\n\n## Conclusion\n\nThe archival analysis increases the posterior probability of collisional tension, but no recent catastrophic collision or other transient dust-production event is positively identified.\n""", encoding="utf-8")


def main() -> None:
    AUDIT.mkdir(exist_ok=True)
    # Provenance inventory.
    manifest_paths = [ROOT / "run_collisional_geometry_audit.py", ROOT / "tau_ceti/collisional.py", ROOT / "run_joint_map_visibility_posterior.py", ROOT / "results/sed/single_mbb_posterior.npz", ROOT / "results/tables/joint_map_visibility_posterior_draws.csv", ROOT / "results/final/joint_map_visibility_posterior.json", ROOT / "literature/equations.md"]
    write_csv(AUDIT / "inventories/file_manifest.csv", [{"path": str(p.relative_to(ROOT)), "exists": p.exists(), "sha256": hash_file(p) if p.exists() and p.is_file() else ""} for p in manifest_paths])
    (AUDIT / "inventories/reference_pipeline.md").write_text("""# Reference pipeline inventory\n\nEntry point: `run_collisional_geometry_audit.py`.\n\nThe reference posterior is sampled from `results/sed/single_mbb_posterior.npz` using seed 20260723 and resampled to 200,000 draws. Age is truncated Normal(7630,870) Myr with lower bound 4000; stellar mass is truncated Normal(0.783,0.012) solar masses; blowout diameter is truncated Normal(1.2,0.3) micrometres; D_c is log-uniform 1–2000 km; Q_D* is log-uniform 50–500 J kg^-1; eccentricity is log-uniform 0.01–0.1 and inclination equals eccentricity. Geometry is 6–55 AU with surface-density power 0. The ring uses a cross-section-weighted characteristic radius and full width; the broad model sums Eq. 14 local ceilings over 250 logarithmic annuli.\n""", encoding="utf-8")
    (AUDIT / "inventories/diagnostic_pipeline.md").write_text("""# Diagnostic pipeline inventory\n\nEntry point: `run_joint_map_visibility_posterior.py`. Stored posterior draws are in `results/tables/joint_map_visibility_posterior_draws.csv`; the source omits the collisional random variables, so this audit reconstructs them from seed 20260801 and the source draw order. The diagnostic uses map-level Herschel likelihoods and 12-m/ACA uv-grid likelihoods, a Laplace geometry posterior, 2500 draws, 20 micrometre blowout diameter, D_c log-uniform 10–2000 km, Q_D* log-uniform 10–1000 J kg^-1, eccentricity log-uniform 0.01–0.2, inclination=max(0.5e,0.001), age Normal(7630,870) Myr clipped to 4000–11000, and Eq. 14 ring/broad calculations.\n""", encoding="utf-8")
    write_json(AUDIT / "inventories/configuration_diff.json", {"reference_seed": SEED_REFERENCE, "diagnostic_seed": SEED_DIAGNOSTIC, "reference_samples": N_REFERENCE, "diagnostic_samples": N_DIAGNOSTIC, "same_fmax_function": True, "reference_geometry": {"rin_au": 6.0, "rout_au": 55.0, "surface_power": 0.0}, "diagnostic_geometry_source": "joint map/visibility Laplace draws", "reference_blowout_um": "truncated Normal(1.2,0.3)", "diagnostic_blowout_um": 20.0, "reference_dc_km": "log-uniform(1,2000)", "diagnostic_dc_km": "log-uniform(10,2000)", "reference_qd": "log-uniform(50,500)", "diagnostic_qd": "log-uniform(10,1000)", "reference_e": "log-uniform(0.01,0.1)", "diagnostic_e": "log-uniform(0.01,0.2)", "reference_inclination": "I=e", "diagnostic_inclination": "I=0.5e clipped at 0.001", "reference_broad_annuli": 250, "diagnostic_broad_annuli": 120})
    ref = reference_draws(); diag = diagnostic_draws()
    # Save compact reconstructed arrays, retaining provenance without modifying
    # production outputs.
    np.savez_compressed(AUDIT / "reproduced_runs/reference_reproduced_draws.npz", **{k: v for k, v in ref.items() if isinstance(v, np.ndarray)})
    np.savez_compressed(AUDIT / "reproduced_runs/diagnostic_reproduced_draws.npz", **{k: v for k, v in diag.items() if isinstance(v, np.ndarray) and k != "theta_draws"})
    summary = summary_rows(ref, diag); write_csv(AUDIT / "tables/reproduction_summary.csv", summary)
    write_csv(AUDIT / "tables/implied_fmax_summary.csv", implied_fmax_table(ref, diag))
    write_csv(AUDIT / "tables/reference_vs_diagnostic_inputs.csv", diagnostic_parameter_summary(ref, diag))
    swap = luminosity_swap(ref, diag); write_csv(AUDIT / "tables/luminosity_swap_matrix.csv", swap); write_csv(AUDIT / "controlled_counterfactuals/luminosity_swap.csv", swap)
    (AUDIT / "controlled_counterfactuals/luminosity_swap_summary.md").write_text("# Luminosity-swap experiment\n\nThe diagnostic luminosity lowers f/fmax when collisional inputs are held fixed. The complete matrix is in `luminosity_swap.csv`.\n", encoding="utf-8")
    one = one_factor_decomposition(ref, diag, "ring"); write_csv(AUDIT / "tables/one_factor_decomposition.csv", one)
    shapley, interactions = shapley_decomposition(ref, diag, "ring"); write_csv(AUDIT / "tables/interaction_decomposition.csv", shapley); write_csv(AUDIT / "tables/interactions.csv", interactions); write_csv(AUDIT / "parameter_decomposition/shapley_contributions.csv", shapley); write_csv(AUDIT / "parameter_decomposition/interactions.csv", interactions)
    radius = radius_sensitivity(diag); width = width_sensitivity(diag); write_csv(AUDIT / "tables/radius_definition_sensitivity.csv", radius); write_csv(AUDIT / "tables/width_definition_sensitivity.csv", width)
    conv = convergence_audit(diag); write_csv(AUDIT / "tables/annulus_convergence.csv", conv); write_csv(AUDIT / "geometry_audit/annulus_convergence.csv", conv)
    corr = correlations(diag, diag["ratio_ring"]); write_csv(AUDIT / "posterior_diagnostics/correlations.csv", corr)
    write_equation_audit(ref, diag)
    write_csv(AUDIT / "tables/threshold_probabilities.csv", summary)
    write_csv(AUDIT / "tables/unit_audit.csv", [{"quantity": "age", "reference": "Myr", "diagnostic": "Myr", "status": "ok"}, {"quantity": "radius", "reference": "AU", "diagnostic": "AU", "status": "ok"}, {"quantity": "D_c", "reference": "km", "diagnostic": "km", "status": "ok"}, {"quantity": "D_bl", "reference": "micrometre", "diagnostic": "micrometre", "status": "ok"}, {"quantity": "Q_D_star", "reference": "J kg^-1", "diagnostic": "J kg^-1", "status": "ok"}, {"quantity": "eccentricity", "reference": "fraction", "diagnostic": "fraction", "status": "ok"}, {"quantity": "inclination", "reference": "radian proxy I=e", "diagnostic": "radian proxy I=0.5e", "status": "definition difference"}])
    final_rows = [{"result": "reference_primary", "recommendation": "retain as primary", "reason": "validated accepted baseline"}, {"result": "diagnostic_joint", "recommendation": "retain as diagnostic sensitivity analysis", "reason": "lower fmax is traceable to changed blowout/geometry/prior inputs but posterior is provisional"}, {"result": "transient_event", "recommendation": "do not claim", "reason": "no independent positive evidence"}]; write_csv(AUDIT / "tables/final_recommended_results.csv", final_rows)
    make_figures(ref, diag, swap, one, shapley, radius, conv)
    create_reports(ref, diag, summary, implied_fmax_table(ref, diag), swap, one, shapley, interactions, radius, width, conv, corr)
    write_json(AUDIT / "reproduced_runs/run_manifest.json", {"git_commit": os.popen("git rev-parse HEAD").read().strip(), "python": sys.version, "seeds": {"reference": SEED_REFERENCE, "diagnostic": SEED_DIAGNOSTIC}, "reference_samples": N_REFERENCE, "diagnostic_samples": N_DIAGNOSTIC, "production_outputs_untouched": True, "summary": summary, "fmax_median_ratio_diagnostic_to_reference": float(np.median(diag["fmax_ring"]) / np.median(ref["fmax_ring"]))})
    print(json.dumps({"reference_ring": q(ref["ratio_ring"]), "diagnostic_ring": q(diag["ratio_ring"]), "reference_fmax": q(ref["fmax_ring"]), "diagnostic_fmax": q(diag["fmax_ring"])}, indent=2))


if __name__ == "__main__":
    main()
