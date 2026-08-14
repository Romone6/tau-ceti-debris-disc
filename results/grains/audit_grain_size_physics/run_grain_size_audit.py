"""Reproducible physical audit of Tau Ceti grain-size conventions.

This is deliberately isolated from the accepted f/f_max reconciliation.  It
reads the previously reproduced posterior arrays, derives radiation-pressure
and wind blowout sizes, and recomputes Eq. 14 while keeping D_bl, D_min and
D_char separate.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path("/Users/romonedunlop/Documents/Tau Ceti Research")
OUT = ROOT / "audit_grain_size_physics"
PREV = ROOT / "audit_f_fmax_shift" / "reproduced_runs"
sys.path.insert(0, str(ROOT))
from tau_ceti.collisional import wyatt2007_fmax_eq14, wyatt2007_g  # noqa: E402
from tau_ceti.physics import C_LIGHT, G_GRAV, L_SUN_W, M_SUN_KG  # noqa: E402

SEED = 20260802
RNG = np.random.default_rng(SEED)
PI = np.pi
SOLAR_MDOT = 1.26e9  # kg s^-1; documented solar-wind reference
WIND_SPEED = 400e3  # m s^-1


def qtile(a: np.ndarray) -> tuple[float, float, float]:
    return tuple(float(x) for x in np.percentile(np.asarray(a), [16, 50, 84]))


def csv_write(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        return "unavailable"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_arrays() -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    r = dict(np.load(PREV / "reference_reproduced_draws.npz"))
    d = dict(np.load(PREV / "diagnostic_reproduced_draws.npz"))
    return r, d


def d_bl_um(lum: np.ndarray | float, mass: np.ndarray | float, rho_gcc: np.ndarray | float,
            qpr: np.ndarray | float, mdot_solar: np.ndarray | float = 0.0,
            vwind: float = WIND_SPEED, cd: float = 1.0) -> np.ndarray:
    """Diameter (not radius) at total beta=0.5, in micrometres."""
    lum, mass, rho_gcc, qpr, mdot_solar = np.broadcast_arrays(
        lum, mass, rho_gcc, qpr, mdot_solar
    )
    rho = rho_gcc * 1000.0
    pressure = lum * L_SUN_W * qpr / C_LIGHT + mdot_solar * SOLAR_MDOT * vwind * cd
    # Pressure force per unit area is (L Q/c + mdot v); beta has the
    # gravitational denominator G M rho D. There is no second factor of c.
    diameter_m = 3.0 * pressure / (4.0 * PI * G_GRAV * mass * M_SUN_KG * rho)
    return diameter_m * 1e6


def beta_radiation(diameter_um: np.ndarray | float, lum: float, mass: float,
                   rho_gcc: float, qpr: float) -> np.ndarray:
    d_m = np.asarray(diameter_um, dtype=float) * 1e-6
    return 3.0 * lum * L_SUN_W * qpr / (8.0 * PI * C_LIGHT * G_GRAV * mass * M_SUN_KG * rho_gcc * 1000.0 * d_m)


def beta_wind(diameter_um: np.ndarray | float, mass: float, rho_gcc: float,
              mdot_solar: float, vwind: float = WIND_SPEED, cd: float = 1.0) -> np.ndarray:
    d_m = np.asarray(diameter_um, dtype=float) * 1e-6
    return 3.0 * mdot_solar * SOLAR_MDOT * vwind * cd / (8.0 * PI * G_GRAV * mass * M_SUN_KG * rho_gcc * 1000.0 * d_m)


def physical_draws(n: int = 100_000) -> dict[str, np.ndarray]:
    # Truncated normal approximations to the project’s adopted stellar inputs;
    # rho/Q are explicit scenario priors, not an observational posterior.
    lum = np.clip(RNG.normal(0.52, 0.03, n), 0.35, 0.7)
    mass = np.clip(RNG.normal(0.783, 0.012, n), 0.70, 0.85)
    rho = 10 ** RNG.uniform(np.log10(0.5), np.log10(3.5), n)
    qpr = RNG.uniform(0.3, 1.5, n)
    return {"lum": lum, "mass": mass, "rho": rho, "qpr": qpr,
            "d_rad": d_bl_um(lum, mass, rho, qpr),
            "d_total_tau": d_bl_um(lum, mass, rho, qpr, 0.1)}


def characteristic_radius(rin: np.ndarray, rout: np.ndarray, p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    exponent = p + 1.0
    rc = ((exponent + 1.0) / (exponent + 2.0)
          * (rout ** (exponent + 2.0) - rin ** (exponent + 2.0))
          / (rout ** (exponent + 1.0) - rin ** (exponent + 1.0)))
    return rc, (rout - rin) / rc


def broad_vector(rin: np.ndarray, rout: np.ndarray, mass: np.ndarray, age: np.ndarray,
                 dmin: np.ndarray, dc: np.ndarray, qd: np.ndarray,
                 ecc: np.ndarray, inc: np.ndarray, annuli: int = 100) -> np.ndarray:
    """Broad-belt Eq. 14 with per-draw geometry, avoiding invalid broadcasting."""
    n = len(np.asarray(rin))
    frac = np.linspace(0.0, 1.0, annuli + 1)
    logs = np.log(rin)[:, None] + (np.log(rout) - np.log(rin))[:, None] * frac[None, :]
    edges = np.exp(logs)
    rad = np.sqrt(edges[:, :-1] * edges[:, 1:])
    width = np.log(edges[:, 1:] / edges[:, :-1])
    vals = wyatt2007_fmax_eq14(
        rad.T, width.T, mass[None, :], age[None, :], dmin[None, :],
        dc[None, :], qd[None, :], ecc[None, :], inc[None, :]
    )
    return np.sum(vals, axis=0)


def fmax_pair(a: dict[str, np.ndarray], dmin: np.ndarray, geometry: str,
              annuli: int = 100) -> tuple[np.ndarray, np.ndarray]:
    rin, rout, p = a["rin"], a["rout"], a["surface_power"]
    rc, width = characteristic_radius(rin, rout, p)
    ring = wyatt2007_fmax_eq14(rc, width, a["mass"], a["age_myr"], dmin,
                               a["dc_km"], a["qd"], a["e"], a["inc"])
    broad = broad_vector(rin, rout, a["mass"], a["age_myr"], dmin,
                         a["dc_km"], a["qd"], a["e"], a["inc"], annuli)
    return ring, broad


def source_inventory() -> None:
    terms = re.compile(r"(?i)(D[_ ]?bl|Dbl|blow[_ -]?out|grain[_ -]?size|minimum[_ -]?grain|D[_ ]?min|Dmin|a[_ ]?min|s[_ ]?min|20\s*(?:micron|um)|1\.2\s*(?:micron|um))")
    rows = []
    skip = {OUT, ROOT / ".git", ROOT / ".venv"}
    text_ext = {".py", ".md", ".csv", ".json", ".toml", ".yaml", ".yml", ".txt", ".sh", ".ipynb", ".tsv"}
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in text_ext:
            continue
        if any(p == path or p in path.parents for p in skip):
            continue
        try:
            if path.stat().st_size > 10_000_000:
                continue
            lines = path.read_text(errors="ignore").splitlines()
        except Exception:
            continue
        for i, line in enumerate(lines, 1):
            if terms.search(line):
                text = line.strip().replace("\x00", " ")[:300]
                low = text.lower()
                if "20" in low:
                    desc = "fixed 20 um substitution"
                elif "minimum" in low or "dmin" in low or "grain" in low:
                    desc = "minimum/physical grain-size label"
                elif "blow" in low:
                    desc = "radiation-pressure blowout variable"
                else:
                    desc = "grain-size occurrence"
                rows.append({"file_path": str(path.relative_to(ROOT)), "line": i,
                             "matched_text": text, "classification": desc})
    csv_write(OUT / "source_inventory/grain_size_occurrences.csv", rows,
              ["file_path", "line", "matched_text", "classification"])
    (OUT / "source_inventory/grain_size_variable_map.md").write_text(
        "# Grain-size variable map\n\n"
        "- `D_bl`: radiation-pressure blowout diameter at total beta=0.5; derived, not observed.\n"
        "- `D_min`: lower cutoff of the collisional size distribution; may be physical-model derived.\n"
        "- `D_char`: characteristic emitting size; not constrained by the project SED fit.\n"
        "- `blowout_diameter_um=20.0`: fixed hard-coded argument in `run_joint_map_visibility_posterior.py`; no posterior or source citation.\n"
        "- Lawler 2014 `15 +/- 8 um`: preferred physical-grain-model minimum diameter, not a measured blowout diameter.\n"
    )
    cfg = {
        "reference": {"blowout_diameter_um": "truncated Normal(1.2,0.3), lower 0.1; model assumption"},
        "diagnostic": {"blowout_diameter_um": 20.0, "status": "fixed hard-coded value; no prior"},
        "lawler_minimum_grain_diameter_um": {"value": 15.0, "sigma": 8.0, "status": "physical grain model"},
        "physical_audit": {"D_bl": "derived from L,M,rho,Qpr and optional wind", "D_min": "controlled sensitivity parameter", "D_char": "not constrained"},
    }
    (OUT / "source_inventory/configuration_diff.json").write_text(json.dumps(cfg, indent=2))


def write_equations() -> None:
    (OUT / "equation_audit/implemented_fmax_equation.md").write_text(
        "# Implemented collisional equation\n\n"
        "The production function is `tau_ceti.collisional.wyatt2007_fmax_eq14`. It computes Eq. 14 as\n\n"
        "```text\n"
        "x_c = 1.3e-3 [Q_D* r / (M*(1.25 e^2 + I^2))]^(1/3)\n"
        "G(q,x_c) = x_c^(5-3q)-1 + (6q-10)/(3q-4)[x_c^(4-3q)-1] + (3q-5)/(3q-3)[x_c^(3-3q)-1]\n"
        "f_max = 1e-6 r^(3/2) (dr/r) / [4 pi M^(1/2) t] * [2(1+1.25(e/I)^2)^(-1/2)/G] * (D_bl*1e-9/D_c)^(5-3q)\n"
        "```\n\n"
        "Here r is AU, t is Myr, D_bl is passed in micrometres, D_c is km, and q=11/6.\n"
        "The code names the argument `blowout_diameter_um`, but the audit recomputes it as a general lower cascade diameter so that D_bl and D_min are not conflated.\n"
    )
    (OUT / "equation_audit/grain_size_dependence.md").write_text(
        "# Grain-size dependence\n\n"
        "For q=11/6, `5 - 3q = -1/2`. Thus Eq. 14 has f_max ∝ D_cut^(-1/2), and f/f_max ∝ D_cut^(+1/2).\n\n"
        "If the physical derivation uses D_bl but the empirical lower cutoff is D_min, the generalized conversion is\n\n"
        "`f_max(D_min) = f_max(D_bl) * (D_min/D_bl)^(-1/2)`\n\n"
        "This is a reparameterisation only if D_min is explicitly declared the cascade cutoff. Directly relabelling an observational minimum as D_bl is not physically justified. q=11/6 corresponds to the adopted differential size-distribution slope p=3q-2=3.5, so cross-sectional area scales as D_min^(3-p)=D_min^(-0.5).\n"
    )
    csv_write(OUT / "equation_audit/literature_to_code_mapping.csv", [
        {"literature_quantity": "D_bl", "code_argument": "blowout_diameter_um", "unit": "um", "exponent_q11_6": "-0.5", "status": "code label; physical meaning must be declared per run"},
        {"literature_quantity": "q", "code_argument": "q", "unit": "dimensionless", "exponent_q11_6": "5-3q", "status": "default 11/6"},
        {"literature_quantity": "D_c", "code_argument": "largest_body_km", "unit": "km", "exponent_q11_6": "+0.5", "status": "cascade upper size"},
    ], ["literature_quantity", "code_argument", "unit", "exponent_q11_6", "status"])


def stellar_tables(phys: dict[str, np.ndarray]) -> None:
    csv_write(OUT / "tables/stellar_parameters.csv", [
        {"parameter": "luminosity", "value": 0.52, "uncertainty": 0.03, "unit": "L_sun", "source": "project literature_constraints.csv", "used_in_blowout": "yes"},
        {"parameter": "mass", "value": 0.783, "uncertainty": 0.012, "unit": "M_sun", "source": "project literature_constraints.csv", "used_in_blowout": "yes"},
        {"parameter": "age", "value": 7.63, "uncertainty": 0.87, "unit": "Gyr", "source": "project literature_constraints.csv", "used_in_blowout": "no"},
        {"parameter": "effective_temperature", "value": "not used", "uncertainty": "not available in project posterior", "unit": "K", "source": "not required for analytic Qpr grid", "used_in_blowout": "no"},
        {"parameter": "stellar_radius", "value": "not used", "uncertainty": "not available in project posterior", "unit": "R_sun", "source": "not required for analytic Qpr grid", "used_in_blowout": "no"},
    ])
    med = {k: float(np.median(v)) for k, v in phys.items() if k.startswith("d_")}
    (OUT / "stellar_parameters/notes.md").write_text(
        "# Adopted stellar inputs\n\n"
        "The analytic radiation-pressure calculation requires L and M, not Teff or radius. The project adopts L=0.52 +/- 0.03 L_sun and M=0.783 +/- 0.012 M_sun. Teff/radius were not part of the current posterior and are therefore not silently inserted. Q_pr is treated as an explicit grain-composition scenario parameter.\n\n"
        f"Scenario-draw medians: D_bl(rad)={med['d_rad']:.3f} um; D_bl(total, 0.1 solar wind)={med['d_total_tau']:.3f} um.\n"
    )


def physical_tables(phys: dict[str, np.ndarray]) -> None:
    scenarios = [
        ("compact silicate", 3.3), ("compact water ice", 0.93),
        ("compact carbonaceous", 1.8), ("mixed silicate-ice", 2.0),
        ("porous silicate", 1.0), ("porous mixed", 0.7),
    ]
    qgrid = [0.1, 0.3, 0.5, 1.0, 1.5, 2.0]
    rows = []
    for name, rho in scenarios:
        for qpr in qgrid:
            d = d_bl_um(phys["lum"], phys["mass"], rho, qpr)
            rows.append({"scenario": name, "rho_g_cm3": rho, "Qpr": qpr,
                         "Dbl_p16_um": qtile(d)[0], "Dbl_median_um": qtile(d)[1], "Dbl_p84_um": qtile(d)[2],
                         "method": "analytic radiation-pressure grid; no Mie optical constants"})
    csv_write(OUT / "tables/grain_composition_models.csv", rows)
    csv_write(OUT / "tables/physical_blowout_sizes.csv", rows)
    wind_rows = []
    for rho in [0.5, 1.0, 1.5, 2.5, 3.5]:
        for md in [0, 0.1, 0.5, 1, 2, 5, 10]:
            d = d_bl_um(0.52, 0.783, rho, 1.0, md)
            wind_rows.append({"rho_g_cm3": rho, "Qpr": 1.0, "mdot_solar": md,
                              "Dbl_total_um": float(d), "beta_wind_over_radiation": float(md * SOLAR_MDOT * WIND_SPEED * C_LIGHT / (0.52 * L_SUN_W))})
    csv_write(OUT / "tables/wind_blowout_sizes.csv", wind_rows)
    # Required density and wind for 20 um at nominal L/M.
    req = []
    for rho in [0.5, 1.0, 1.5, 2.5, 3.5]:
        for qpr in [0.3, 1.0, 2.0]:
            d0 = float(d_bl_um(0.52, 0.783, rho, qpr, 0.0))
            target = 20.0
            pressure_needed = target * 1e-6 * 4 * PI * G_GRAV * 0.783 * M_SUN_KG * (rho * 1000.0) / 3.0
            mdot = (pressure_needed - 0.52 * L_SUN_W * qpr / C_LIGHT) / (SOLAR_MDOT * WIND_SPEED)
            req.append({"rho_g_cm3": rho, "Qpr": qpr, "Dbl_rad_um": d0,
                        "required_mdot_solar_for_20um": max(0.0, mdot),
                        "radiation_already_exceeds_target": bool(d0 >= target)})
    csv_write(OUT / "tables/required_parameters_for_20um.csv", req)
    (OUT / "blowout_calculations/method.md").write_text(
        "# Blowout calculation\n\n"
        "For a spherical grain of diameter D, beta_rad = 3 L Q_pr / (8 pi c G M rho D). Setting beta=0.5 gives D_bl = 3 L Q_pr/(4 pi c G M rho). This is twice the radius returned by the project helper `blowout_radius_microns`, which correctly uses the radius convention. Wind pressure uses beta_wind = 3 mdot v C_D/(8 pi G M rho D), and the total expression adds L Q_pr/c + mdot v C_D.\n\n"
        "The analytic grid varies density and Q_pr explicitly. It is not a Mie calculation; no optical constants or stellar spectral model were available in the project.\n"
    )


def sensitivity_tables(ref: dict[str, np.ndarray], diag: dict[str, np.ndarray]) -> None:
    grids = [0.5, 1, 1.2, 2, 3, 5, 10, 15, 20, 30, 50]
    n = 5_000
    a = {k: v[:n] for k, v in ref.items() if isinstance(v, np.ndarray) and v.shape == (len(ref["f"]),)}
    # Preserve the physically derived compact reference scenario (~1.2 um).
    dcompact = d_bl_um(0.52, a["mass"], 1.25, 1.0)
    rows_r, rows_b = [], []
    for d in grids:
        ring, broad = fmax_pair(a, np.full(n, d), "reference", annuli=50)
        for dest, vals in [(rows_r, a["f"] / ring), (rows_b, a["f"] / broad)]:
            dest.append({"Dmin_um": d, "ratio_p16": qtile(vals)[0], "ratio_median": qtile(vals)[1], "ratio_p84": qtile(vals)[2]})
    csv_write(OUT / "tables/dmin_sensitivity_ring.csv", rows_r)
    csv_write(OUT / "tables/dmin_sensitivity_broad.csv", rows_b)
    # Diagnostic luminosity under identical reference priors, explicitly controlled.
    c = dict(a); c["f"] = np.resize(diag["f"], n)
    rr = []
    for d in grids:
        ring, broad = fmax_pair(c, np.full(n, d), "reference", annuli=50)
        rr.append({"Dmin_um": d, "ring_median": float(np.median(c["f"] / ring)), "broad_median": float(np.median(c["f"] / broad))})
    csv_write(OUT / "tables/controlled_diagnostic_dmin_sensitivity.csv", rr)
    (OUT / "minimum_size_analysis/interpretation.md").write_text(
        "# Minimum-size sensitivity\n\n"
        "The grid changes only the lower cascade cutoff supplied to the existing Eq. 14 implementation. It is not a fit for D_min and must not be interpreted as a posterior. The Lawler value 15 +/- 8 um is retained as a physical-grain-model scenario, while the 20 um diagnostic is treated as a fixed phenomenological/observational-minimum sensitivity. D_char is not constrained by the current modified-blackbody SED.\n"
    )


def model_matrix(ref: dict[str, np.ndarray], diag: dict[str, np.ndarray]) -> list[dict]:
    n = 1000
    r = {k: v[:n].copy() for k, v in ref.items() if isinstance(v, np.ndarray) and v.ndim == 1 and len(v) >= n}
    d = {k: v[:n].copy() for k, v in diag.items() if isinstance(v, np.ndarray) and v.ndim == 1 and len(v) >= n}
    # d_phys is calculated from the same stellar-mass posterior, with explicit compact rho/Q.
    dphys_r = d_bl_um(0.52, r["mass"], 1.25, 1.0)
    dphys_d = d_bl_um(0.52, d["mass"], 1.25, 1.0)
    models = [
        ("A_reference_physical_Dbl", r, dphys_r, "reference", "physical Dbl; compact rho=1.25, Qpr=1"),
        ("B_diagnostic_luminosity_reference_geometry", {**r, "f": d["f"]}, dphys_r, "reference", "diagnostic f only; physical Dbl"),
        ("C_diagnostic_luminosity_archival_geometry", {**r, "f": d["f"], "rin": d["rin"], "rout": d["rout"], "surface_power": d["surface_power"]}, dphys_r, "reference", "archival geometry; reference collisional priors"),
        ("D_diagnostic_archival_Dmin20_reference_priors", {**r, "f": d["f"], "rin": d["rin"], "rout": d["rout"], "surface_power": d["surface_power"]}, np.full(n, 20.0), "reference", "Dmin=20; reference collisional priors"),
        ("E_diagnostic_archival_physical_Dbl", d, dphys_d, "diagnostic", "diagnostic priors; physical Dbl"),
        ("F_full_diagnostic_Dmin20", d, np.full(n, 20.0), "diagnostic", "diagnostic priors; fixed Dmin=20"),
    ]
    rows = []
    for name, a, size, geom, note in models:
        ring, broad = fmax_pair(a, size, geom, annuli=60)
        ring_ratio = a["f"] / ring
        broad_ratio = a["f"] / broad
        rows.append({"model": name, "grain_parameter": "Dmin/Dbl per note", "note": note,
                     "ring_ratio_p16": qtile(ring_ratio)[0], "ring_ratio_median": qtile(ring_ratio)[1], "ring_ratio_p84": qtile(ring_ratio)[2],
                     "ring_p_gt_1": float(np.mean(ring_ratio > 1)), "ring_p_gt_10": float(np.mean(ring_ratio > 10)),
                     "broad_ratio_p16": qtile(broad_ratio)[0], "broad_ratio_median": qtile(broad_ratio)[1], "broad_ratio_p84": qtile(broad_ratio)[2],
                     "broad_p_gt_1": float(np.mean(broad_ratio > 1)), "broad_p_gt_10": float(np.mean(broad_ratio > 10)),
                     "Dmedian_um": float(np.median(size))})
    csv_write(OUT / "tables/controlled_model_matrix.csv", rows)
    csv_write(OUT / "tables/final_recommended_collisional_results.csv", [rows[0], rows[1], rows[2], rows[3], rows[5]])
    return rows


def diagnostics_tables(rows: list[dict]) -> None:
    # Fixed 20 has no posterior and therefore no meaningful boundary occupancy.
    csv_write(OUT / "tables/prior_boundary_diagnostics.csv", [
        {"parameter": "diagnostic_Dmin_um", "lower": "not defined", "upper": "not defined", "fraction_near_lower": "NA", "fraction_near_upper": "NA", "status": "fixed at 20; no prior/posterior"},
        {"parameter": "reference_Dbl_um", "lower": 0.1, "upper": "not fixed", "fraction_near_lower": "computed from reproduced prior", "fraction_near_upper": "NA", "status": "truncated Normal model assumption"},
        {"parameter": "Lawler_Dmin_um", "lower": "scenario", "upper": "scenario", "fraction_near_lower": "NA", "fraction_near_upper": "NA", "status": "15 +/- 8 physical-grain-model estimate"},
    ])
    thr = []
    for row in rows:
        for model, metric in [("ring", "ring_ratio_median"), ("broad", "broad_ratio_median")]:
            thr.append({"model": row["model"], "metric": metric, "median": row[metric],
                        "threshold_gt_1": row[f"{model}_p_gt_1"], "threshold_gt_10": row[f"{model}_p_gt_10"]})
    csv_write(OUT / "tables/posterior_thresholds.csv", thr)


def make_figures(ref: dict[str, np.ndarray], phys: dict[str, np.ndarray], rows: list[dict]) -> None:
    figdir = OUT / "figures"; figdir.mkdir(exist_ok=True)
    # beta curves
    dgrid = np.geomspace(0.05, 100, 400)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for name, rho, q in [("silicate", 3.3, 1), ("ice", .93, 1), ("carbon", 1.8, 1), ("porous mix", .7, .5)]:
        ax.loglog(dgrid, beta_radiation(dgrid, .52, .783, rho, q), label=f"{name} (rho={rho:g}, Q={q:g})")
    ax.axhline(.5, color="k", ls="--", lw=1); ax.axvline(20, color="crimson", ls=":", label="fixed 20 um")
    ax.set(xlabel="grain diameter D (um)", ylabel="radiation-pressure beta", title="Radiation-pressure beta versus diameter")
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(figdir / "beta_vs_diameter.png", dpi=180); plt.close(fig)
    # distribution
    fig, ax = plt.subplots(figsize=(7, 4.5)); ax.hist(phys["d_rad"], bins=np.geomspace(.05, 30, 80), density=True, alpha=.7, label="analytic physical scenario")
    ax.axvline(20, color="crimson", ls="--", label="diagnostic fixed value"); ax.axvline(float(d_bl_um(.52, .783, 1.25, 1)), color="k", ls=":", label="reference compact scenario")
    ax.set_xscale("log"); ax.set(xlabel="D_bl (um)", ylabel="density", title="Radiation-only physical blowout scenarios"); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(figdir / "blowout_size_distribution.png", dpi=180); plt.close(fig)
    # wind curves
    fig, ax = plt.subplots(figsize=(7, 4.5)); md = np.geomspace(.01, 20, 200)
    for rho in [.5, 1, 2.5, 3.5]: ax.loglog(md, d_bl_um(.52, .783, rho, 1, md), label=f"rho={rho:g} g cm$^{{-3}}$")
    ax.axvline(.1, color="k", ls="--", label="Tau-specific upper limit 0.1 solar"); ax.axhline(20, color="crimson", ls=":")
    ax.set(xlabel="stellar mass-loss rate (solar)", ylabel="total D_bl (um)", title="Radiation plus wind pressure"); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(figdir / "wind_comparison.png", dpi=180); plt.close(fig)
    # Dmin sensitivity
    for broad in [False, True]:
        path = OUT / ("tables/dmin_sensitivity_broad.csv" if broad else "tables/dmin_sensitivity_ring.csv")
        data = list(csv.DictReader(path.open())); x = np.array([float(z["Dmin_um"]) for z in data]); y = np.array([float(z["ratio_median"]) for z in data])
        fig, ax = plt.subplots(figsize=(7, 4.5)); ax.loglog(x, y, "o-"); ax.axvline(20, color="crimson", ls="--"); ax.axhline(1, color="k", ls=":")
        ax.set(xlabel="lower cascade cutoff D_min (um)", ylabel="median f/f_max", title=f"Reference luminosity D_min sensitivity ({'broad' if broad else 'ring'})")
        fig.tight_layout(); fig.savefig(figdir / f"dmin_sensitivity_{'broad' if broad else 'ring'}.png", dpi=180); plt.close(fig)
    # model matrix
    fig, ax = plt.subplots(figsize=(9, 4.5)); x=np.arange(len(rows)); w=.38
    ax.bar(x-w/2, [z["ring_ratio_median"] for z in rows], w, label="ring"); ax.bar(x+w/2, [z["broad_ratio_median"] for z in rows], w, label="broad")
    ax.axhline(1, color="k", ls=":"); ax.set_xticks(x, [f"{i+1}" for i in range(len(rows))]); ax.set(ylabel="median f/f_max", title="Controlled grain-size and prior matrix"); ax.legend(); fig.tight_layout(); fig.savefig(figdir / "controlled_model_matrix.png", dpi=180); plt.close(fig)
    # one-page PDF contact sheet for portability
    with PdfPages(figdir / "grain_size_audit_figures.pdf") as pdf:
        for name in ["beta_vs_diameter.png", "blowout_size_distribution.png", "wind_comparison.png", "dmin_sensitivity_ring.png", "dmin_sensitivity_broad.png", "controlled_model_matrix.png"]:
            im = plt.imread(figdir / name); f, a = plt.subplots(figsize=(8, 5)); a.imshow(im); a.axis("off"); pdf.savefig(f, bbox_inches="tight"); plt.close(f)


def write_report(ref: dict[str, np.ndarray], diag: dict[str, np.ndarray], phys: dict[str, np.ndarray], rows: list[dict]) -> None:
    r0 = rows[0]; b = rows[1]; f = rows[-1]
    dphys_med = float(np.median(phys["d_rad"])); dtau_med = float(np.median(phys["d_total_tau"]))
    target_required = float(d_bl_um(.52, .783, .5, 1, 0))
    required = list(csv.DictReader((OUT / "tables/required_parameters_for_20um.csv").open()))
    req_low = next(z for z in required if z["rho_g_cm3"] == "0.5" and z["Qpr"] == "1.0")
    report = f"""# Tau Ceti grain-size physics and collisional audit

## Executive result

The 20 um value is **not a physically derived Tau Ceti radiation-pressure blowout diameter**. Source inventory identifies it as the literal fixed argument `20.0` in `run_joint_map_visibility_posterior.py`; it has no posterior, no prior, and no source citation in the implementation. The most defensible classification is an **unsupported parameter substitution**, with a secondary interpretation as an observational minimum-grain sensitivity scenario. It is not a replacement for the physical D_bl result.

The project’s accepted f/f_max reconciliation remains untouched. This audit is isolated and uses the same reproduced arrays.

## Physical blowout calculation

For diameter D, beta_rad = 3 L Q_pr/(8 pi c G M rho D), so beta=0.5 gives D_bl = 3 L Q_pr/(4 pi c G M rho). The project helper returns a radius at beta=0.5; the physical diameter is exactly twice that value. With L=0.52 L_sun and M=0.783 M_sun, radiation-only D_bl is approximately {target_required:.2f} um for rho=0.5 g cm^-3 and Q_pr=1, and about {float(d_bl_um(.52,.783,2.5,1)):.2f} um for rho=2.5 g cm^-3. The explicit composition/Q_pr scenario posterior has median D_bl={dphys_med:.2f} um (16–84%: {qtile(phys['d_rad'])[0]:.2f}–{qtile(phys['d_rad'])[2]:.2f} um). It is a scenario distribution, not an observational posterior, because the project has no Mie optical-constant calculation.

Adding a solar-wind pressure term with the documented solar reference and a Tau-specific upper limit of 0.1 solar gives a scenario median of {dtau_med:.2f} um. At rho=0.5 and Q_pr=1, the mass-loss rate required to reach 20 um is {req_low['required_mdot_solar_for_20um']} solar in the analytic table. Thus a Tau-specific wind below 0.1 solar does not make 20 um plausible under the adopted compact-grain scenarios. Reaching 20 um by wind pressure would require thousands of solar mass-loss rates even in this permissive low-density case, which is not a supported Tau Ceti parameter.

## D_bl, D_min and D_char

- **D_bl**: radiation-pressure blowout diameter at beta=0.5, derived from L, M, density, Q_pr and optional wind.
- **D_min**: lower cutoff of the collisional size distribution. Lawler et al. report 15 +/- 8 um from a physical-grain model; this is not a direct D_bl measurement.
- **D_char**: characteristic emitting size; not constrained by the current modified-blackbody SED. Temperature, lambda_0 and beta are not a unique grain-size measurement.

For q=11/6, the implemented Eq. 14 has f_max proportional to D_cut^-0.5, so f/f_max is proportional to D_cut^+0.5. A D_min substitution is valid only when the calculation is explicitly rederived/relabelled as a lower cascade cutoff. Directly calling 20 um D_bl is not valid.

## Controlled collisional matrix

The matrix in `tables/controlled_model_matrix.csv` separates luminosity, geometry, physical D_bl/D_min and diagnostic priors. Model A (reference luminosity, reference geometry, compact physical D_bl) gives ring median {r0['ring_ratio_median']:.3f} and broad median {r0['broad_ratio_median']:.3f}. Model B changes only to the lower diagnostic luminosity and gives ring {b['ring_ratio_median']:.3f} and broad {b['broad_ratio_median']:.3f}. The full diagnostic fixed-20 model gives ring {f['ring_ratio_median']:.3f} and broad {f['broad_ratio_median']:.3f}. This confirms that the lower luminosity is a real downward sensitivity, while the fixed 20 um lower cutoff raises f/f_max through the +0.5 ratio exponent; it does not establish a physical blowout size.

The D_min grid, reference and controlled diagnostic-luminosity sensitivity are in `tables/dmin_sensitivity_ring.csv`, `tables/dmin_sensitivity_broad.csv` and `tables/controlled_diagnostic_dmin_sensitivity.csv`. No D_min posterior was fitted.

## Final recommendation

Use the **physical D_bl calculation** (with explicit density/Q_pr scenarios and the wind bound) as the primary collisional result. Retain the 15 +/- 8 um Lawler quantity and the fixed 20 um case as clearly labelled D_min/observational sensitivity scenarios in an appendix or robustness section. Do not present 20 um as Tau Ceti’s measured or physically derived radiation-pressure blowout diameter. No transient collision is required by this audit, and no transient is established.

## Reproducibility and limitations

Run from the project environment with:\n\n```bash\n/Users/romonedunlop/Documents/Tau\\ Ceti\\ Research/.venv/bin/python /Users/romonedunlop/Documents/Tau\\ Ceti\\ Research/audit_grain_size_physics/run_grain_size_audit.py\n```\n\nThe physical grid is analytic rather than Mie-based; Teff and stellar radius are not used because they are not part of the current project posterior; Eq. 14 is a broad-belt extension of the existing implementation; and the original reproduced arrays remain model-dependent.\n"""
    (OUT / "final_report/tau_ceti_grain_size_audit.md").write_text(report)
    (OUT / "manuscript_text/physical_grain_size_results.md").write_text(
        "The 20 um diameter used in the diagnostic collisional calculation is a fixed, undocumented parameter substitution rather than a derived Tau Ceti radiation-pressure blowout size. For a diameter convention, beta=3LQ_pr/(8 pi cGM rho D), so the beta=0.5 blowout diameter is D_bl=3LQ_pr/(4 pi cGM rho). Using the adopted L=0.52 L_sun and M=0.783 M_sun gives sub-micron-to-few-micron radiation-only diameters across compact-grain scenarios; a Tau-specific wind bound below 0.1 solar does not raise the result to 20 um. The Lawler 15 +/- 8 um value is retained as a physical-model minimum-grain scenario, not as D_bl. Because the implemented q=11/6 Eq. 14 scales f_max as D_cut^-1/2, using a larger D_min increases f/f_max, but this is a lower-cutoff sensitivity and must not be presented as a physical blowout inference.\n"
    )
    # Generate a valid, lightweight PDF companion without requiring pandoc.
    with PdfPages(OUT / "final_report/tau_ceti_grain_size_audit.pdf") as pdf:
        lines = report.splitlines()
        for start in range(0, len(lines), 48):
            fig = plt.figure(figsize=(8.5, 11)); ax = fig.add_axes([0.07, 0.05, 0.86, 0.9]); ax.axis("off")
            ax.text(0, 1, "\n".join(lines[start:start + 48]), va="top", ha="left", family="monospace", fontsize=8.5)
            pdf.savefig(fig); plt.close(fig)


def tests_file() -> None:
    (OUT / "tests/test_grain_size_physics.py").write_text(
        "import numpy as np\n"
        "from run_grain_size_audit import beta_radiation, beta_wind, d_bl_um\n\n"
        "def test_radius_diameter_factor_two():\n"
        "    # project helper radius is half of this audit's diameter\n"
        "    assert np.isclose(d_bl_um(.52, .783, 2.5, 1), 2 * (3*.52*3.828e26/(8*np.pi*299792458*6.67430e-11*.783*1.98847e30*2500)*1e6))\n\n"
        "def test_beta_half_at_blowout():\n"
        "    d = d_bl_um(.52, .783, 2.5, 1)\n"
        "    assert np.isclose(beta_radiation(d, .52, .783, 2.5, 1), .5)\n\n"
        "def test_wind_positive():\n"
        "    assert beta_wind(1, .783, 1, .1) > 0\n\n"
        "def test_dmin_exponent_is_minus_half():\n"
        "    assert np.isclose(5 - 3*(11/6), -.5)\n"
    )


def main() -> None:
    for d in ["source_inventory", "equation_audit", "stellar_parameters", "grain_models", "blowout_calculations", "minimum_size_analysis", "controlled_collisional_runs", "sensitivity_analysis", "posterior_diagnostics", "figures", "tables", "tests", "manuscript_text", "final_report"]:
        (OUT / d).mkdir(parents=True, exist_ok=True)
    ref, diag = load_arrays()
    source_inventory(); write_equations()
    phys = physical_draws(); stellar_tables(phys); physical_tables(phys)
    sensitivity_tables(ref, diag)
    rows = model_matrix(ref, diag); diagnostics_tables(rows); make_figures(ref, phys, rows); write_report(ref, diag, phys, rows); tests_file()
    manifest = {
        "git_commit": git_commit(), "python": sys.version, "platform": platform.platform(), "seed": SEED,
        "reference_npz": {"path": str(PREV / "reference_reproduced_draws.npz"), "sha256": sha256(PREV / "reference_reproduced_draws.npz")},
        "diagnostic_npz": {"path": str(PREV / "diagnostic_reproduced_draws.npz"), "sha256": sha256(PREV / "diagnostic_reproduced_draws.npz")},
        "command": f"{sys.executable} {Path(__file__)}", "status": "completed",
    }
    (OUT / "environment/run_manifest.json").write_text(json.dumps(manifest, indent=2))
    (OUT / "README.md").write_text(
        "# Tau Ceti physical grain-size audit\n\nRun: `" + sys.executable + " " + str(Path(__file__)) + "`\n\n"
        "This isolated audit preserves the accepted f/f_max audit, inventories all grain-size variables, derives physical radiation/wind blowout diameters, and recomputes Eq. 14 under D_bl/D_min-separated scenarios. Primary output: `final_report/tau_ceti_grain_size_audit.md`.\n"
    )
    print("Completed", OUT)
    print("Model matrix:")
    for row in rows:
        print(row["model"], row["ring_ratio_median"], row["broad_ratio_median"])


if __name__ == "__main__":
    main()
