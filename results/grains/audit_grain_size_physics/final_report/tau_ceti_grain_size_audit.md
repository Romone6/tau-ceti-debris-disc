# Tau Ceti grain-size physics and collisional audit

## Executive result

The 20 um value is **not a physically derived Tau Ceti radiation-pressure blowout diameter**. Source inventory identifies it as the literal fixed argument `20.0` in `run_joint_map_visibility_posterior.py`; it has no posterior, no prior, and no source citation in the implementation. The most defensible classification is an **unsupported parameter substitution**, with a secondary interpretation as an observational minimum-grain sensitivity scenario. It is not a replacement for the physical D_bl result.

The project’s accepted f/f_max reconciliation remains untouched. This audit is isolated and uses the same reproduced arrays.

## Physical blowout calculation

For diameter D, beta_rad = 3 L Q_pr/(8 pi c G M rho D), so beta=0.5 gives D_bl = 3 L Q_pr/(4 pi c G M rho). The project helper returns a radius at beta=0.5; the physical diameter is exactly twice that value. With L=0.52 L_sun and M=0.783 M_sun, radiation-only D_bl is approximately 3.05 um for rho=0.5 g cm^-3 and Q_pr=1, and about 0.61 um for rho=2.5 g cm^-3. The explicit composition/Q_pr scenario posterior has median D_bl=0.95 um (16–84%: 0.46–2.06 um). It is a scenario distribution, not an observational posterior, because the project has no Mie optical-constant calculation.

Adding a solar-wind pressure term with the documented solar reference and a Tau-specific upper limit of 0.1 solar gives a scenario median of 0.95 um. At rho=0.5 and Q_pr=1, the mass-loss rate required to reach 20 um is 7319.216376895384 solar in the analytic table. Thus a Tau-specific wind below 0.1 solar does not make 20 um plausible under the adopted compact-grain scenarios. Reaching 20 um by wind pressure would require thousands of solar mass-loss rates even in this permissive low-density case, which is not a supported Tau Ceti parameter.

## D_bl, D_min and D_char

- **D_bl**: radiation-pressure blowout diameter at beta=0.5, derived from L, M, density, Q_pr and optional wind.
- **D_min**: lower cutoff of the collisional size distribution. Lawler et al. report 15 +/- 8 um from a physical-grain model; this is not a direct D_bl measurement.
- **D_char**: characteristic emitting size; not constrained by the current modified-blackbody SED. Temperature, lambda_0 and beta are not a unique grain-size measurement.

For q=11/6, the implemented Eq. 14 has f_max proportional to D_cut^-0.5, so f/f_max is proportional to D_cut^+0.5. A D_min substitution is valid only when the calculation is explicitly rederived/relabelled as a lower cascade cutoff. Directly calling 20 um D_bl is not valid.

## Controlled collisional matrix

The matrix in `tables/controlled_model_matrix.csv` separates luminosity, geometry, physical D_bl/D_min and diagnostic priors. Model A (reference luminosity, reference geometry, compact physical D_bl) gives ring median 0.187 and broad median 0.231. Model B changes only to the lower diagnostic luminosity and gives ring 0.062 and broad 0.077. The full diagnostic fixed-20 model gives ring 0.664 and broad 0.726. This confirms that the lower luminosity is a real downward sensitivity, while the fixed 20 um lower cutoff raises f/f_max through the +0.5 ratio exponent; it does not establish a physical blowout size.

The D_min grid, reference and controlled diagnostic-luminosity sensitivity are in `tables/dmin_sensitivity_ring.csv`, `tables/dmin_sensitivity_broad.csv` and `tables/controlled_diagnostic_dmin_sensitivity.csv`. No D_min posterior was fitted.

## Final recommendation

Use the **physical D_bl calculation** (with explicit density/Q_pr scenarios and the wind bound) as the primary collisional result. Retain the 15 +/- 8 um Lawler quantity and the fixed 20 um case as clearly labelled D_min/observational sensitivity scenarios in an appendix or robustness section. Do not present 20 um as Tau Ceti’s measured or physically derived radiation-pressure blowout diameter. No transient collision is required by this audit, and no transient is established.

## Reproducibility and limitations

Run from the project environment with:

```bash
/Users/romonedunlop/Documents/Tau\ Ceti\ Research/.venv/bin/python /Users/romonedunlop/Documents/Tau\ Ceti\ Research/audit_grain_size_physics/run_grain_size_audit.py
```

The physical grid is analytic rather than Mie-based; Teff and stellar radius are not used because they are not part of the current project posterior; Eq. 14 is a broad-belt extension of the existing implementation; and the original reproduced arrays remain model-dependent.
