# Collisional geometry audit

## Finding

The previous broad-annulus result was numerically invalid. It divided a sum of width-dependent local ceilings by an area normalisation. Because Eq. 14 contains `dr/r`, that operation makes the reported total decrease as annuli are made narrower. It explains the former broad-belt median near 6 and is not a defensible physical ceiling.

## Equation and units

`literature/equations.md` now maps Eq. 14 directly to code. The source defines `r` and `dr` in au, age in Myr, `D_c` as a *diameter* in km, `D_bl` as a diameter in μm, and `Q_D*` in J kg^-1. The implementation converts μm/km before taking the diameter ratio. Eq. 14 is a narrow belt expression; its width is `dr/r`. It contains no independent stellar-luminosity exponent except through `D_bl` (Eq. 5).

## Formulations

A uses a cross-section-weighted characteristic radius R_c=37.06 au and full Δr/R_c. B integrates local narrow-belt ceilings over 6–55 au. C assumes independently evolving zones but sums their luminosity ceilings, so it is algebraically identical to B when the same Eq. 14 local prescription is imposed. A different surface-density profile changes allocation among zones, but not Σfmax_i.

## Convergence and diagnosis

Across binning and midpoint conventions, the largest difference for at least 100 annuli is 0.010%. The corrected continuous integral and independent-annuli sum agree by construction; any remaining difference from A is geometry (characteristic radius and full-width approximation), together with interactions in Eq. 14's radius-dependent collision threshold. Rank correlations with Δ=log10[(f/fmax)_B/(f/fmax)_A] are: age_myr=0.004, stellar_mass_solar=0.009, blowout_diameter_um=-0.001, largest_body_km=0.003, disruption_energy_j_per_kg=-0.423, eccentricity=0.900. Shared f_obs and age cancel algebraically from Δ, so they are not drivers.

## Classification

- **A_characteristic_radius_ring**: median f/fmax=0.19; 68% 0.037–0.94; 95% 0.0096–3.6; P(>1)=0.150; P(>10)=0.001; **Consistent with steady state**.
- **B_continuous_broad_belt**: median f/fmax=0.23; 68% 0.045–1.2; 95% 0.012–4.4; P(>1)=0.183; P(>10)=0.003; **Marginal or assumption-dependent**.
- **C_independent_annuli**: median f/fmax=0.23; 68% 0.045–1.2; 95% 0.012–4.4; P(>1)=0.183; P(>10)=0.003; **Marginal or assumption-dependent**.

The evidence permits steady-state evolution. It does not establish a transient event under any valid formulation. The broad-belt conclusion remains assumption-dependent because the narrow-belt analytic ceiling has been extrapolated over a wide resolved belt; it must not be compressed into one geometry-independent number.
