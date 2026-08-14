# Tau Ceti f/fmax reconciliation report

## Direct answers

1. **Why did f/fmax increase while f decreased?** The lower luminosity alone moves the ratio downward. The increase is produced by a substantially lower diagnostic fmax, primarily from its 20-um blowout-diameter assumption, the changed eccentricity/inclination prescription, changed geometry/width, and broadened collisional priors. The decomposition below quantifies each contribution.
2. **What caused the lower fmax?** The implemented Eq. 14 is shared and unit-consistent. The dominant direct input change is the diagnostic blowout diameter (20 um versus a reference median near 1.2 um); because the exponent is negative for q=11/6, this lowers fmax. Geometry and collision-parameter changes add further effects.
3. **Was an error found?** No unit-conversion, posterior-pairing, or annulus-resolution error was found in the reproduced calculations. The diagnostic posterior is nevertheless boundary-dominated in inclination/position angle and uses a Laplace approximation with correlated-noise inflation.
4. **Which result belongs in the manuscript?** Retain the accepted baseline as primary. Report the archival joint result only as a diagnostic sensitivity analysis, with its provisional calibration and prior dependence stated explicitly.
5. **Reference, diagnostic, or both?** Both should be retained as model-dependent alternatives; the reference remains the headline result until a validated HIPE/original-CASA joint posterior exists.
6. **Steady state or transient?** The archival diagnostic raises the probability of collisional tension but does not identify a transient event. No independent clump, asymmetry, or impact signature establishes a recent collision.
7. **Justified language:** The archival analysis increases the posterior probability of collisional tension, but the magnitude depends on geometry and physical priors. No recent catastrophic collision is positively identified.

## Reproduction

Reference ring median f/fmax = 0.18608; diagnostic ring median = 0.673886; observed shift = 3.621x.
Reference ring median fmax = 4.51051e-05; diagnostic ring median fmax = 4.21138e-06; diagnostic/reference fmax = 0.0933682x.
The production diagnostic ratio is reproduced from the stored diagnostic draws after reconstructing its omitted collisional random variables from the recorded seed and source order.

## Interaction-aware multiplicative decomposition

The permutation-averaged Shapley decomposition uses matched 400-draw samples. The factors below multiply to the matched total shift; the residual is numerical roundoff.
- f: ×0.3462 (Δlog10=-0.4606)
- geometry: ×1.669 (Δlog10=0.2225)
- age_myr: ×0.9908 (Δlog10=-0.0040)
- mass: ×0.9961 (Δlog10=-0.0017)
- blowout_um: ×4.17 (Δlog10=0.6202)
- dc_km: ×0.5814 (Δlog10=-0.2355)
- qd: ×1.311 (Δlog10=0.1175)
- e: ×1.834 (Δlog10=0.2635)
- inc: ×1.131 (Δlog10=0.0534)
- width: ×1 (Δlog10=0.0000)
- Matched total: ×3.76; Shapley residual: 1.11e-16 dex.

## Dominant variables

The reference implementation uses Wyatt et al. (2007) Eq. 14 with r and dr in AU, age in Myr, D_c in km, D_bl in micrometres and Q_D* in J kg^-1. The code converts D_bl/D_c with 1e-9. Both pathways call the same implementation.

Top posterior correlations with diagnostic log10(f/fmax):

- e: Spearman rho=0.698, Pearson r=0.695
- inc: Spearman rho=0.698, Pearson r=0.695
- qd: Spearman rho=-0.529, Pearson r=-0.541
- dc_km: Spearman rho=-0.413, Pearson r=-0.425
- rout: Spearman rho=-0.183, Pearson r=-0.195

## Numerical and validation boundaries

Broad-belt convergence was evaluated for 10–1000 annuli, linear/log spacing, and arithmetic/geometric midpoints. The maximum relative error at production-scale resolution is 0.0026% in the median-input test.

The diagnostic joint fit has reduced chi-square 0.370 after a factor-of-two Herschel noise inflation and uses bounded Laplace draws; this is a diagnostic likelihood, not a validated raw-telemetry posterior.

Machine-readable tables and figures are in the sibling audit directories.
