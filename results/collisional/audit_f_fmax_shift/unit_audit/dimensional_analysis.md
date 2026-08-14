# Dimensional and range audit

The implementation is dimensionless at the output: `fmax` and `f/fmax` are fractional luminosities. Unit-bearing inputs are converted at the API boundary by convention, not by implicit Astropy coercion.

| Input | Required convention | Reference | Diagnostic |
|---|---|---|---|
| age | Myr | 4000–11000 effective support | 4000–11000 clipped Normal |
| radius | AU | 6–55 | posterior-derived AU |
| D_c | km | 1–2000 | 10–2000 |
| D_bl | micrometres | truncated Normal near 1.2 | fixed 20 |
| Q_D* | J kg^-1 | 50–500 | 10–1000 |
| M* | solar masses | near 0.783 | fixed 0.78 |
| e | fraction | 0.01–0.1 | 0.01–0.2 |
| I | dimensionless radian proxy | I=e | I=0.5e, clipped at 0.001 |

Automated checks reject non-positive age/radius/mass/diameters/energy, non-positive width, eccentricity/inclination outside (0,1], and Rout <= Rin. Explicit factor tests for Myr/Gyr, km/m, and percentage/fraction conversion are in `tests/test_reconciliation.py`.
