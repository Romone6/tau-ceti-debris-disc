# Reference pipeline inventory

Entry point: `run_collisional_geometry_audit.py`.

The reference posterior is sampled from `results/sed/single_mbb_posterior.npz` using seed 20260723 and resampled to 200,000 draws. Age is truncated Normal(7630,870) Myr with lower bound 4000; stellar mass is truncated Normal(0.783,0.012) solar masses; blowout diameter is truncated Normal(1.2,0.3) micrometres; D_c is log-uniform 1–2000 km; Q_D* is log-uniform 50–500 J kg^-1; eccentricity is log-uniform 0.01–0.1 and inclination equals eccentricity. Geometry is 6–55 AU with surface-density power 0. The ring uses a cross-section-weighted characteristic radius and full width; the broad model sums Eq. 14 local ceilings over 250 logarithmic annuli.
