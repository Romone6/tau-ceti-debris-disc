# Blowout calculation

For a spherical grain of diameter D, beta_rad = 3 L Q_pr / (8 pi c G M rho D). Setting beta=0.5 gives D_bl = 3 L Q_pr/(4 pi c G M rho). This is twice the radius returned by the project helper `blowout_radius_microns`, which correctly uses the radius convention. Wind pressure uses beta_wind = 3 mdot v C_D/(8 pi G M rho D), and the total expression adds L Q_pr/c + mdot v C_D.

The analytic grid varies density and Q_pr explicitly. It is not a Mie calculation; no optical constants or stellar spectral model were available in the project.
