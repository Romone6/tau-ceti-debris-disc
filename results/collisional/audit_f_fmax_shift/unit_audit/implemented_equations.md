# Implemented collisional equation

Both pathways call `tau_ceti.collisional.wyatt2007_fmax_eq14`, with the broad pathway calling `broad_belt_fmax_continuous` over local Eq. 14 annuli.

\[f_{\max}=\frac{10^{-6}r^{1.5}(dr/r)}{4\pi M_*^{0.5}t}\frac{2[1+1.25(e/I)^2]^{-0.5}}{G(q,X_c)}\left(\frac{D_{bl}}{D_c}\right)^{5-3q}.\]

\[X_c=1.3\times10^{-3}\left[\frac{Q_D^*rM_*^{-1}}{1.25e^2+I^2}\right]^{1/3}.\]

The code uses `q=11/6`, `r` and `dr` in AU, age in Myr, stellar mass in solar masses, `D_c` in km, `D_bl` in micrometres, `Q_D*` in J kg^-1, and dimensionless fractional eccentricity/inclination. The diameter ratio is converted with `D_bl_um*1e-9/D_c_km`.

Reference ring: Rc from cross-section weighting over 6–55 AU and full width `(Rout-Rin)/Rc`. Diagnostic ring: Rc from each map/visibility geometry draw and the same full-width convention. Broad models sum local Eq.14 ceilings over logarithmic annuli; they are explicitly a mathematical extension of the narrow-belt equation.
