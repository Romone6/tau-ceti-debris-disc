# Implemented collisional equation

The production function is `tau_ceti.collisional.wyatt2007_fmax_eq14`. It computes Eq. 14 as

```text
x_c = 1.3e-3 [Q_D* r / (M*(1.25 e^2 + I^2))]^(1/3)
G(q,x_c) = x_c^(5-3q)-1 + (6q-10)/(3q-4)[x_c^(4-3q)-1] + (3q-5)/(3q-3)[x_c^(3-3q)-1]
f_max = 1e-6 r^(3/2) (dr/r) / [4 pi M^(1/2) t] * [2(1+1.25(e/I)^2)^(-1/2)/G] * (D_bl*1e-9/D_c)^(5-3q)
```

Here r is AU, t is Myr, D_bl is passed in micrometres, D_c is km, and q=11/6.
The code names the argument `blowout_diameter_um`, but the audit recomputes it as a general lower cascade diameter so that D_bl and D_min are not conflated.
