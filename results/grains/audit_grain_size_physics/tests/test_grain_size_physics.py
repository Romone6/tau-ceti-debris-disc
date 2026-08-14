import numpy as np
from run_grain_size_audit import beta_radiation, beta_wind, d_bl_um

def test_radius_diameter_factor_two():
    # project helper radius is half of this audit's diameter
    assert np.isclose(d_bl_um(.52, .783, 2.5, 1), 2 * (3*.52*3.828e26/(8*np.pi*299792458*6.67430e-11*.783*1.98847e30*2500)*1e6))

def test_beta_half_at_blowout():
    d = d_bl_um(.52, .783, 2.5, 1)
    assert np.isclose(beta_radiation(d, .52, .783, 2.5, 1), .5)

def test_wind_positive():
    assert beta_wind(1, .783, 1, .1) > 0

def test_dmin_exponent_is_minus_half():
    assert np.isclose(5 - 3*(11/6), -.5)
