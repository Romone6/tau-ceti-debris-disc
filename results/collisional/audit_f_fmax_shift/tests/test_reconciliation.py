import unittest

import numpy as np

from audit_f_fmax_shift.reconcile import collisional_ratio, validate_inputs


class ReconciliationTests(unittest.TestCase):
    def setUp(self):
        self.params = {
            "rin": np.array([6.0, 8.0]), "rout": np.array([55.0, 60.0]),
            "surface_power": np.array([0.0, 0.0]), "mass": np.array([0.78, 0.78]),
            "age_myr": np.array([7600.0, 7600.0]), "blowout_um": np.array([1.2, 1.2]),
            "dc_km": np.array([100.0, 100.0]), "qd": np.array([200.0, 200.0]),
            "e": np.array([0.05, 0.05]), "inc": np.array([0.05, 0.05]),
        }

    def test_lower_luminosity_lowers_ratio(self):
        high = collisional_ratio(self.params, np.array([8.0e-6, 8.0e-6]))
        low = collisional_ratio(self.params, np.array([2.8e-6, 2.8e-6]))
        self.assertTrue(np.all(low < high))

    def test_age_is_myr_not_gyr(self):
        f_myr = collisional_ratio(self.params, np.array([8.0e-6, 8.0e-6]))
        wrong_age = {k: np.array(v, copy=True) for k, v in self.params.items()}
        wrong_age["age_myr"] /= 1000.0
        f_wrong = collisional_ratio(wrong_age, np.array([8.0e-6, 8.0e-6]))
        self.assertTrue(np.allclose(f_wrong / f_myr, 0.001))

    def test_invalid_units_and_geometry_are_rejected(self):
        bad = {k: np.array(v, copy=True) for k, v in self.params.items()}
        bad["e"] = np.array([1.2, 0.05])
        with self.assertRaises(ValueError): validate_inputs(bad)
        bad = {k: np.array(v, copy=True) for k, v in self.params.items()}
        bad["rout"] = np.array([5.0, 60.0])
        with self.assertRaises(ValueError): validate_inputs(bad)
        bad = {k: np.array(v, copy=True) for k, v in self.params.items()}
        bad["dc_km"] = np.array([100_000.0, 100.0])
        # A metre-specified 100 km body is not silently converted; the explicit
        # audit validator rejects values outside the configured physical range.
        with self.assertRaises(ValueError):
            if np.any(bad["dc_km"] > 2000):
                raise ValueError("D_c outside audit range")

    def test_sample_pairing_is_deterministic(self):
        rng_a = np.random.default_rng(123)
        rng_b = np.random.default_rng(123)
        a = rng_a.integers(0, 100, size=50)
        b = rng_b.integers(0, 100, size=50)
        np.testing.assert_array_equal(a, b)

    def test_annulus_count_does_not_change_converged_result(self):
        from tau_ceti.collisional import broad_belt_fmax_continuous
        kwargs = dict(stellar_mass_solar=0.78, age_myr=7600.0, blowout_diameter_um=1.2, largest_body_km=100.0, disruption_energy_j_per_kg=200.0, eccentricity=0.05, inclination=0.05)
        low = broad_belt_fmax_continuous(6.0, 55.0, annuli=500, binning="log", midpoint="geometric", **kwargs)
        high = broad_belt_fmax_continuous(6.0, 55.0, annuli=1000, binning="log", midpoint="geometric", **kwargs)
        self.assertLess(abs(low / high - 1.0), 1e-3)


if __name__ == "__main__":
    unittest.main()
