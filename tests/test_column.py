import unittest

import numpy as np

from column import potential_temperature, pressure_change, simulate


class ColumnTests(unittest.TestCase):
    def test_dry_adiabatic_displacement_preserves_theta(self):
        temperature = 300.0 * (80000 / 100000) ** (287.0 / 1004.0)
        self.assertAlmostEqual(potential_temperature(temperature, 80000), 300)

    def test_hydrostatic_pressure_tendency(self):
        self.assertEqual(pressure_change(0, 1e10, 21600), 0)
        self.assertAlmostEqual(pressure_change(5e6, 1e10, 21600) / 100, 1.05948)

    def test_budget_and_grid_refinement(self):
        for k in (2, 20):
            errors = []
            for n in (20, 40, 80):
                _, initial, numerical, exact, residual = simulate(n, k)
                self.assertLess(abs(residual), 1e-8)
                self.assertTrue(np.all(numerical >= 0))
                self.assertLess(numerical.sum(), initial.sum())
                errors.append(np.linalg.norm(numerical - exact) / np.sqrt(n))
            self.assertGreater(np.log2(errors[-2] / errors[-1]), 1.8)
