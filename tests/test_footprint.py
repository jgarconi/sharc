# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 08:53:15 2017

@author: Calil

Updated on Sun Jul 29 23:17:07 2024

- Added unit tests for LEO satellites at 600 km and 1200 km heights.
- Created tests for plotting footprints, elevation vs. area, and area print.
- TODO: Update 'plot_footprints_print_areas' and 'plot_elevation_vs_area' to
        handle different beam widths.
- TODO: Resolve 'RuntimeWarning: invalid value encountered in arcsin:
        beta_n = arcsin((1/self.sigma)*sin(gamma_n)) - gamma_n'

@author: Juliana Garçoni
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.support.footprint import Footprint


class FootprintAreaTest(unittest.TestCase):

    def setUp(self):
        # Initialize footprints with different parameters
        self.fa1 = Footprint(0.1, bore_lat_deg=0, bore_subsat_long_deg=0.0)
        self.fa2 = Footprint(0.325, bore_lat_deg=0)
        self.fa3 = Footprint(0.325, elevation_deg=20)
        self.fa4 = Footprint(0.325, bore_lat_deg=5, bore_subsat_long_deg=10)
        self.fa5 = Footprint(5, elevation_deg=20, sat_height=1200000)
        self.fa6 = Footprint(5, bore_lat_deg=0, bore_subsat_long_deg=0, sat_height=600000)

    def test_construction(self):
        # Verify test properties with only beam width and boresight coordinates provided
        self.assertEqual(self.fa1.sat_height, 35786000)
        self.assertEqual(self.fa1.bore_lat_deg, 0)
        self.assertEqual(self.fa1.bore_subsat_long_deg, 0)
        self.assertEqual(self.fa1.beam_width_deg, 0.1)
        self.assertEqual(self.fa1.elevation_deg, 90)
        self.assertEqual(self.fa1.bore_lat_rad, 0)
        self.assertEqual(self.fa1.bore_subsat_long_rad, 0)
        self.assertEqual(self.fa1.beam_width_rad, np.pi / 1800)
        self.assertEqual(self.fa1.beta, 0)
        self.assertAlmostEqual(self.fa1.sigma, 0.15, delta=0.01)
        self.assertEqual(self.fa1.bore_tilt, 0)
        self.assertAlmostEqual(self.fa1.elevation_rad, 1.57, delta=0.01)

        # Verify test properties with only beam width and latitude provided
        self.assertEqual(self.fa2.bore_lat_deg, 0)
        self.assertEqual(self.fa2.bore_subsat_long_deg, 0)
        self.assertEqual(self.fa2.bore_lat_rad, 0)
        self.assertEqual(self.fa2.bore_subsat_long_rad, 0)

        # Verify test properties with only beam width and elevation provided
        self.assertEqual(self.fa3.sat_height, 35786000)
        self.assertEqual(self.fa3.elevation_deg, 20)
        self.assertEqual(self.fa3.bore_lat_deg, 0)
        self.assertAlmostEqual(self.fa3.bore_subsat_long_deg, 61.84, delta=0.01)
        self.assertAlmostEqual(self.fa3.elevation_rad, 0.35, delta=0.01)
        self.assertAlmostEqual(self.fa3.bore_subsat_long_rad, 1.08, delta=0.01)

        # Verify test properties with non-zero boresight coordinates
        self.assertEqual(self.fa4.bore_lat_deg, 5)
        self.assertEqual(self.fa4.bore_subsat_long_deg, 10)
        self.assertEqual(self.fa4.beam_width_deg, 0.325)
        self.assertAlmostEqual(self.fa4.beta, 0.2, delta=0.01)
        self.assertAlmostEqual(self.fa4.sigma, 0.15, delta=0.01)
        self.assertAlmostEqual(self.fa4.bore_tilt, 0.034, delta=0.001)
        self.assertAlmostEqual(self.fa4.elevation_deg, 76.86, delta=0.01)

        # Verify test properties with non-geostationary satellite (elevation provided)
        self.assertEqual(self.fa5.sat_height, 1200000)
        self.assertEqual(self.fa5.bore_lat_deg, 0)
        self.assertAlmostEqual(self.fa5.bore_subsat_long_deg, 17.75, delta=0.01)
        self.assertAlmostEqual(self.fa5.sigma, 0.84, delta=0.01)

        # Verify test properties with non-geostationary satellite (elevation not provided)
        self.assertEqual(self.fa6.sat_height, 600000)
        self.assertAlmostEqual(self.fa6.sigma, 0.91, delta=0.01)

    def test_set_elevation(self):
        self.fa2.set_elevation(20)
        self.assertEqual(self.fa2.bore_lat_deg, 0)
        self.assertAlmostEqual(self.fa2.bore_subsat_long_deg, 61.84, delta=0.01)

    def test_calc_footprint(self):
        fp_long, fp_lat = self.fa1.calc_footprint(4)
        npt.assert_allclose(fp_long, np.array([0.0, 0.487, -0.487, 0.0]), atol=1e-2)
        npt.assert_allclose(fp_lat, np.array([-0.562, 0.281, 0.281, -0.562]), atol=1e-2)

    def test_calc_area(self):
        # Test area calculation with 0.3% accuracy
        a1 = self.fa2.calc_area(1000)
        self.assertAlmostEqual(a1, 130000, delta=130000 * 0.003)
        a2 = self.fa3.calc_area(1000)
        self.assertAlmostEqual(a2, 486300, delta=486300 * 0.003)
        a3 = self.fa6.calc_area(1000)
        self.assertAlmostEqual(a3, 8700, delta=8700 * 0.003)

    def test_plot_footprints_print_areas(self):
        # Test footprint plotting and area printing for different heights and elevations
        self.fa6.plot_footprints_print_areas(
            beam_width=5,
            heights=[600000, 1200000],
            elevations=[90, 45, 30, 20, 10],
            n_el=100,
            n_poly=1000
        )

    def test_plot_elevation_vs_area(self):
        # Test elevation vs area plotting for different heights
        self.fa6.plot_elevation_vs_area(
            beam_width=5,
            heights=[600000, 1200000],
            n_el=100,
            n_poly=1000
        )


if __name__ == '__main__':
    unittest.main()
