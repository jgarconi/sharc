# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 08:52:28 2017

@author: Calil

Updated on Sun Jul 29 23:17:07 2024

- Updated '__init__' to accept different heights when elevation isn't provided.
- Created 'plot_footprints_print_areas' and 'plot_elevation_vs_area' functions.
- Removed individual tests; all tests are now in 'test_footprint.py'.
- TODO: Update 'plot_footprints_print_areas' and 'plot_elevation_vs_area' to
        handle different beam widths.
- TODO: Combine 'plot_footprints_print_areas' and 'plot_elevation_vs_area' into
        a single function.
- TODO: Resolve 'RuntimeWarning: invalid value encountered in arcsin:
        beta_n = arcsin((1/self.sigma)*sin(gamma_n)) - gamma_n'

@author: Juliana Garçoni
"""

from area import area as earth_area
from numpy import (cos, sin, tan, arctan, deg2rad, rad2deg, arccos,
                   pi, linspace, arcsin, vstack, arctan2, where, zeros_like)
import matplotlib.pyplot as plt
from sharc.parameters.constants import EARTH_RADIUS


class Footprint(object):
    """
    Defines a satellite footprint region and calculates its area.
    Method for generating footprints (Siocos, 1973) is found in the book
    "Satellite Communication Systems" by M. Richharia ISBN 0-07-134208-7

    Construction:
        Footprint(bore_lat_deg, bore_subsat_long_deg, beam)
            beam_deg (float): half of beam width in degrees
            elevation_deg (float): optional. Satellite elevation at
            boresight bore_lat_deg (float): optional, default = 0.
                Latitude of boresight point. If elevation is given this
                parameter is not used. Default = 0
            bore_subsat_long_deg (float): longitude of boresight with respect
                to sub-satellite point, taken positive when to the west of the
                sub-satellite point. If elevation is given this
                parameter is not used. Default = 0
            sat_height (int): optional, Default = 35786000.
                Height of satellite in meters. If none are given, it is assumed
                that it is a geostationary satellite.
    """
    def __init__(self, beam_deg: float, **kwargs):
        # Initialize attributes
        self.bore_lat_deg = 0.0
        self.bore_subsat_long_deg = 0.0
        self.sat_height = 35786000

        if 'elevation_deg' in kwargs and 'sat_height' in kwargs:
            self.elevation_deg = kwargs['elevation_deg']
            self.sat_height = kwargs['sat_height']
            self.sigma = EARTH_RADIUS / (EARTH_RADIUS + self.sat_height)
            self.bore_subsat_long_deg = self.calc_beta(self.elevation_deg)
        elif 'elevation_deg' in kwargs:
            self.elevation_deg = kwargs['elevation_deg']
            self.sigma = EARTH_RADIUS / (EARTH_RADIUS + self.sat_height)
            self.bore_subsat_long_deg = self.calc_beta(self.elevation_deg)
        else:
            if 'sat_height' in kwargs:
                self.sat_height = kwargs['sat_height']
            if 'bore_lat_deg' in kwargs:
                self.bore_lat_deg = kwargs['bore_lat_deg']
            if 'bore_subsat_long_deg' in kwargs:
                self.bore_subsat_long_deg = kwargs['bore_subsat_long_deg']
            self.sigma = EARTH_RADIUS / (EARTH_RADIUS + self.sat_height)
            self.elevation_deg = self.calc_elevation(
                self.bore_lat_deg, self.bore_subsat_long_deg
            )

        self.beam_width_deg = beam_deg

        # Convert to radians
        self.elevation_rad = deg2rad(self.elevation_deg)
        self.bore_lat_rad = deg2rad(self.bore_lat_deg)
        self.bore_subsat_long_rad = deg2rad(self.bore_subsat_long_deg)
        self.beam_width_rad = deg2rad(self.beam_width_deg)

        # Calculate tilt
        self.beta = arccos(
            cos(self.bore_lat_rad) * cos(self.bore_subsat_long_rad)
        )
        self.bore_tilt = arctan2(
            sin(self.beta), (1 / self.sigma - cos(self.beta))
        )

        # Maximum tilt and latitude coverage
        self.max_beta_rad = arccos(self.sigma)
        self.max_gamma_rad = pi / 2 - self.max_beta_rad

    def calc_beta(self, elev_deg: float) -> float:
        """
        Calculates elevation angle based on given elevation. Beta is the
        subpoint to earth station great-circle distance

        Input:
            elev_deg (float): elevation in degrees

        Output:
            beta (float): beta angle in degrees
        """
        elev_rad = deg2rad(elev_deg)
        beta = 90 - elev_deg - rad2deg(arcsin(
            cos(elev_rad) * self.sigma
        ))
        return beta

    def calc_elevation(self, lat_deg: float, long_deg: float) -> float:
        """
        Calculates elevation for given latitude of boresight point and
        longitude of boresight with respect to sub-satellite point.

        Inputs:
            lat_deg (float): latitude of boresight point in degrees
            long_deg (float): longitude of boresight with respect
                to sub-satellite point, taken positive when to the west of the
                sub-satellite point, in degrees

        Output:
            elev (float): elevation in degrees
        """
        lat_rad = deg2rad(lat_deg)
        long_rad = deg2rad(long_deg)
        beta = arccos(cos(lat_rad) * cos(long_rad))
        elev = arctan2(cos(beta) - self.sigma, sin(beta))

        return rad2deg(elev)

    def set_elevation(self, elev: float):
        """
        Resets elevation angle to given value
        """
        self.elevation_deg = elev
        self.bore_lat_deg = 0.0
        self.bore_subsat_long_deg = self.calc_beta(self.elevation_deg)

        # Convert to radians
        self.elevation_rad = deg2rad(self.elevation_deg)
        self.bore_lat_rad = deg2rad(self.bore_lat_deg)
        self.bore_subsat_long_rad = deg2rad(self.bore_subsat_long_deg)

        # Calculate tilt
        self.beta = arccos(
            cos(self.bore_lat_rad) * cos(self.bore_subsat_long_rad)
        )
        self.bore_tilt = arctan2(
            sin(self.beta), (1 / self.sigma - cos(self.beta))
        )

    def calc_footprint(self, n: int):
        """
        Defines footprint polygonal approximation

        Input:
            n (int): number of vertices on polygonal

        Outputs:
            pt_long (np.array): longitude of vertices in degrees
            pt_lat (np.array): latitude of vertices in degrees
        """
        # Projection angles
        phi = linspace(0, 2 * pi, num=n)

        cos_gamma_n = (cos(self.bore_tilt) * cos(self.beam_width_rad) +
                       sin(self.bore_tilt) * sin(self.beam_width_rad) *
                       cos(phi))

        gamma_n = arccos(cos_gamma_n)
        phi_n = arctan2(
            sin(phi),
            (sin(self.bore_tilt) * self.cot(self.beam_width_rad) -
             cos(self.bore_tilt) * cos(phi))
        )

        eps_n = arctan2(sin(self.bore_subsat_long_rad), tan(self.bore_lat_rad)) + phi_n

        beta_n = arcsin((1 / self.sigma) * sin(gamma_n)) - gamma_n
        beta_n[where(gamma_n > self.max_gamma_rad)] = self.max_beta_rad

        pt_lat = arcsin(sin(beta_n) * cos(eps_n))
        pt_long = arctan(tan(beta_n) * sin(eps_n))

        return rad2deg(pt_long), rad2deg(pt_lat)

    def calc_area(self, n: int) -> float:
        """
        Returns footprint area in km^2

        Input:
            n (int): number of vertices on polygonal approximation

        Output:
            a (float): footprint area in km^2
        """
        long, lat = self.calc_footprint(n)
        long_lat = vstack((long, lat)).T

        obj = {'type': 'Polygon', 'coordinates': [long_lat.tolist()]}

        return earth_area(obj) * 1e-6

    def cot(self, angle):
        return tan(pi / 2 - angle)

    def arccot(self, x):
        return pi / 2 - arctan(x)

    def plot_footprints_print_areas(self, beam_width, heights, elevations, n_el, n_poly):
        colors = ['k', 'b', 'r', 'g', 'y']

        # Plot footprints for different elevations
        for height in heights:
            plt.figure(figsize=(15, 2))

            for i, elevation in enumerate(elevations):
                footprint = Footprint(beam_width, elevation_deg=elevation, sat_height=height)
                long, lat = footprint.calc_footprint(n_el)
                plt.plot(long, lat, colors[i], label=f"{elevation}º")

                # Print areas
                area = footprint.calc_area(n_poly)
                print(f"Sat at {height // 1000}km elevation {elevation}º: area = {area}")

            plt.title(f"Footprints at {height // 1000}km")
            plt.legend(loc='upper right')
            plt.xlabel('Longitude [deg]')
            plt.ylabel('Latitude [deg]')
            plt.grid()
            plt.show()

    def plot_elevation_vs_area(self, beam_width, heights, n_el, n_poly):
        el_range = linspace(0, 90, num=n_el)
        plt.figure(figsize=(15, 4))

        # Plot area vs elevation
        for height in heights:
            area = zeros_like(el_range)
            footprint = Footprint(beam_width, elevation_deg=0, sat_height=height)
            for k in range(len(el_range)):
                footprint.set_elevation(el_range[k])
                area[k] = footprint.calc_area(n_poly)
            plt.plot(el_range, area, label=f"{height // 1000}km")

        plt.xlabel('Elevation [deg]')
        plt.ylabel('Footprint area [$km^2$]')
        plt.legend(loc='upper right')
        plt.xlim([0, 90])
        plt.grid()
        plt.show()
