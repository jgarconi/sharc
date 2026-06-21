# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 15:35:51 2017
Versão Sincronizada com MATLAB - Média Rigorosa Vetorizada
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.io

from sharc.antenna.antenna_element_imt_m2101 import AntennaElementImtM2101
from sharc.antenna.antenna_element_imt_f1336 import AntennaElementImtF1336
from sharc.antenna.antenna_element_imt_const import AntennaElementImtConst
from sharc.antenna.antenna_subarray_imt import AntennaSubarrayIMT
from sharc.antenna.antenna import Antenna
from sharc.parameters.imt.parameters_antenna_imt import ParametersAntennaImt

class AntennaBeamformingImt(Antenna):
    def __init__(self, par: ParametersAntennaImt, azimuth: float, elevation: float):
        super().__init__()
        self.param = par
        self.subarray = None

        # Seleção do Elemento
        if (par.element_pattern).upper() == "M2101":
            self.element = AntennaElementImtM2101(par)
        elif (par.element_pattern).upper() == "F1336":
            self.element = AntennaElementImtF1336(par)
        elif (par.element_pattern).upper() == "FIXED":
            self.element = AntennaElementImtConst(par)

        # Configuração do Subarray
        if par.subarray.is_enabled:
            self.subarray = AntennaSubarrayIMT(
                element=self.element,
                eletrical_downtilt=par.subarray.eletrical_downtilt,
                n_rows=par.subarray.n_rows,
                element_vert_spacing=par.subarray.element_vert_spacing
            )

        self.azimuth = azimuth
        self.elevation = elevation
        self._calculate_rotation_matrix()
        self.minimum_array_gain = par.minimum_array_gain
        self.n_rows = par.n_rows
        self.n_cols = par.n_columns
        self.dh = par.element_horiz_spacing
        self.dv = par.element_vert_spacing
        
        self.beams_list = []
        self.w_vec_list = []

    def add_beam(self, phi_etilt: float, theta_etilt: float):
        # Converte para local e gera pesos (W)
        phi, theta = self.to_local_coord(phi_etilt, theta_etilt)
        self.beams_list.append((phi.item(), theta.item() - 90))
        self.w_vec_list.append(self._weight_vector(phi, theta - 90))

    def calculate_gain(self, phi_vec, theta_vec, beam_idx=0) -> np.array:
        """
        Cálculo de Ganho VETORIZADO (Igual ao MATLAB)
        """
        # 1. Coordenadas Locais
        lo_phi, lo_theta = self.to_local_coord(phi_vec, theta_vec)
        
        # 2. Ganho do Elemento ou Subarray (Vetorizado)
        if self.subarray is None:
            g_elem = self.element.element_pattern(lo_phi, lo_theta)
        else:
            g_elem = self.subarray.calculate_gain(lo_phi, lo_theta)
        
        # 3. Fator de Arranjo (AF) Vetorizado
        v_vec = self._super_position_vector_vectorized(lo_phi, lo_theta)
        w_vec = self.w_vec_list[beam_idx]
        
        af_linear = np.abs(np.sum(v_vec * w_vec, axis=(1, 2)))**2
        g_arr = 10 * np.log10(af_linear + 1e-12)
        
        # 4. Ganho Total
        g_total = g_elem + g_arr
        return np.maximum(g_total, self.minimum_array_gain)

    def _super_position_vector_vectorized(self, phi, theta):
        r_phi, r_theta = np.deg2rad(phi), np.deg2rad(theta)
        n = (np.arange(self.n_rows) + 1).reshape(1, -1, 1)
        m = (np.arange(self.n_cols) + 1).reshape(1, 1, -1)
        
        exp_arg = (n - 1) * self.dv * np.cos(r_theta[:, None, None]) + \
                  (m - 1) * self.dh * np.sin(r_theta[:, None, None]) * np.sin(r_phi[:, None, None])
        return np.exp(2j * np.pi * exp_arg)

    def _weight_vector(self, phi_tilt, theta_tilt):
        r_phi, r_theta = np.deg2rad(phi_tilt), np.deg2rad(theta_tilt)
        n = (np.arange(self.n_rows) + 1).reshape(-1, 1)
        m = (np.arange(self.n_cols) + 1).reshape(1, -1)
        exp_arg = (n - 1) * self.dv * np.sin(r_theta) - (m - 1) * self.dh * np.cos(r_theta) * np.sin(r_phi)
        return (1 / np.sqrt(self.n_rows * self.n_cols)) * np.exp(2j * np.pi * exp_arg)

    def to_local_coord(self, phi, theta):
        phi_rad, theta_rad = np.deg2rad(np.atleast_1d(phi)), np.deg2rad(np.atleast_1d(theta))
        pts = np.array([np.sin(theta_rad)*np.cos(phi_rad), np.sin(theta_rad)*np.sin(phi_rad), np.cos(theta_rad)])
        rot = np.asarray(self.rotation_mtx @ pts)
        return np.rad2deg(np.arctan2(rot[1], rot[0])), np.rad2deg(np.arccos(np.clip(rot[2], -1, 1)))

    def _calculate_rotation_matrix(self):
        a, b = np.deg2rad(self.azimuth), np.deg2rad(self.elevation)
        ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
        rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
        self.rotation_mtx = ry @ rz.T

    def reset_beams(self):
        self.beams_list, self.w_vec_list = [], []


class PlotAntennaPattern:
    def __init__(self, figs_dir): self.figs_dir = figs_dir

    def plot_element_pattern(self, antenna, sta_type, plot_type):
        theta_obs = np.linspace(0, 180, 720)
        
        if plot_type == "PAIR_ARRAY":
            theta_tilts = theta_tilts = [
                92.58, 92.60, 92.63]
            phi_tilts = np.linspace(-60, 60, 48)
            phi_scans = np.linspace(-180, 172, 45)

            sum_linear = np.zeros_like(theta_obs)
            total_iterations = len(theta_tilts) * len(phi_tilts) * len(phi_scans)

            for t_tilt in theta_tilts:
                for p_tilt in phi_tilts:
                    antenna.reset_beams()
                    antenna.add_beam(p_tilt, t_tilt)
                    
                    for azim in phi_scans:
                        phi_vec = azim * np.ones_like(theta_obs)
                        gains = antenna.calculate_gain(phi_vec, theta_obs, beam_idx=0)
                        sum_linear += 10**(gains / 10.0)

            mean_gain_db = 10 * np.log10(sum_linear / total_iterations)
            g_plot = mean_gain_db - np.max(mean_gain_db)

            pd.DataFrame({'elevation_axis': 90-theta_obs, 'g_plot': g_plot}).to_csv('resultado_sincronizado.csv', index=False)

            plt.figure(figsize=(10, 6))
            plt.plot(90-theta_obs, g_plot, 'b-', lw=2)
            plt.grid(True)
            plt.xlabel("Elevation [deg]")
            plt.ylabel("Gain [dB]")
            plt.title("Média Rigorosa - Python Sincronizado com MATLAB")
            plt.xlim(-90, 90)
            plt.ylim(-60, 5)
            plt.show()

#  ==========================================================
# MAIN (INALTERADO)
# ==========================================================
if __name__ == "__main__":
    figs_dir = "figs/"

    bs_param = ParametersAntennaImt()
    bs_param.adjacent_antenna_model = "BEAMFORMING"
    bs_param.element_pattern = "M2101"
    bs_param.element_max_g = 6.4
    bs_param.element_phi_3db = 90
    bs_param.element_theta_3db = 65
    bs_param.element_am = 30
    bs_param.element_sla_v = 30
    bs_param.n_rows = 8
    bs_param.n_columns = 16
    bs_param.element_horiz_spacing = 0.5
    bs_param.element_vert_spacing = 2.1
    bs_param.multiplication_factor = 12
    bs_param.downtilt = 6

    bs_param.subarray.is_enabled = True
    bs_param.subarray.n_rows = 3
    bs_param.subarray.element_vert_spacing = 0.7
    bs_param.subarray.eletrical_downtilt = 3

    par = bs_param.get_antenna_parameters()
    bs_array = AntennaBeamformingImt(par, 0, -bs_param.downtilt)

    plot = PlotAntennaPattern(figs_dir)
    plot.plot_element_pattern(bs_array, "TX", "PAIR_ARRAY")