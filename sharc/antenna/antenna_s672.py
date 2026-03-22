# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 17:18:59 2017

@author: edgar
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_single_space_station import ParametersSingleSpaceStation


import numpy as np
import sys


class AntennaS672(Antenna):
    """
    Implements the satellite antenna pattern in the fixed-satellite service
    according to Recommendation ITU-R S.672-4 Annex 1
    """

    def __init__(self, param: ParametersSingleSpaceStation):
        super().__init__()
        self.peak_gain = param.antenna_gain
        self.l_s = param.antenna_l_s

        #Cutting opening to be deleted 
        self.cutting_angle = getattr(param, "antenna_cutting_angle", None) #If the "antenna_cutting_angle" attribute has not been defined


        if self.l_s == -20:
            self.a = 2.58
        elif self.l_s == -25:
            self.a = 2.88
        elif self.l_s == -30:
            self.a = 3.16
        else:
            sys.stderr.write(
                "ERROR\nInvalid AntennaS672 L_s parameter: " + self.l_s,
            )
            sys.exit(1)

        self.b = 6.32

        self.psi_0 = param.antenna_3_dB / 2
        self.psi_1 = self.psi_0 * \
            np.power(10, (self.peak_gain + self.l_s + 20) / 25)

    def calculate_gain(self, *args, **kwargs) -> np.array:
        psi = np.absolute(kwargs["off_axis_angle_vec"])
        gain = np.zeros(len(psi))

        idx_0 = np.where(psi < self.psi_0)
        gain[idx_0] = self.peak_gain            

        idx_1 = np.where((self.psi_0 <= psi) & (psi <= self.a * self.psi_0))[0]
        gain[idx_1] = self.peak_gain - 3 * np.power(psi[idx_1] / self.psi_0, 2)

        idx_2 = np.where((self.a * self.psi_0 < psi) &
                         (psi <= self.b * self.psi_0))[0]
        gain[idx_2] = self.peak_gain + self.l_s

        idx_3 = np.where((self.b * self.psi_0 < psi) & (psi <= self.psi_1))[0]
        gain[idx_3] = self.peak_gain + self.l_s + \
            20 - 25 * np.log10(psi[idx_3] / self.psi_0)
        
        #Applying crop factor
        if self.cutting_angle : 
            idx_cut = np.where((psi > self.cutting_angle))[0]
            gain[idx_cut] = -500 #very low arbitrary gain

        return gain


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    # initialize antenna parameters Cenario 1
    param = ParametersSingleSpaceStation()
    param.antenna_gain = 55
    param.antenna_pattern = "ITU-R S.672-4"
    param.antenna_3_dB = 0.29
    param.antenna_cutting_angle = 0.4
    psi = np.linspace(0.1, 90, num=10000)

    param.antenna_l_s = -20
    antenna = AntennaS672(param)
    gain = antenna.calculate_gain(off_axis_angle_vec=psi)

    # initialize antenna parameters Cenario 2
    param2 = ParametersSingleSpaceStation()
    param2.antenna_gain = 30
    param2.antenna_pattern = "ITU-R S.672-4"
    param2.antenna_3_dB = 5.0

    param2.antenna_l_s = -20
    antenna2 = AntennaS672(param2)
    gain_2 = antenna2.calculate_gain(off_axis_angle_vec=psi)

    # initialize antenna parameters Cenario 3
    param3 = ParametersSingleSpaceStation()
    param3.antenna_gain = 38
    param3.antenna_pattern = "ITU-R S.672-4"
    param3.antenna_3_dB = 2.2

    param3.antenna_l_s = -20
    antenna_3 = AntennaS672(param3)
    gain_3 = antenna_3.calculate_gain(off_axis_angle_vec=psi)
    """
    # initialize antenna parameters (17.2)
    param2 = ParametersSingleSpaceStation()
    param2.antenna_gain = 35
    param2.antenna_pattern = "ITU-R S.672-4"
    param2.antenna_3_dB = 17.2
    
    param2.antenna_l_s = -20
    antenna = AntennaS672(param2)
    gain35 = antenna.calculate_gain(off_axis_angle_vec=psi)
    """
   
    fig = plt.figure(
        figsize=(12, 7), facecolor='w',
        edgecolor='k',
    )  # create a figure object
   
    plt.semilogx(
        psi, gain ,
        "-b", label=f"BW = {param.antenna_3_dB}º, gain = {param.antenna_gain} dBi",
    )
    
    """
    plt.semilogx(
        psi, gain_2 ,
        "-r", label=f"BW = {param2.antenna_3_dB }º, gain = {param2.antenna_gain} dBi",
    ) 
    plt.semilogx(
        psi, gain_3 ,
        "-g", label=f"BW = {param3.antenna_3_dB }º, gain = {param3.antenna_gain} dBi",
    ) 
    """

    plt.ylim((-0.5, 56))
    plt.xlim((0.1, 30))
    plt.title("ITU-R S.672-4 antenna radiation pattern ($L_S = -20$ dB)")
    plt.xlabel(r"off-axis angle [°]")
    plt.ylabel("Gain [dBi]")
    plt.legend(loc="upper right")

    #Ante3s do ponto
    
    ax = plt.gca()
    ax.set_yticks([0,10,20,30,40,50])

    # Potências de 10 (0.1, 1, 10, 100)
    ticks_potencias = np.logspace(-1, 2, 4)

    # Marcações intermediárias (0.2, 0.3, ..., 0.9 e 2, 3, ..., 9)
    ticks_intermediarios = np.concatenate([
        np.linspace(0.2, 0.9, 8),  # Entre 0.1 e 1
        np.linspace(2, 9, 8),      # Entre 1 e 10
        np.linspace(20, 90, 8)     # Entre 10 e 100
    ])

    # Unir os dois conjuntos de ticks
    xticks = np.concatenate([ticks_potencias, ticks_intermediarios])

    # Definir os xticks no gráfico
    ax.set_xticks(xticks)

    #ax.set_xticks(np.logspace(-1, 2, 4))  # De 10^-1 até 10^2
    """

    #Ponto
    ax = plt.gca()
    ax.set_yticks([0, 10, 20, 30, 40, 50])

    # Exemplo de ponto
    ponto = float(psi[gain == 35][0])

    #fig, ax = plt.subplots()
    ax.set_yticks([0, 10, 20, 30, 40, 50])

    # Linha vertical em ponto
    plt.vlines(ponto, 0, 60, color="red")

    # Potências de 10
    ticks_potencias = np.logspace(-1, 2, 4)  # [0.1, 1, 10, 100]

    # Intermediários (sem label)
    ticks_intermediarios = np.concatenate([
        np.linspace(0.2, 0.9, 8),
        np.linspace(2, 9, 8),
        np.linspace(20, 90, 8)
    ])

    # Todos os xticks
    xticks = np.concatenate([ticks_potencias, ticks_intermediarios, [ponto]])
    xticks = np.unique(np.round(xticks, 6))
    ax.set_xticks(xticks)

    # Criar labels: só para potências de 10 e ponto
    labels = []
    for x in xticks:
        if np.isclose(x, ponto, atol=1e-6):
            labels.append(f'{ponto:.4f}')
        elif x in ticks_potencias:
            labels.append(f'{x:g}')
        else:
            labels.append('')  # Sem label para intermediários

    # Aplicar os labels
    ax.set_xticklabels(labels)

    # Mudar cor apenas da label do ponto
    for label, x in zip(ax.get_xticklabels(), xticks):
        if np.isclose(x, ponto, atol=1e-6):
            label.set_color('red')

    plt.tight_layout()
    plt.show()

    """
    plt.grid()
    plt.show()
