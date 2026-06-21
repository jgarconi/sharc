import numpy as np
import random
from sharc.parameters.imt.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_antenna import ParametersAntenna
from antenna_beamforming_imt import AntennaBeamformingImt 
from antenna_f1245_fs import Atenna_f1245_fs

# Potência de Transmissão da BS IMT [dBm]
POT_BS_IMT = 30

# Configuração do IMT
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
bs_param.normalization = False
bs_param.minimum_array_gain = -200

bs_param.subarray.is_enabled = True
bs_param.subarray.n_rows = 3
bs_param.subarray.element_vert_spacing = 0.7
bs_param.subarray.eletrical_downtilt = 3

# Configuração do Setor 1
bs_azimuth_fisico = 0.0 

# Inicializa a antena. Note que o elevation aqui é o tilt físico (x-axis ref)
bs_antenna = AntennaBeamformingImt(bs_param.get_antenna_parameters(), bs_azimuth_fisico, -bs_param.downtilt)

# Sorteio do UE dentro dos 120 graus deste setor específico
ue_azimuth = (bs_azimuth_fisico + random.uniform(-60, 60)) % 360
ue_elevation = random.uniform(90, 100)

# Configuração FS 
fs_par = ParametersAntenna()
fs_par.gain = 36.0
fs_par.frequency = 8000
fs_par.diameter = 2.0
fs_antenna = Atenna_f1245_fs(fs_par)
fs_antenna.add_beam(0, 0) 

# Configuração da Topologia
casos = [
    {"caso": "EP01-Site10", "AzEP": -5.9, "ElevEP":  0.8, "AzSite": 131.0 , "ElevSite": -0.8, "PathLoss": 148.8},
    {"caso": "EP01-Site15", "AzEP": 22.6, "ElevEP": -1.5, "AzSite": 102.5 , "ElevSite":  1.5, "PathLoss": 153.3},
    {"caso": "EP02-Site05", "AzEP": 34.3, "ElevEP": -0.9, "AzSite": 90.8  , "ElevSite":  0.9, "PathLoss": 172.2},
    {"caso": "EP02-Site15", "AzEP": 41.7, "ElevEP": -0.9, "AzSite": 83.4  , "ElevSite":  0.9, "PathLoss": 172.4}

]

print(f"UE posicionado em: Az={ue_azimuth:.2f}°, Elev={ue_elevation:.2f}°\n")

for c in casos:
    # BS IMT cria o feixe eletrônico na direção do UE
    bs_antenna.reset_beams()
    bs_antenna.add_beam(ue_azimuth, ue_elevation) # Adiciona ao beams_list

    target_az = c["AzEP"] % 360
    target_el = 90 - c["ElevEP"]

    # CALCULO IMT: Correção do nome do argumento para 'beams_l' e uso de arrays
    imt_gain = bs_antenna.calculate_gain(
        phi_vec=np.array([target_az]),
        theta_vec=np.array([target_el]), 
        beams_l=np.array([0]), # Índice do feixe apontado para o UE
        co_channel=True
    )[0]

    # CALCULO FS:
    off_axis = fs_antenna.calculate_off_axis_angle(
        Az=c["AzSite"],
        b=c["ElevSite"]
    )
    fs_gain = fs_antenna.calculate_gain(
        off_axis_angle_vec=np.array([off_axis])
    )[0]

    # Interferência Final
    interferencia = POT_BS_IMT + imt_gain + fs_gain - c["PathLoss"]
    
    print(f"Caso: {c['caso']}")
    print(f"Ganho IMT (vazamento): {imt_gain:.2f} dBi")
    print(f"Ganho FS (recepção): {fs_gain:.2f} dBi")
    print(f"Interferência: {interferencia:.2f} dBm")
    print("-" * 35)