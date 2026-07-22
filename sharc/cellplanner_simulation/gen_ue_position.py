# -*- coding: utf-8 -*-
"""
Script para gerar e visualizar:
 - 57 estações base IMT (TopologyMacrocell, num_clusters=1 -> 19 sites x 3 setores)
 - 171 UEs (3 por BS)
 - 20 estações receptoras FS ao redor da rede IMT

Rode isto DENTRO do seu ambiente com o pacote sharc instalado.
Ajuste os valores marcados com # <<< AJUSTE >>> conforme os parâmetros
reais do seu estudo de compartilhamento.
"""
import numpy as np
import matplotlib.pyplot as plt

from sharc.topology.topology_macrocell import TopologyMacrocell
from sharc.parameters.imt.parameters_imt import ParametersImt
from sharc.parameters.parameters_single_earth_station import ParametersSingleEarthStation
from sharc.station_factory import StationFactory
from sharc.support.enumerations import StationType

rand_gen = np.random.RandomState(101)

# ------------------------------------------------------------------
# 1) Topologia IMT: 19 sites x 3 setores = 57 BS
# ------------------------------------------------------------------
intersite_distance = 1500  # <<< AJUSTE >>> [m]
num_clusters = 1           # 1 -> 57 BS (19 sites x 3 setores)

topology = TopologyMacrocell(intersite_distance, num_clusters)
topology.calculate_coordinates(rand_gen)
print("Num BS:", topology.num_base_stations)

# ------------------------------------------------------------------
# 2) Parâmetros IMT mínimos para gerar 3 UEs por BS (57*3 = 171 UEs)
#    Ajuste os campos abaixo para bater com os parâmetros do seu
#    cenário (arquivo de config .yaml que você normalmente usa).
# ------------------------------------------------------------------
parameters = ParametersImt()
parameters.topology.type = "MACROCELL"
parameters.topology.macrocell.intersite_distance = intersite_distance

parameters.ue.k = 3                 # <<< AJUSTE >>> UEs ativos por BS
parameters.ue.k_m = 1
parameters.ue.azimuth_range = (-60, 60)
parameters.ue.distribution_type = "ANGLE_AND_DISTANCE"
parameters.ue.distribution_distance = "UNIFORM"       # ou "RAYLEIGH" / "SQRT(UNIFORM)"
parameters.ue.distribution_azimuth = "NORMAL"          # ou "UNIFORM"
parameters.ue.height = 1.5
parameters.ue.indoor_percent = 0
parameters.ue.noise_figure = 9                         # <<< AJUSTE >>>
parameters.bs.height = 6                                # <<< AJUSTE >>>
parameters.bandwidth = 100                               # <<< AJUSTE >>> [MHz]
parameters.frequency = 26000                             # <<< AJUSTE >>> [MHz]
parameters.minimum_separation_distance_bs_ue = 1

imt_ue = StationFactory.generate_imt_ue_outdoor(
    parameters,
    parameters.ue.antenna.array,
    rand_gen,
    topology,
)
print("Num UEs:", len(imt_ue.x))

# ------------------------------------------------------------------
# 3) 20 estações FS ao redor da rede
#    location.type = "UNIFORM_DIST" distribui a estação de forma
#    uniforme (em área) num anel entre min_dist_to_center e
#    max_dist_to_center, com origem no centro da topologia IMT.
# ------------------------------------------------------------------
network_span = np.sqrt(np.max(topology.x**2 + topology.y**2))
r_min_fs = network_span * 1.6   # <<< AJUSTE >>> distância mínima FS-rede
r_max_fs = network_span * 2.3   # <<< AJUSTE >>> distância máxima FS-rede

fs_x_list, fs_y_list = [], []
for i in range(20):
    fs_param = ParametersSingleEarthStation()
    fs_param.geometry.location.type = "UNIFORM_DIST"
    fs_param.geometry.location.uniform_dist.min_dist_to_center = r_min_fs
    fs_param.geometry.location.uniform_dist.max_dist_to_center = r_max_fs
    fs_param.geometry.height = 15               # <<< AJUSTE >>> altura antena FS [m]
    fs_param.geometry.azimuth.type = "POINTING_AT_IMT_CENTER"
    fs_param.geometry.elevation.type = "POINTING_AT_IMT_CENTER"
    fs_param.antenna_pattern = "OMNI"            # <<< AJUSTE >>> conforme padrão de antena real da FS
    fs_param.antenna_gain = 0                    # <<< AJUSTE >>>
    fs_param.bandwidth = 40                      # <<< AJUSTE >>> [MHz]
    fs_param.frequency = parameters.frequency
    fs_param.tx_power_density = -60              # não usado (é Rx), só placeholder
    fs_param.noise_temperature = 290             # <<< AJUSTE >>>
    fs_param.adjacent_ch_emissions = "OFF"

    fs_station = StationFactory.generate_single_earth_station(
        fs_param, rand_gen, StationType.SINGLE_EARTH_STATION, topology,
    )
    fs_x_list.append(fs_station.x[0])
    fs_y_list.append(fs_station.y[0])

fs_x = np.array(fs_x_list)
fs_y = np.array(fs_y_list)

# ------------------------------------------------------------------
# 4) Plot
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 9), facecolor='w')
topology.plot(ax)
ax.scatter(imt_ue.x, imt_ue.y, s=10, color='crimson', alpha=0.75,
           label=f"UE ({len(imt_ue.x)})", zorder=4)
ax.scatter(fs_x, fs_y, s=90, color='blue', marker='*',
           label=f"FS ({len(fs_x)})", zorder=6)

ax.set_aspect('equal')
ax.set_title(
    f"Topologia IMT ({topology.num_base_stations} BS / {len(imt_ue.x)} UEs) "
    f"+ {len(fs_x)} estações FS ao redor"
)
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.legend(loc='upper right')
plt.tight_layout()
plt.savefig("topologia_imt_fs.png", dpi=150)
plt.show()