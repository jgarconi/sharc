import numpy as np
import matplotlib.pyplot as plt

version = "2024-09-16_02"

# Carregar os dados já calculados dos arquivos
INR_DL_LG = np.loadtxt(f'/home/ju/Documents/Github/sharc/sharc/campaigns/imt_hotspot_eess_active/output/imt_dl_hotspot_eess_active_beam_large_{version}/SYS_CDF_of_system_INR.csv', delimiter=',', skiprows=2)
INR_UL_LG = np.loadtxt(f'/home/ju/Documents/Github/sharc/sharc/campaigns/imt_hotspot_eess_active/output/imt_ul_hotspot_eess_active_beam_large_{version}/SYS_CDF_of_system_INR.csv', delimiter=',', skiprows=2)
INR_DL_SM = np.loadtxt(f'/home/ju/Documents/Github/sharc/sharc/campaigns/imt_hotspot_eess_active/output/imt_dl_hotspot_eess_active_beam_small_{version}/SYS_CDF_of_system_INR.csv', delimiter=',', skiprows=2)
INR_UL_SM = np.loadtxt(f'/home/ju/Documents/Github/sharc/sharc/campaigns/imt_hotspot_eess_active/output/imt_ul_hotspot_eess_active_beam_small_{version}/SYS_CDF_of_system_INR.csv', delimiter=',', skiprows=2)

# Carregar o valor de referência
EESS_REFERENCE = np.loadtxt('/home/ju/Downloads/EESS_ACTIVE_18.csv', delimiter=',', skiprows=1)

# Converter os dados de INR para escala linear
INR_DL_SM_lin = 10 ** (INR_DL_SM[:, 0] / 10)
INR_DL_LG_lin = 10 ** (INR_DL_LG[:, 0] / 10)
INR_UL_SM_lin = 10 ** (INR_UL_SM[:, 0] / 10)
INR_UL_LG_lin = 10 ** (INR_UL_LG[:, 0] / 10)

# Calcular o INR agregado
INR_DL = INR_DL_SM_lin + INR_DL_LG_lin
INR_UL = INR_UL_SM_lin + INR_UL_LG_lin

# Gerar índices aleatórios para o agregado
indices = [np.random.randint(0, len(INR_DL), len(INR_DL)) for _ in range(6)]

# Calcular o INR agregado
INR_Agg_1 = 10 * np.log10(0.75 * (INR_DL_SM_lin[indices[2]] + INR_DL_LG_lin[indices[3]]) + 
                          0.25 * (INR_UL_SM_lin[indices[4]] + INR_UL_LG_lin[indices[5]]))
INR_Agg_2 = 10 * np.log10(0.75 * INR_DL[indices[0]] + 0.25 * INR_UL[indices[1]])

# Definir uma função para plotar as CDFs
def plot_cdfs(data_dict):
    """
    Plot multiple CDFs from a dictionary of data.
    """
    plt.figure()
    for label, (X, Y) in data_dict.items():
        plt.plot(X, Y, linewidth=2, label=label)

    plt.legend()
    plt.xlabel('Interference to Noise Ratio [dB]')
    plt.ylabel('CDF [%]')
    plt.title('CDF of Interference to Noise Ratio for 18deg Nadir - EESS Active')
    plt.grid(True)
    plt.show()

# Preparar os dados de CDF para o gráfico, incluindo a referência
cdf_data = {
    #"DL - large beam": (INR_DL_LG[:, 0], INR_DL_LG[:, 1]),
    #"UL - large beam": (INR_UL_LG[:, 0], INR_UL_LG[:, 1]),
    "DL - 3dB spotbeam": (INR_DL_SM[:, 0], INR_DL_SM[:, 1]),
    "UL - 3dB spotbeam": (INR_UL_SM[:, 0], INR_UL_SM[:, 1]),
    "Aggregate CDF 1": (np.sort(INR_Agg_1), np.arange(1, len(INR_Agg_1) + 1) / len(INR_Agg_1)),
    #"Aggregate CDF 2": (np.sort(INR_Agg_2), np.arange(1, len(INR_Agg_2) + 1) / len(INR_Agg_2)),
    "EESS Reference": (EESS_REFERENCE[:, 0], EESS_REFERENCE[:, 1])
}

# Plotar todas as CDFs
plot_cdfs(cdf_data)
