import numpy as np
import matplotlib.pyplot as plt

def cdf_empirical(X):
    """
    Compute the Empirical Cumulative Distribution Function (CDF) of the input array X.
    The empirical CDF y=F(x) is defined as the proportion of X values less than or equal to x.
    If input X is a matrix, it flattens the array and computes the CDF of all values.
    """
    Xplot = np.sort(X.flatten())
    Yplot = np.arange(1, len(Xplot) + 1) / len(Xplot)

    # Salvar Xplot em um arquivo de texto
    #np.savetxt('/home/ju/Downloads/Xplot_output_python.txt', Xplot, delimiter='\t', fmt='%.10f')

    return np.column_stack((Xplot, Yplot))

def process_inr_data(INR_DL_LG, INR_UL_LG, INR_DL_SM, INR_UL_SM):
    """
    Process INR data, convert to linear scale, and calculate aggregate values.
    """
    # Adjust INR data
    INR_DL_LG = INR_DL_LG[:, 1] + 10 * np.log10(2047.6)
    INR_UL_LG = INR_UL_LG[:, 1] + 10 * np.log10(2047.6)
    INR_DL_SM = INR_DL_SM[:, 1] + 10 * np.log10(1.14)
    INR_UL_SM = INR_UL_SM[:, 1] + 10 * np.log10(1.14)

    # Convert to linear values
    INR_DL_SM_lin = 10 ** (INR_DL_SM / 10)
    INR_DL_LG_lin = 10 ** (INR_DL_LG / 10)
    INR_UL_SM_lin = 10 ** (INR_UL_SM / 10)
    INR_UL_LG_lin = 10 ** (INR_UL_LG / 10)

    # Calculate aggregated INR values (DL and UL)
    INR_DL = INR_DL_SM_lin + INR_DL_LG_lin
    INR_UL = INR_UL_SM_lin + INR_UL_LG_lin

    return INR_DL_SM_lin, INR_DL_LG_lin, INR_UL_SM_lin, INR_UL_LG_lin, INR_DL, INR_UL

def calculate_aggregate_inr(INR_DL, INR_UL, INR_DL_SM_lin, INR_DL_LG_lin, INR_UL_SM_lin, INR_UL_LG_lin):
    """
    Generate random indices and compute the aggregated INR values.
    """
    # Generate random indices for aggregation
    indices = [np.random.randint(0, len(INR_DL), len(INR_DL)) for _ in range(6)]

    # Compute aggregated INR values
    INR_Agg_1 = 10 * np.log10(0.75 * (INR_DL_SM_lin[indices[2]] + INR_DL_LG_lin[indices[3]]) + 
                              0.25 * (INR_UL_SM_lin[indices[4]] + INR_UL_LG_lin[indices[5]]))
    INR_Agg_2 = 10 * np.log10(0.75 * INR_DL[indices[0]] + 0.25 * INR_UL[indices[1]])

    return INR_Agg_1, INR_Agg_2

def plot_cdfs(data_dict):
    """
    Plot multiple CDFs from a dictionary of data.
    """
    plt.figure(figsize=(12, 8))
    for label, (X, Y) in data_dict.items():
        plt.plot(X, Y, linewidth=2, label=label)

    plt.legend()
    plt.xlabel('Interference to Noise Ratio [dB]')
    plt.ylabel('CDF')
    plt.title('CDF of Interference to Noise Ratio for 18deg Nadir - EESS Active')
    plt.grid(True)
    plt.show()

version = "2024-09-17_01"

# Leitura dos arquivos CSV para as variáveis INR_DL_LG, INR_UL_LG, INR_DL_SM, INR_UL_SM
INR_DL_LG = np.loadtxt(f'/home/ju/Documents/Github/sharc/sharc/campaigns/imt_hotspot_eess_active/output/imt_dl_hotspot_eess_active_beam_large_{version}/SYS_INR_samples.csv', delimiter=',', skiprows=2)
INR_UL_LG = np.loadtxt(f'/home/ju/Documents/Github/sharc/sharc/campaigns/imt_hotspot_eess_active/output/imt_ul_hotspot_eess_active_beam_large_{version}/SYS_INR_samples.csv', delimiter=',', skiprows=2)
INR_DL_SM = np.loadtxt(f'/home/ju/Documents/Github/sharc/sharc/campaigns/imt_hotspot_eess_active/output/imt_dl_hotspot_eess_active_beam_small_{version}/SYS_INR_samples.csv', delimiter=',', skiprows=2)
INR_UL_SM = np.loadtxt(f'/home/ju/Documents/Github/sharc/sharc/campaigns/imt_hotspot_eess_active/output/imt_ul_hotspot_eess_active_beam_small_{version}/SYS_INR_samples.csv', delimiter=',', skiprows=2)

# Carregar o valor de referência
EESS_REFERENCE = np.loadtxt('/home/ju/Downloads/EESS_ACTIVE_18.csv', delimiter=',', skiprows=1)

# Process INR data
INR_DL_SM_lin, INR_DL_LG_lin, INR_UL_SM_lin, INR_UL_LG_lin, INR_DL, INR_UL = process_inr_data(INR_DL_LG, INR_UL_LG, INR_DL_SM, INR_UL_SM)

# Calculate aggregate INR values
INR_Agg_1, INR_Agg_2 = calculate_aggregate_inr(INR_DL, INR_UL, INR_DL_SM_lin, INR_DL_LG_lin, INR_UL_SM_lin, INR_UL_LG_lin)

# Calculate empirical CDFs
cdf_data = {
#    "DL - large beam": cdf_empirical(INR_DL_LG),
#    "UL - large beam": cdf_empirical(INR_UL_LG),
#    "DL - 3dB spotbeam": cdf_empirical(INR_DL_SM),
#    "UL - 3dB spotbeam": cdf_empirical(INR_UL_SM),
    "Aggregate CDF 1": cdf_empirical(INR_Agg_1),
    "Aggregate CDF 2": cdf_empirical(INR_Agg_2),
    "EESS Reference": (EESS_REFERENCE[:, 0], EESS_REFERENCE[:, 1])
}

# Prepare data for plotting
cdf_dict = {}
for label, cdf in cdf_data.items():
    if isinstance(cdf, tuple):
        # Se já for uma tupla, não precisa fazer nada
        cdf_dict[label] = cdf
    else:
        # Caso contrário, acesse as colunas 0 e 1
        cdf_dict[label] = (cdf[:, 0], cdf[:, 1])

# Plot all CDFs
plot_cdfs(cdf_dict)
