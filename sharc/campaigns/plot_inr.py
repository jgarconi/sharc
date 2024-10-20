import numpy as np
import matplotlib.pyplot as plt

# Carregar os dados dos dois arquivos
# Carregar os dados ignorando a primeira linha (cabeçalho)
dados1 = np.loadtxt('/home/ju/Downloads/EESS_ACTIVE_18.csv', delimiter=',', skiprows=1)
dados2 = np.loadtxt('/home/ju/Documents/Github/sharc/sharc/campaigns/imt_hotspot_eess_active/output/imt_dl_hotspot_eess_active_beam_small_2024-09-16_01/SYS_CDF_of_system_INR.csv', delimiter=',', skiprows=2)

# Extrair as colunas necessárias
INR1 = dados1[:, 0]
prob_interferencia1 = dados1[:, 1]

INR2 = dados2[:, 0]
prob_interferencia2 = dados2[:, 1]

# Criar o gráfico
plt.figure(figsize=(10, 5))

# Plotar a primeira curva
plt.plot(INR1, prob_interferencia1, label="EESS-Luciano: 50°")

# Plotar a segunda curva
plt.plot(INR2, prob_interferencia2, label="EESS-Segundo: 18°")
# Adicionar o título e os eixos
plt.title("Análise de sensibilidade")
plt.xlabel("INR (dB)")
plt.ylabel("Probabilidade de INR < X")

# Adicionar grid ao gráfico
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Exibir a legenda
plt.legend(loc='center left', bbox_to_anchor=(0.8, 1.05))

# Mostrar o gráfico
plt.show()
