import numpy as np
import matplotlib.pyplot as plt
import csv

versions = ['1743isd_omni','1743isd_itursa509_azimute0', '1743isd_itursa509_azimute45', '1743isd_itursa509_azimute-90']  # Lista de versões

# Função para carregar e processar os dados do arquivo CSV
def carregar_dados_csv(file_path, skip_rows=2):
    interferencia_dBW = []
    prob_interferencia = []
    
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        # Ignorar as primeiras linhas, se necessário
        for _ in range(skip_rows):
            next(reader)
        
        # Ler os dados restantes
        for row in reader:
            if len(row) >= 2:  # Garantir que há pelo menos duas colunas
                interferencia_dBW.append(float(row[0]) - 30)  # Subtrair 30 para converter para dBW
                prob_interferencia.append(1 - float(row[1]))   # Calcular a CCDF (1 - CDF)
    
    return interferencia_dBW, prob_interferencia

# Plotar os dados de cada versão
for version in versions:
    # Caminho do arquivo CSV para cada versão
    file_path1 = f'/home/ju/Documents/Github/sharc/sharc/campaigns/imt_dl_ras_10000_MHz/output/imt_dl_ras_{version}/SYS_CDF_of_system_interference_power_from_IMT_DL.csv'
    
    # Carregar os dados
    interferencia_dBW_1, prob_interferencia_1 = carregar_dados_csv(file_path1)
    
    # Criar o gráfico para cada versão
    plt.plot(interferencia_dBW_1, prob_interferencia_1, label=f"IMT Downlink - {version}")

# Adicionar os dados da contribuição 22 (RAS_DL_50Km.csv)
file_path_contrib_22 = '/home/ju/Downloads/RAS_DL_50Km.csv'

# Carregar e processar os dados da contribuição 22
interferencia_dBW_contrib_22, prob_interferencia_contrib_22 = carregar_dados_csv(file_path_contrib_22)

# Aplique 1 - probabilidade para cada ponto de probabilidade na contribuição 22 (CCDF)
prob_interferencia_contrib_22 = [1 - p for p in prob_interferencia_contrib_22]
interferencia_dBW_contrib_22 = [i + 30 for i in interferencia_dBW_contrib_22]

# Plotar os dados da contribuição 22 (aplicando 1 - probabilidade para CCDF)
plt.plot(interferencia_dBW_contrib_22, prob_interferencia_contrib_22, label="Contribuição 22 - RAS DL 50 km")

# Adicionar linha vertical para o critério de proteção
plt.axvline(x=-202, color='black', linestyle='--', linewidth=1, label='Critério de Proteção (-202 dBW)')

# Configurar título e eixos
plt.title("Comparação de Interferência IMT Downlink (CCDF)")
plt.xlabel("Interferência (dBW)")
plt.ylabel("Probabilidade de Interferência > X")
plt.yscale('log')

# Adicionar grid ao gráfico
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Exibir a legenda
plt.legend(loc='lower left')

# Mostrar o gráfico
plt.show()
