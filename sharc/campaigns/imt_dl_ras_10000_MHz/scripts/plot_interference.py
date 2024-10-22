import numpy as np
import matplotlib.pyplot as plt
import csv

# Caminhos dos arquivos CSV
file_path1 = '/home/ju/Documents/Github/sharc/sharc/campaigns/imt_dl_ras_10000_MHz/output/imt_dl_ras_2024-10-20_01/SYS_CDF_of_system_interference_power_from_IMT_DL.csv'
file_path2 = '/home/ju/Downloads/contribuicao22.csv'

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

# Carregar os dados dos dois arquivos
interferencia_dBW_1, prob_interferencia_1 = carregar_dados_csv(file_path1)
interferencia_dBW_2, prob_interferencia_2 = carregar_dados_csv(file_path2)

# Criar o gráfico
plt.figure(figsize=(10,5))
plt.plot(interferencia_dBW_1, prob_interferencia_1, label="IMT Downlink - Calculado")
#plt.plot(interferencia_dBW_2, prob_interferencia_2, label="Resultado Contribuição 22")
plt.axvline(x=-202, color='red', linestyle='--', linewidth=1, label='Critério de Proteção (-202 dBW)')  # Adicionar linha vertical

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
