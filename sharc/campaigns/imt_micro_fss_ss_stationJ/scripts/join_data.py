import os
from pathlib import Path
import pandas as pd
from collections import defaultdict

# Caminho base da campanha
campaign_base_dir = str((Path(__file__) / "../.." ).resolve())
dl_dir = os.path.join(campaign_base_dir, "output_dl")

# Novo diretório para salvar arquivos concatenados
output_dir = os.path.join(dl_dir, "output_imt_micro_fss_ss_stationJ_small_dl_2025-05-31_01")
os.makedirs(output_dir, exist_ok=True)

# Dicionário para agrupar arquivos por nome
arquivos_por_nome = defaultdict(list)

# Apenas as subpastas de primeiro nível em dl_dir (exceto a de saída)
for subdir in os.listdir(dl_dir):
    subdir_path = os.path.join(dl_dir, subdir)

    if not os.path.isdir(subdir_path) or subdir_path == output_dir:
        continue  # Ignora se não for pasta ou for a pasta de saída

    for file in os.listdir(subdir_path):
        if file.endswith(".csv"):
            caminho_arquivo = os.path.join(subdir_path, file)
            arquivos_por_nome[file].append(caminho_arquivo)

# Para cada grupo de arquivos com o mesmo nome
for nome_arquivo, lista_caminhos in arquivos_por_nome.items():
    list_df = []
    for caminho in lista_caminhos:
        try:
            df = pd.read_csv(caminho)
            list_df.append(df)
        except Exception as e:
            print(f"Erro ao ler {caminho}: {e}")

    if list_df:
        df_concat = pd.concat(list_df, ignore_index=True)
        output_path = os.path.join(output_dir, nome_arquivo)
        df_concat.to_csv(output_path, index=False)
        print(f"Salvo: {output_path}")
