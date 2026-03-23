import os
import shutil
from pathlib import Path
import pandas as pd
from collections import defaultdict

# Caminho base da campanha
campaign_base_dir = "/home/juliana/Documentos/projetos/sharc/sharc/campaigns/resim_study_A.1_FSS/output"

# Novos diretórios de saída
output_final = os.path.join(campaign_base_dir, "output_resim_study_A.1_FSS")

os.makedirs(output_final, exist_ok=True)

def processar_base(base_dir, output_dir, prefixo):
    arquivos_por_nome = defaultdict(list)
    yaml_files = []

    # Apenas as subpastas de primeiro nível com o prefixo desejado
    for subdir in os.listdir(base_dir):
        if not subdir.startswith(prefixo):
            continue  # pula pastas que não têm o prefixo correto

        subdir_path = os.path.join(base_dir, subdir)
        if not os.path.isdir(subdir_path) or subdir_path == output_dir:
            continue

        for file in os.listdir(subdir_path):
            caminho_arquivo = os.path.join(subdir_path, file)
            if file.endswith(".csv"):
                arquivos_por_nome[file].append(caminho_arquivo)
            elif file.endswith(".yaml"):
                yaml_files.append(caminho_arquivo)

    # Concatena CSVs
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

    # Copia arquivos YAML
    for yaml_file in yaml_files:
        try:
            shutil.copy(yaml_file, output_dir)
            print(f"Copiado: {yaml_file}")
        except Exception as e:
            print(f"Erro ao copiar {yaml_file}: {e}")

# Processa DL e UL separadamente com filtro por prefixo
processar_base(campaign_base_dir, output_final, prefixo="out_resim")