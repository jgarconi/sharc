import pandas as pd
import os
from pathlib import Path
from sharc.results import Results, SampleList
from sharc.post_processor import PostProcessor
import time
import numpy as np

inicio = time.time()

# Define o diretório base da campanha
campaign_base_dir = str((Path(__file__) / ".." / "..").resolve())
dl_dir = os.path.join(campaign_base_dir, "output2")
# ul_dir = os.path.join(campaign_base_dir, "output_ul")

# Inicializa o pós-processador
post_processor = PostProcessor()

# Função para gerar legendas com base no nome do diretório
import re

def legend_gen(dir_name):
    link = re.search("_(dl|ul)", dir_name)
    if link is not None:
        link = link.group(1)
    else:
        return "None"
    
    return f"{link.upper()}"

post_processor.add_plot_legend_generator(legend_gen)

# Atributos a serem plotados
attributes_to_plot = [
    "system_inr",
    # "system_dl_interf_power_per_mhz",
    # "system_ul_interf_power_per_mhz",
    # "system_imt_antenna_gain",
    # "imt_system_antenna_gain",
]


# Carrega os resultados para diferentes cenários
dl_urban_results = Results.load_many_from_dir(
    dl_dir, only_latest=False,
    only_samples=attributes_to_plot
)

# ul_urban_results = Results.load_many_from_dir(
#     ul_dir, only_latest=True,
#     only_samples=attributes_to_plot
# )

# Combina todos os resultados em uma única lista
all_results = [
    *dl_urban_results,
    # *ul_urban_results,
]

# for result in [*dl_urban_results]:
#     result.system_inr = SampleList(
#         np.array(result.system_inr) - 10*np.log10(1.25)
#     )

#Salva um Sample List( nesse caso o agregado) como .csv , usando o nome e caminhos escolhido
def save_samplelist_as_csv(data, name: str, path: str):
    os.makedirs(path, exist_ok=True)
    df = pd.DataFrame({"samples": data})
    df.to_csv(os.path.join(path, name + ".csv"), index=False)

    # Acessa os resultados diretamente
dl_urb_r = dl_urban_results[0]
# ul_urb_r = ul_urban_results[0]

if None in [dl_urb_r, dl_urb_r]:
    raise Exception("Cannot aggregate results")

n_bs_sim = 19 * 3 * 3 * 7

print("Gerando resultados agregados")
aggregated_results_ra1rb1 = PostProcessor.aggregate_results(
    dl_samples=dl_urb_r.system_inr,
    ul_samples=dl_urb_r.system_inr,
    ul_tdd_factor=0,
    n_bs_sim=n_bs_sim,
    n_bs_actual=503000,
    n_aggregate=10000,
    )

# aggregated_results_ra2rb1 = PostProcessor.aggregate_results(
#     dl_samples=dl_urb_r.system_inr,
#     ul_samples=ul_urb_r.system_inr,
#     ul_tdd_factor=0.25,
#     n_bs_sim=n_bs_sim,
#     n_bs_actual=303300,
#     n_drops=10000
#     )

# save_samplelist_as_csv(aggregated_results_ra2rb1, f"fss_satD_ra2rb1", os.path.join(campaign_base_dir, "csv_agg"))

save_samplelist_as_csv(aggregated_results_ra1rb1, f"ressim_fss_stationD_BR_DL_ra2rb1", os.path.join(campaign_base_dir, "csv_files"))

# print(dl_dir)

fim = time.time()
tempo_execucao = fim - inicio
print(f"Tempo de execução: {tempo_execucao} segundos")
