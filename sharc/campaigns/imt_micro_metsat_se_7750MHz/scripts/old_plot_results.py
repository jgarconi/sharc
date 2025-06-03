import os
from pathlib import Path
from sharc.results import Results, SampleList
from sharc.post_processor import PostProcessor
import plotly.graph_objects as go
from sharc.parameters.parameters import Parameters
import glob
import numpy as np
from sharc.antenna.antenna_s465 import AntennaS465
from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt, PlotAntennaPattern

# Define o diretório base da campanha
campaign_base_dir = str((Path(__file__) / ".." / "..").resolve())
dl_dir = os.path.join(campaign_base_dir, "output_dl")
ul_dir = os.path.join(campaign_base_dir, "output_ul")

# Inicializa o pós-processador
post_processor = PostProcessor()

# Função para gerar legendas com base no nome do diretório
import re

def legend_gen(dir_name):
    sat = re.search("_(C|E)", dir_name)
    if sat is not None:
        sat = sat.group(1)
    else:
        return "None"
    
    sat2 = re.search("_(S|G)", dir_name)
    if sat2 is not None:
        sat2 = sat2.group(1)
    else:
        return "None"
    
    dist = re.search("_(100|200)", dir_name)
    if dist is not None:
        dist = dist.group(1)
    else:
        return "None"
    
    return f"Sat {sat.upper()}&{sat2.upper()} (DL - b = {dist.capitalize()} km)"

post_processor.add_plot_legend_generator(legend_gen)

# Atributos a serem plotados
attributes_to_plot = [
    "system_dl_interf_power_per_mhz"
]

# Função para filtrar resultados com base no tipo de ambiente (urbano/suburbano)
def filter_fn(result_dir: str, is_satP: bool) -> bool:
    sub = "_C" if is_satP else "_E"
    return sub in result_dir

# Carrega os resultados para diferentes cenários
dl_satP_results = Results.load_many_from_dir(
    dl_dir, only_latest=False,
    only_samples=attributes_to_plot,
    filter_fn=lambda x: filter_fn(x, True)
)

dl_satQ_results = Results.load_many_from_dir(
    dl_dir, only_latest=False,
    only_samples=attributes_to_plot,
    filter_fn=lambda x: filter_fn(x, False)
)

# Combina todos os resultados em uma única lista
all_results = [
    *dl_satP_results,
    *dl_satQ_results,
]

# transforming dBm / MHz to dBW / MHz
for result in all_results:
    result.system_dl_interf_power_per_mhz = SampleList(
        np.array(result.system_dl_interf_power_per_mhz) - 30
    )

# Adiciona os resultados ao pós-processador
post_processor.add_results(all_results)

post_processor.add_plots(
    post_processor.generate_ccdf_plots_from_results(
        all_results,
        cutoff_percentage=0.000016
    )
)

# Lista de atributos para adicionar linhas de critério de proteção
plots_to_add_vline = [
    "system_dl_interf_power_per_mhz"
]

# Critérios de proteção: linha horizontal, linha vertical, estilo tracejado
# interf_protection_criteria = {
#     "Protection criterion [-158 dBW/MHz, 20%]": [0.2, -158, "dash"],
#     "Protection criterion [-137 dBW/MHz, 0.0016%]": [0.000016, -137, "dot"],
# }

interf_protection_criteria = {
    "Protection criterion [-148 dBW/10MHz, 20%]": [0.2, -148, "dash"],
    "Protection criterion [-127 dBW/10MHz, 0.0016%]": [0.000016, -127, "dot"]
}

def add_protection_criteria(fig: go.Figure, interf_protection_criteria: dict) -> go.Figure:
    """
    Adiciona linhas de critério de proteção ao gráfico.
    
    Parâmetros:
    - fig: go.Figure -> Gráfico Plotly onde as linhas serão adicionadas.
    - interf_protection_criteria: dict -> Dicionário com critérios de proteção.
    """
    for legend_crite, val_crite in interf_protection_criteria.items():
        # Adiciona a linha vertical
        fig.add_trace(
            go.Scatter(
                x=[val_crite[1], val_crite[1]],
                y=[0, 1],
                mode='lines',
                line=dict(dash=val_crite[2], color="black"),
                name=legend_crite,
                showlegend=True
            )
        )

        # Adiciona a linha horizontal, se aplicável
        if val_crite[0] is not None:
            fig.add_hline(
                y=val_crite[0],
                line_dash=val_crite[2],
                line_color="gray",
            )

    return fig

def adjust_range_x(fig: go.Figure) -> go.Figure:
    """
    Ajusta automaticamente o eixo X do gráfico para melhor visualização.
    
    Parâmetros:
    - fig: go.Figure -> Gráfico Plotly a ser ajustado.
    """
    lim = fig.full_figure_for_development(warn=False)
    min_x_auto = lim.layout.xaxis.range[0] if lim.layout.xaxis.range else None
    max_x_auto = lim.layout.xaxis.range[1] if lim.layout.xaxis.range else None

    if min_x_auto is not None and max_x_auto is not None:
        fig.update_layout(xaxis=dict(range=[min_x_auto - 1, max_x_auto+1]))

    return fig

# Adiciona critérios de proteção aos gráficos selecionados
for prop_name in plots_to_add_vline:
    for plot_type in ["cdf", "ccdf"]:
        plt = post_processor.get_plot_by_results_attribute_name(prop_name, plot_type=plot_type)
        if plt:
            plt = add_protection_criteria(plt, interf_protection_criteria)
            plt = adjust_range_x(plt)

plots = [*post_processor.plots]

#                 azul       laranja    amarelo    roxo
custom_colors = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"]

# Aplica estilo personalizado às curvas
for plot in plots:
    curve_idx = 0  # Índice apenas para curvas reais (exclui critérios de proteção)
    for trace in plot.data:
        if isinstance(trace, go.Scatter) and "Protection criterion" not in trace.name:
            # Define o estilo da linha
            trace.line.dash = "solid" if curve_idx < 2 else "dash"
            # Aplica cor cíclica com base no índice
            trace.line.color = custom_colors[curve_idx % len(custom_colors)]
            curve_idx += 1

PostProcessor.save_plots(
    os.path.join(campaign_base_dir, "output"),
    plots,
    width=1200, height=1200,
)

# for plot in plots:
#     plot.show()