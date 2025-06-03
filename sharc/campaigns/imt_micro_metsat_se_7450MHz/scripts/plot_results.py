import os
from pathlib import Path
from sharc.results import Results, SampleList
from sharc.post_processor import PostProcessor
import plotly.graph_objects as go
from sharc.parameters.parameters import Parameters
import glob
import numpy as np

# Define o diretório base da campanha
campaign_base_dir = str((Path(__file__) / ".." / "..").resolve())
dl_dir = os.path.join(campaign_base_dir, "output_dl")

# Inicializa o pós-processador
post_processor = PostProcessor()

# Função para gerar legendas com base no nome do diretório
import re

def legend_gen(dir_name):
    
    sat = re.search("_(P|Q)_", dir_name)
    if sat is not None:
        sat = sat.group(1)
    else:
        return "None"
    
    dist = re.search("_(100km|200km)", dir_name)
    if dist is not None:
        dist = dist.group(1)
    else:
        return "None"
    
    return f"Sat {sat} (DL - b = {dist})"

post_processor.add_plot_legend_generator(legend_gen)

# Atributos a serem plotados
attributes_to_plot = [
    "system_dl_interf_power_per_mhz"
]

# Função para filtrar resultados com base no tipo de ambiente (urbano/suburbano)
def filter_fn(result_dir: str, is_satP: bool) -> bool:
    sub = "_P" if is_satP else "_Q"
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
# dBm = dBW - 30
for result in all_results:
    result.system_dl_interf_power_per_mhz = SampleList(
        np.array(result.system_dl_interf_power_per_mhz) - 30
    )
    
# Adiciona os resultados ao pós-processador
post_processor.add_results(all_results)

# Critérios de proteção: linha horizontal, linha vertical, estilo tracejado
interf_protection_criteria = {
    "Protection criterion [-154 dBW/MHz, 1%]": [0.01, -154, "dash"]
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
                line=dict(dash=val_crite[2], color="gray"),
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
        fig.update_layout(xaxis=dict(range=[-190, -140]))

    return fig

def legend_result(result: Results, pos_cenario: int = 2):
    """
    Retorna a string com a parte relevante da legenda.

    Parâmetros:
    - result : Results -> Objeto de resultados.
    - pos_cenario : int -> Índice da posição onde está a informação relevante na legenda.
      O padrão é 2, pois geralmente a informação importante para o cenário está na terceira posição
      da string gerada (por exemplo: 'Link Urbano/Suburbano diferença_cenario').

    Retorna:
    - Parte relevante da legenda como string.
    """
    # Usando o método que encontra as possíveis legendas atribuídas a esse resultado
    legends = post_processor.get_results_possible_legends(result)
    
    # Supondo que a legenda está separada por espaços
    return list(legends[0].values())[0].split()[pos_cenario]

cutoff_percentage = 0.001
next_tick = 1

ticks_major = []
ticks_minor = []

current_tick = next_tick
while current_tick > cutoff_percentage:
    ticks_major.append(current_tick)
    # Generate minor ticks for the current major interval:
    # They range from 10% to 90% of the current major value (step 10%)
    minor_ticks_for_interval = [current_tick * i for i in np.arange(1, .1, -0.1)]
    ticks_minor.extend(minor_ticks_for_interval)
    
    # Divide the current major tick by 10 for the next iteration
    current_tick /= 10  

ticks_major.append(cutoff_percentage)
ticks_major.reverse()
ticks_minor.append(cutoff_percentage)
ticks_minor.reverse()
# Create tick labels so that only major ticks are labeled
all_ticks = np.sort(np.unique(np.concatenate((ticks_major, ticks_minor))))
ticktext = [str(tick) if tick in ticks_major else "" for tick in all_ticks]

# Create tick labels so that only major ticks are labeled
all_ticks = np.sort(np.unique(np.concatenate((ticks_major, ticks_minor))))
ticktext = [f'10<sup><span style="font-size: 1.2em;">{int(np.log10(tick))}</span></sup>' if tick in ticks_major else "" for tick in all_ticks]

plots = go.Figure()
plots.update_layout(
            margin=dict( #Afastas os ticks dos eixos do grafico
                pad=10
            ),
            #  title={
            #     'text': "(b) Micro Cell - CCDF of Interference from IMT-DL",
            #     'x': 0.5,         # Centraliza horizontalmente
            #     'y': 0.947,
            #     'xanchor': 'center',
            #     'yanchor': 'top',  # Mantém a posição vertical
            #     'font': {
            #         'size': 32,
            #         'family': 'Arial',
            #         'color': 'black',
            #         'weight': 'bold'
            #     }
            # },
            xaxis_title="Interference Power (dBW/MHz)",
            yaxis_title="P (I>x)",
            yaxis=dict(tickmode="array", tickvals=all_ticks, type="log",
                    range=[np.log10(cutoff_percentage), 0],
                    ticktext=ticktext,
                    gridcolor="lightgray",
                    gridwidth=.5,
                    griddash="dot",                  
                    ),
            xaxis=dict(tickmode="linear",
                    dtick=5,
                    gridcolor="lightgray",
                    gridwidth=.5,
                    griddash="dot",
                    ),
            meta={"related_results_attribute": "Aggregated", "plot_type": "ccdf"},
            plot_bgcolor="white",
            paper_bgcolor="white",
            font=dict(
                family="Arial, sans-serif",
                size=16,         # Base font size for all text
                color="black"    # Text color
            ),
            shapes=[
                dict(
                    type="rect",
                    xref="paper",
                    yref="paper",
                    x0=0,
                    y0=cutoff_percentage,
                    x1=1,
                    y1=1,
                    line=dict(
                        color="black",
                        width=1
                    ),
                    fillcolor="rgba(0,0,0,0)"  # transparent fill
                )
            ],
            legend=dict(
                x=0.34,          # x position (95% from the left)
                y=0.235,          # y position (95% from the bottom)
                xanchor='right', # anchor the legend's right side at x=0.95
                yanchor='top',   # anchor the legend's top at y=0.95
                bgcolor='rgba(255,255,255,0.5)',  # Optional: semi-transparent white background
                bordercolor='gray',              # Optional: border color for better separation
                borderwidth=1                     # Optional: border width in pixels
            )
        )


#         azul       laranja    amarelo    roxo
color = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"]
for i,dl_urb_r in enumerate(all_results):

    legenda = post_processor.get_results_possible_legends(dl_urb_r)[0]['legend']
    sat = legend_result(dl_urb_r,1)
    # dist = legend_result(dl_urb_r,6)
    # print(f'{dist} --- {sat}')

    line_style = {"P" : "solid","Q" : "dash"}

    x, y = PostProcessor.ccdf_from(dl_urb_r.system_dl_interf_power_per_mhz,n_bins=200)
    plots.add_trace(
        go.Scatter(x=x, y=y, mode='lines',name=legenda, line_width=3, line_dash= line_style[sat], line_color = color[i])
    )

plots = add_protection_criteria(plots, interf_protection_criteria)
plots = adjust_range_x(plots)

#Gambiarra(Exite o minor no ploty que talze possa funcionar )
# -------- Adicionando linhas para simular os subtick (Jeito tosco)
lim_infe = -190 #limite inferior no eixo x( e o valor que esta na adjust_range_x)
lim_supe = -140 #limite superior no eixo x( e o valor que esta na adjust_range_x)

tamanho_maior = 0.4    
tamanho_menor= 0.20  

# Adiciona as linhas horizontais (a Esquerda e a Direita)
for tick in all_ticks :
    infe = 1
    for lim in (lim_infe,lim_supe):
        plots.add_trace(
            go.Scatter(
                x=[lim, lim + infe*(tamanho_maior if tick in ticks_major else tamanho_menor)],
                y=[tick,tick],
                mode='lines',
                line=dict(dash="solid", width=1.5,color="gray"),
                name = "subtick",
                showlegend=False
            )
        )
        infe = -1

#Linhas Verticais
for i in range(lim_infe,lim_supe,5):
    ymin = min(all_ticks)
    plots.add_trace(
            go.Scatter(
                x=[i,i],
                y=[0,ymin + ymin*tamanho_menor*0.35], #ultimo numero e apenas um fator de tamanho que escolho(reduz) 
                mode='lines',
                line=dict(dash="solid", width=1.5,color="gray"),
                name = "subtick",
                showlegend=False
            )
        )


# Salva os gráficos 

PostProcessor.save_plots(
    os.path.join(campaign_base_dir, "output"),
    [plots],
    width=1200, height=800,
)

# Exibe o gráfico CCDF agregado
# plots.show()