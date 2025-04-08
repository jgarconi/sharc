import os
from pathlib import Path
from sharc.results import Results, SampleList
from sharc.post_processor import PostProcessor
import plotly.graph_objects as go

import re
import numpy as np
from sharc.antenna.antenna_beamforming_imt import  PlotAntennaPattern

campaign_base_dir = str((Path(__file__) / ".." / "..").resolve())
dl_dir = os.path.join(campaign_base_dir, "output_dl")
ul_dir = os.path.join(campaign_base_dir, "output_ul")

post_processor = PostProcessor()

def legend_gen(dir_name):
    link = re.search("_(dl|ul)", dir_name)
    if link is not None:
        link = link.group(1)
    else:
        return "None"
    
    t = re.search("_(100km|550km)_", dir_name)
    if t is not None:
        t = t.group(1)
    else:
        return "None"

    return f"{link.upper()} {t.capitalize()} "

post_processor.add_plot_legend_generator(legend_gen)

attributes_to_plot = [
    # "system_imt_antenna_gain",
    # "imt_system_path_loss",
    # "imt_system_antenna_gain",
    # "imt_dl_tx_power",
    # "imt_dl_tx_power_density",
    # "imt_ul_inr",
    # "imt_dl_inr",
    "system_dl_interf_power_per_mhz",
    "system_ul_interf_power_per_mhz"
    # "system_inr"
]

def filter_fn(x, is_100km):
    return ("100km" in x if is_100km else "550km" in x)

def load_results(base_dir, is_100km):
    return Results.load_many_from_dir(
        base_dir, only_latest=True,
        only_samples=attributes_to_plot,
        filter_fn=lambda x: filter_fn(x, is_100km)
    )

dl_results_100km, ul_results_100km = load_results(dl_dir, True), load_results(ul_dir, True)
dl_results_550km, ul_results_550km = load_results(dl_dir, False), load_results(ul_dir, False)

all_results = [
    *dl_results_100km,
    *ul_results_100km,
    *dl_results_550km,
    *ul_results_550km,
]

# dBm -> - 30 -> MHz
#  dB -> - 30 -> kHz
for result in all_results:
    result.system_dl_interf_power_per_mhz = SampleList(
        np.array(result.system_dl_interf_power_per_mhz) - 30 - 30
    )
    result.system_ul_interf_power_per_mhz = SampleList(
        np.array(result.system_ul_interf_power_per_mhz) - 30 - 30
    )

post_processor.add_results(all_results)

post_processor.add_plots(
    post_processor.generate_ccdf_plots_from_results(
        all_results,
        cutoff_percentage=0.001
    )
)
post_processor.add_plots(
    post_processor.generate_cdf_plots_from_results(
        all_results,
    )
)

plots_to_add_vline = [
    "system_inr",
    "imt_ul_inr",
    "imt_dl_inr",
]

interf_protection_criteria = {
    "Protection criterion [-177 dBW/kHz, 0,1%]": [None, -177, "dash"],
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
                line_color="black",
                annotation_text=f"{legend_crite}",
                annotation_position="top left"
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
        fig.update_layout(xaxis=dict(range=[min_x_auto - 1, max_x_auto + 4]))

    return fig

# Adiciona critérios de proteção aos gráficos selecionados
for prop_name in plots_to_add_vline:
    for plot_type in ["cdf", "ccdf"]:
        plt = post_processor.get_plot_by_results_attribute_name(prop_name, plot_type=plot_type)
        if plt:
            plt = add_protection_criteria(plt, interf_protection_criteria)
            plt = adjust_range_x(plt)

system_dl_interf_power_plot = post_processor.get_plot_by_results_attribute_name("system_dl_interf_power_per_mhz")
system_ul_interf_power_plot = post_processor.get_plot_by_results_attribute_name("system_ul_interf_power_per_mhz")
aggregated_plot_100km = None
aggregated_ccdf_plot_100km = None
aggregated_plot_550km = None
aggregated_ccdf_plot_550km = None

def create_aggregated_plot(title, plot_type):
    return go.Figure().update_layout(
        title=title,
        title_font_size=22,
        yaxis=dict(title=f"$\\text{{P }} (X {'>' if plot_type == 'ccdf' else '<'} x)$",
                   title_font_size=18, tickmode="array", tickvals=ticks_at,
                   type="log", range=[np.log10(cutoff_percentage - cutoff_percentage / 2), 0],
                   tickfont=dict(size=18)),
        xaxis=dict(title="Interference [dBW/kHz]", title_font_size=18,
                   tickmode="linear", dtick=5, tickfont=dict(size=18)),
        legend=dict(x=0.55, y=0.1, bgcolor="rgba(255, 255, 255, 0.5)",
                    bordercolor="rgba(0, 0, 0, 0.5)", borderwidth=1,
                    font=dict(size=15), xanchor="center", yanchor="auto"),
        meta={"plot_type": plot_type},
    )

if system_ul_interf_power_plot and system_dl_interf_power_plot:
    aggregated_plot_100km = go.Figure()
    aggregated_ccdf_plot_100km = go.Figure()
    aggregated_plot_550km = go.Figure()
    aggregated_ccdf_plot_550km = go.Figure()
    comparison_ccdf_plot = go.Figure()

    cutoff_percentage = 0.001
    next_tick = 1
    ticks_at = []

    while next_tick > cutoff_percentage:
        ticks_at.append(next_tick)
        next_tick /= 10

    ticks_at.append(cutoff_percentage)
    ticks_at.reverse()

    aggregated_plot_100km = create_aggregated_plot(
        "CDF Aggregated Plot for SRS Space Station received interference (100 km²)", "cdf"
        )
    aggregated_ccdf_plot_100km = create_aggregated_plot(
        "CCDF Aggregated Plot for SRS Space Station received interference (100 km²)", "ccdf"
        )
    aggregated_plot_550km = create_aggregated_plot(
        "CDF Aggregated Plot for SRS Space Station received interference (550 km²)", "cdf"
        )
    aggregated_ccdf_plot_550km = create_aggregated_plot(
        "CCDF Aggregated Plot for SRS Space Station received interference (550 km²)", "ccdf"
        )
    comparison_ccdf_plot = create_aggregated_plot(
        "Comparison of CCDF Aggregated Plot for SRS Space Station receveid interference from Micro IMT in 7182 MHz", "ccdf"
        )
    
    post_processor.add_plots(
    post_processor.generate_ccdf_plots_from_results(
        all_results,
        cutoff_percentage=0.001
        )
    )
    post_processor.add_plots(
        post_processor.generate_cdf_plots_from_results(
            all_results,
        )
    )

    # Como há apenas um resultado de cada tipo, você pode acessá-los diretamente
    dl_100km = dl_results_100km[0]
    ul_100km = ul_results_100km[0]
    
    dl_550km = dl_results_550km[0]
    ul_550km = ul_results_550km[0]  # Pega o único resultado de UL urbano

    n_bs_sim = 19 * 3 * 3 * 7

    # NOTE: From Table 13 Annex 4.15 for micro cells
    ra_urban = np.array([.05, .1])  # ra1, ra2 urbano micro
    
    area = np.array([1517697.32, 13238154.59])

    # Deployment density for IMT BS (BS/km²)
    # NOTE: From Table 13 Annex 4.15 for micro cells
    ds_urb = 30

    for i in range(len(area)):  # Percorre cada área
        if area[i] > 3500000:
            rb_values = np.array([0.01, 0.03])  # rb1, rb2
        else:
            rb_values = np.array([0.01, 0.05])  # rb1, rb2

        for j in range(len(ra_urban)):  # Percorre cada combinação de ra e rb
            n_bs_actual_urban = int(area[i] * ds_urb * ra_urban[j] * rb_values[j])

            print(f"Área {i+1} = {area[i]}, Ra{j+1}Rb{j+1} : N_bs_urban = {n_bs_actual_urban}")

            aggregated_results_100km = PostProcessor.aggregate_results(
                dl_samples=dl_100km.system_dl_interf_power_per_mhz,
                ul_samples=ul_100km.system_ul_interf_power_per_mhz,
                ul_tdd_factor=0.25,
                n_bs_sim=n_bs_sim,
                n_bs_actual=n_bs_actual_urban,
                n_drops=10000
            )

            aggregated_results_550km = PostProcessor.aggregate_results(
                dl_samples=dl_550km.system_dl_interf_power_per_mhz,
                ul_samples=ul_550km.system_ul_interf_power_per_mhz,
                ul_tdd_factor=0.25,
                n_bs_sim=n_bs_sim,
                n_bs_actual=n_bs_actual_urban,
                n_drops=10000
            )

        #Agregado total	- 100 km²
        x_100km, y_100km = PostProcessor.cdf_from(aggregated_results_100km)
        # aggregated_plot_100km = go.Figure()
        aggregated_plot_100km.add_trace(
            go.Scatter(x=x_100km, y=y_100km, mode='lines', name=f'Ra{i+1}Rb{i+1}, area 100 km²')
        )
        x_ccdf_100km, y_ccdf_100km = PostProcessor.ccdf_from(aggregated_results_100km)
        # aggregated_ccdf_plot_100km = go.Figure()
        aggregated_ccdf_plot_100km.add_trace(
            go.Scatter(x=x_ccdf_100km, y=y_ccdf_100km, mode='lines', name=f'Ra{i+1}Rb{i+1}, area 100 km²')
        )

        #Agregado total	- 550 km²
        x_550km, y_550km = PostProcessor.cdf_from(aggregated_results_550km)
        # aggregated_plot_550km = go.Figure()
        aggregated_plot_550km.add_trace(
            go.Scatter(x=x_550km, y=y_550km, mode='lines', name=f'Ra{i+1}Rb{i+1}, area 550 km²')
        )
        x_ccdf_550km, y_ccdf_550km = PostProcessor.ccdf_from(aggregated_results_550km)
        # aggregated_ccdf_plot_5500km = go.Figure()
        aggregated_ccdf_plot_550km.add_trace(
            go.Scatter(x=x_ccdf_550km, y=y_ccdf_550km, mode='lines', name=f'Ra{i+1}Rb{i+1}, area 550 km²')
        )

        # Comparação entre os dois agregados
        comparison_ccdf_plot.add_trace(
            go.Scatter(x=x_ccdf_100km, y=y_ccdf_100km, mode='lines', name=f'Ra{i+1}Rb{i+1}, area 100 km²')
        )
        comparison_ccdf_plot.add_trace(
            go.Scatter(x=x_ccdf_550km, y=y_ccdf_550km, mode='lines', name=f'Ra{i+1}Rb{i+1}, area 550 km²')
        )
    # Ajusta eixos e critérios de proteção
    aggregated_plot_100km = add_protection_criteria(aggregated_plot_100km, interf_protection_criteria)
    aggregated_ccdf_plot_100km = add_protection_criteria(aggregated_ccdf_plot_100km, interf_protection_criteria)
    aggregated_plot_550km = add_protection_criteria(aggregated_plot_550km, interf_protection_criteria)
    aggregated_ccdf_plot_550km = add_protection_criteria(aggregated_ccdf_plot_550km, interf_protection_criteria)
    comparison_ccdf_plot = add_protection_criteria(comparison_ccdf_plot, interf_protection_criteria)

    # Ajusta o range dos eixos
    aggregated_plot_100km = adjust_range_x(aggregated_plot_100km)
    aggregated_ccdf_plot_100km = adjust_range_x(aggregated_ccdf_plot_100km)
    aggregated_plot_550km = adjust_range_x(aggregated_plot_550km)
    aggregated_ccdf_plot_550km = adjust_range_x(aggregated_ccdf_plot_550km)
    comparison_ccdf_plot = adjust_range_x(comparison_ccdf_plot)

plots = [
    *post_processor.plots,
    aggregated_plot_100km,
    aggregated_ccdf_plot_100km,
    aggregated_plot_550km,
    aggregated_ccdf_plot_550km,
    comparison_ccdf_plot
]

# PostProcessor.save_plots(
#     os.path.join(campaign_base_dir, "output", "figs"),
#     plots,
#     width = 1200,
#     height= 800
# )

plot_antenna_imt = PlotAntennaPattern("")

for plot in plots:
    plot.show()
