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

campaign_base_dir = str((Path(__file__) / ".." / "..").resolve())
dl_dir = os.path.join(campaign_base_dir, "output_dl")
ul_dir = os.path.join(campaign_base_dir, "output_ul")

post_processor = PostProcessor()

# Add a legend to results in folder that match the pattern
# This could easily come from a config file
post_processor\
    .add_plot_legend_pattern(
        dir_name_contains="_eua_dl",
        legend="Hotspot DL US"
    ).add_plot_legend_pattern(
        dir_name_contains="_eua_ul",
        legend="Hotspot UL US"
    )

attributes_to_plot = [
    "system_imt_antenna_gain",
    "imt_system_path_loss",
    "imt_system_antenna_gain",
    "system_dl_interf_power_per_mhz",
    "system_ul_interf_power_per_mhz",
]

def filter_fn(result_dir: str) -> bool:
    # return "10000m" in result_dir
    return True

dl_results = Results.load_many_from_dir(
    dl_dir, only_latest=True,
    only_samples=attributes_to_plot,
    filter_fn=filter_fn
)
ul_results = Results.load_many_from_dir(
    ul_dir, only_latest=True,
    only_samples=attributes_to_plot,
    filter_fn=filter_fn
)
# ^: typing.List[Results]

all_results = [
    *dl_results,
    *ul_results
]

# NOTE: dBm to dBW (-30) and MHz to kHz (-30)
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
        all_results
    )
)

# Add a protection criteria line:
# dB to dBm (+ 30)
# the following conversion makes the criteria more strict, so there may not be a problem
plots_to_add_vline = [
    "system_ul_interf_power_per_mhz",
    "system_dl_interf_power_per_mhz"
]
interf_protection_criteria = -161

for prop_name in plots_to_add_vline:
    for plot_type in ["cdf", "ccdf"]:
        plt = post_processor\
            .get_plot_by_results_attribute_name(prop_name, plot_type=plot_type)
        if plt:
            plt.add_vline(
                interf_protection_criteria, line_dash="dash",
                name="1% criteria"
            )

system_dl_interf_power_plot = post_processor\
    .get_plot_by_results_attribute_name("system_dl_interf_power_per_mhz")

system_ul_interf_power_plot = post_processor\
    .get_plot_by_results_attribute_name("system_ul_interf_power_per_mhz")
aggregated_plot = None

if system_ul_interf_power_plot and system_dl_interf_power_plot:
    aggregated_plot = go.Figure()
    aggregated_plot.update_layout(
        title=f'CDF Plot for aggregated Spectral Power Density',
        xaxis_title="Spectral Power Density (dB/KHz)",
        yaxis_title="CDF",
        yaxis=dict(tickmode="array", tickvals=[0, 0.25, 0.5, 0.75, 1]),
        xaxis=dict(tickmode="linear", dtick=5),
        legend_title="Labels",
        meta={"plot_type": "cdf"},
    )
    
    aggregated_ccdf_plot = go.Figure()
    cutoff_percentage = 0.001
    next_tick = 1
    ticks_at = []
    while next_tick > cutoff_percentage:
        ticks_at.append(next_tick)
        next_tick /= 10
    ticks_at.append(cutoff_percentage)
    ticks_at.reverse()
    aggregated_ccdf_plot.update_layout(
        title=f'CCDF Plot for aggregated Spectral Power Density from Interference',
        xaxis_title="Interference (dBW/kHz)",
        yaxis_title="CCDF",
        yaxis=dict(tickmode="array", tickvals=ticks_at, type="log", range=[np.log10(cutoff_percentage), 0]),
        xaxis=dict(tickmode="linear", dtick=5),
        legend_title="Labels",
        meta={"plot_type": "ccdf"},
    )

    for dl_r in dl_results:
        legend1 = post_processor.get_results_possible_legends(dl_r)[0]
        ul_r = None
        for maybe in ul_results:
            legend2 = post_processor.get_results_possible_legends(maybe)[0]
            if legend1["dir_name_contains"][:-3] == legend2["dir_name_contains"][:-3]:
                ul_r = maybe
                break
        if ul_r is None:
            raise Exception(f"Cannot aggregate {legend1} and {legend2}")
            # continue
            
        n_bs_sim = 19*7*3*3

        rb = np.array([.01, .03])  # rb1, rb2
        ra_urban = np.array([.05, 0.1])  # ra1, ra2 urbano
        ra_suburban = np.array([0, 0])  # ra1, ra2 suburbano
        area = 9867000  # US area (km²)

        ds_urb = 30
        ds_sub = 2.4

        # Lista para armazenar os resultados
        aggregated_results_list = []

        # Loops para todas as combinações de rb com ra_urban e ra_suburban
        for i in range(len(rb)):  
            rb_i = rb[i]
            ra_u = ra_urban[i]
            ra_s = ra_suburban[i]

            # Cálculo do número real de estações base macro
            n_bs_actual = area * rb_i * ((ds_urb * ra_u) + (ds_sub * ra_s))

            # Chamando a função aggregate_results
            aggregated_results = PostProcessor.aggregate_results(
                dl_samples=dl_r.system_dl_interf_power_per_mhz,
                ul_samples=ul_r.system_ul_interf_power_per_mhz,
                ul_tdd_factor=0.25,
                n_bs_sim=n_bs_sim,
                n_bs_actual=n_bs_actual,  # Valor calculado dinamicamente
            )

            # Armazena cada resultado
            aggregated_results_list.append({
                "rb": rb_i,
                "ra_urban": ra_u,
                "ra_suburban": ra_s,
                "n_bs_actual": n_bs_actual,
                "aggregated_results": aggregated_results
            })

        # Agora adicionamos os dois conjuntos de resultados ao gráfico
        for result in aggregated_results_list:
            x, y = PostProcessor.cdf_from(result["aggregated_results"])
            aggregated_plot.add_trace(
                go.Scatter(x=x, y=y, mode='lines', name=f'Aggregated CDF rb={result["rb"]}',),
            )

            x, y = PostProcessor.ccdf_from(result["aggregated_results"])
            aggregated_ccdf_plot.add_trace(
                go.Scatter(x=x, y=y, mode='lines', name=f'Aggregated CCDF rb={result["rb"]}',),
            )

    # Add a protection criteria line:
    # dB to dBm (+ 30)
    # the following conversion makes the criteria more strict, so there may not be a problem
    interf_protection_criteria = -161

    aggregated_plot.add_vline(
        interf_protection_criteria, line_dash="dash",
        name="1% criteria"
    )

    aggregated_ccdf_plot.add_vline(
        interf_protection_criteria, line_dash="dash",
        name="1% criteria"
    )

plots = [*post_processor.plots, aggregated_plot, aggregated_ccdf_plot]
# Plot every plot:
for plot in plots:
    plot.show()

# Save plots
# PostProcessor.save_plots(
#     os.path.join(campaign_base_dir, "output", "figs"),
#     plots,
# )

plot_antenna_imt = PlotAntennaPattern("")