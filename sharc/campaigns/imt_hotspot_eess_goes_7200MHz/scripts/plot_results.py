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

post_processor = PostProcessor()

# Add a legend to results in folder that match the pattern
# This could easily come from a config file
post_processor\
    .add_plot_legend_pattern(
        dir_name_contains="_afr_1_cluster_dl",
        legend="Hotspot DL 1 clusters AFR"
    ).add_plot_legend_pattern(
        dir_name_contains="_br_1_cluster_ul",
        legend="Hotspot UL 1 clusters BR"
    ).add_plot_legend_pattern(
        dir_name_contains="_afr_7_cluster_dl",
        legend="Hotspot DL 7 clusters AFR"
    ).add_plot_legend_pattern(
        dir_name_contains="_br_7_cluster_ul",
        legend="Hotspot UL 7 clusters BR"
    ).add_plot_legend_pattern(
        dir_name_contains="_afr_1_cluster_ul",
        legend="Hotspot UL 1 clusters AFR"
    ).add_plot_legend_pattern(
        dir_name_contains="_sa_1_cluster_dl",
        legend="Hotspot DL 1 clusters SA"
    ).add_plot_legend_pattern(
        dir_name_contains="_afr_7_cluster_ul",
        legend="Hotspot UL 7 clusters AFR"
    ).add_plot_legend_pattern(
        dir_name_contains="_sa_7_cluster_dl",
        legend="Hotspot DL 7 clusters SA"
    ).add_plot_legend_pattern(
        dir_name_contains="_br_1_cluster_dl",
        legend="Hotspot DL 1 clusters BR"
    ).add_plot_legend_pattern(
        dir_name_contains="_sa_1_cluster_ul",
        legend="Hotspot UL 1 clusters SA"
    ).add_plot_legend_pattern(
        dir_name_contains="_br_7_cluster_dl",
        legend="Hotspot DL 7 clusters BR"
    ).add_plot_legend_pattern(
        dir_name_contains="_sa_7_cluster_ul",
        legend="Hotspot UL 7 clusters SA"
    )

attributes_to_plot = [
    # "system_imt_antenna_gain",
    # "imt_system_path_loss",
    # "imt_system_antenna_gain",
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
    os.path.join(campaign_base_dir, "output_ul"), only_latest=True,
    only_samples=attributes_to_plot,
    filter_fn=filter_fn
)
# ^: typing.List[Results]

all_results = [
    *dl_results,
    *ul_results
]

# transforming dBm / MHz to dB / kHz
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
            

        n_bs_sim = 19*3*3
        if "_7_" in legend1["dir_name_contains"]:
            n_bs_sim = n_bs_sim * 7

        aggregated_results = PostProcessor.aggregate_results(
            dl_samples=dl_r.system_dl_interf_power,
            ul_samples=ul_r.system_ul_interf_power,
            ul_tdd_factor=0.25,
            n_bs_sim=n_bs_sim,
            # n_bs_actual=19*3*3*7                                                                                                                           
            # n_bs_actual=363300,                                                                                                                            
            # n_bs_actual=1035900,                                                                                                                           
            # n_bs_actual=723300,                                                                                                                            
            # n_bs_actual=1364850,                                                                                                                           
            n_bs_actual=4391550,                                                                                                                             
            # n_bs_actual=1904850
        )
        x, y = PostProcessor.cdf_from(aggregated_results)

        aggregated_plot.add_trace(
            go.Scatter(x=x, y=y, mode='lines', name=f'{legend1["legend"]}',),
        )

    # Add a protection criteria line:
    # dB to dBm (+ 30)
    # the following conversion makes the criteria more strict, so there may not be a problem
    interf_protection_criteria = -161

    aggregated_plot.add_vline(
        interf_protection_criteria, line_dash="dash",
        name="1% criteria"
    )

# i = 0

# for result in all_results:
#     if "_10000m_" not in result.output_directory:
#         continue

#     params_file = glob.glob(result.output_directory + "/*.yaml")[0]
#     params = Parameters()
#     params.set_file_name(params_file)
#     params.read_params()

#     # TODO: use antenna factory here if it ever exists
#     legend = post_processor.get_results_possible_legends(result)[0]
#     if params.single_earth_station.antenna.pattern == "ITU-R S.465":
#         antenna = AntennaS465(params.single_earth_station.antenna.itu_r_s_465)
#         PostProcessor.generate_antenna_radiation_pattern_plot(antenna, antenna_legends[i]).show()
#     if i == 0:
#         antenna_bs = AntennaBeamformingImt(
#             params.imt.bs.antenna.get_antenna_parameters(),
#             0,
#             0,
#             # -params.imt.bs.antenna.downtilt
#         )
#         antenna_ue = AntennaBeamformingImt(
#             params.imt.ue.antenna.get_antenna_parameters(),
#             0,
#             0
#         )

#     i += 1


# Show a single plot:

PostProcessor.save_plots(
    os.path.join(campaign_base_dir, "output", "figs"),
    [*post_processor.plots, aggregated_plot],
)

# if aggregated_plot:
#     aggregated_plot.show()


plot_antenna_imt = PlotAntennaPattern("")

# Plot BS TX radiation patterns
# f = plot_antenna_imt.plot_element_pattern(antenna_bs, "BS", "ELEMENT")
# # f.savefig(figs_dir + "BS_element.pdf", bbox_inches='tight')
# f = plot_antenna_imt.plot_element_pattern(antenna_bs, "TX", "ARRAY")
# # f.savefig(figs_dir + "BS_array.pdf", bbox_inches='tight')

# # Plot UE TX radiation patterns
# plot_antenna_imt.plot_element_pattern(antenna_ue, "UE", "ELEMENT")
# plot_antenna_imt.plot_element_pattern(antenna_ue, "UE", "ARRAY")

# Plot every plot:
# for plot in plots:
#     plot.show()

# full_results = ""

# for result in all_results:
#     # This generates the mean, median, variance, etc
#     stats = PostProcessor.generate_statistics(
#         result=result
#     ).write_to_results_dir()

#     full_results += str(stats) + "\n"
#     # # do whatever you want here:
#     # if "fspl_45deg" in stats.results_output_dir:
#     #     get some stat and do something

# with open(dl_dir + "/stats.txt", "w") as f:
#     f.write(full_results)
