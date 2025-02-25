import os
from pathlib import Path
from sharc.results import Results
import plotly.graph_objects as go
from sharc.post_processor import PostProcessor
import pandas

post_processor = PostProcessor()

# Add a legend to results in folder that match the pattern
# This could easily come from a config file
post_processor\
    .add_plot_legend_pattern(
        dir_name_contains="small_dl",
        legend="Small Beam DL"
    ).add_plot_legend_pattern(
        dir_name_contains="small_ul",
        legend="Small Beam UL"
    ).add_plot_legend_pattern(
        dir_name_contains="large_dl",
        legend="Large Beam DL"
    ).add_plot_legend_pattern(
        dir_name_contains="large_ul",
        legend="Large Beam UL"
    )

campaign_base_dir = str((Path(__file__) / ".." / "..").resolve())

many_results = Results.load_many_from_dir(os.path.join(campaign_base_dir, "output"), only_latest=True)

post_processor.add_results(many_results)

plots = post_processor.generate_cdf_plots_from_results(
    many_results
)
post_processor.add_plots(plots)

# This function aggregates IMT downlink and uplink
aggregated_results = PostProcessor.aggregate_results(
    sm_dl_samples=post_processor.get_results_by_output_dir("small_dl").system_inr,
    sm_ul_samples=post_processor.get_results_by_output_dir("small_ul").system_inr,
    lg_dl_samples=post_processor.get_results_by_output_dir("large_dl").system_inr,
    lg_ul_samples=post_processor.get_results_by_output_dir("large_ul").system_inr,
    ul_tdd_factor=0.25,
    # SF is not exactly 1, but approx
    sm_n_bs_sim=399,
    sm_n_bs_actual=420,
    lg_n_bs_sim=399,
    lg_n_bs_actual=816992
)

# Add a protection criteria line:
# protection_criteria = int

# post_processor\
#      .get_plot_by_results_attribute_name("system_dl_interf_power")\
    #  .add_vline(protection_criteria, line_dash="dash")

# Show a single plot:
relevant = post_processor\
    .get_plot_by_results_attribute_name("system_inr")

aggr_x, aggr_y = PostProcessor.cdf_from(aggregated_results)

relevant.add_trace(
    go.Scatter(x=aggr_x, y=aggr_y, mode='lines', name='Aggregate interference',),
)

compare_to = pandas.read_csv(
    os.path.join(campaign_base_dir, "comparison", "contribution_20.csv"),
    skiprows=1
)

comp_x, comp_y = (compare_to.iloc[:, 0], compare_to.iloc[:, 1])

relevant.add_trace(
    go.Scatter(x=comp_x, y=comp_y, mode='lines', name='Reference INR',),
)

# Adiciona uma linha vertical em x = -6 usando Scatter para aparecer na legenda
relevant.add_trace(
    go.Scatter(
        x=[-6, -6], y=[0, 1],  # Linha vertical de y=0 até y=1
        mode="lines",
        line=dict(color="red", width=2, dash="dash"),
        name="Protection Criterion"
    )
)

relevant.show()

for result in many_results:
    # This generates the mean, median, variance, etc
    stats = PostProcessor.generate_statistics(
        result=result
    ).write_to_results_dir()
