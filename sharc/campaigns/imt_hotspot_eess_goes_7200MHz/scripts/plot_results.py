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
import re

def legend_gen(dir_name):
    #print(dir_name)
    link = re.search("_(dl|ul)", dir_name)
    if link is not None:
        link = link.group(1)
    else:
        return "None"
    
    t = re.search("_((sub){0,1}urban)", dir_name)
    if t is not None:
        t = t.group(1)
    else:
        return "None"
    
    return f"{link.upper()} {t.capitalize()} "

post_processor.add_plot_legend_generator(legend_gen)

attributes_to_plot = [
    "system_imt_antenna_gain",
    "imt_system_path_loss",
    "imt_system_antenna_gain",
    "system_dl_interf_power_per_mhz",
    "system_ul_interf_power_per_mhz",
]

def filter_fn(result_dir: str, is_suburban: bool) -> bool:
    sub = "_suburban" if is_suburban else "_urban"
    return "usa" in result_dir and sub in result_dir

dl_urban_results = Results.load_many_from_dir(
    dl_dir, only_latest=True,
    only_samples=attributes_to_plot,
    filter_fn=lambda x: filter_fn(x, False)
)

# dl_suburban_results = Results.load_many_from_dir(
#     dl_dir, only_latest=True,
#     only_samples=attributes_to_plot,
#     filter_fn=lambda x: filter_fn(x, True)
# )
ul_urban_results = Results.load_many_from_dir(
    ul_dir, only_latest=True,
    only_samples=attributes_to_plot,
    filter_fn=lambda x: filter_fn(x, False)
)
# ul_suburban_results = Results.load_many_from_dir(
#     ul_dir, only_latest=True,
#     only_samples=attributes_to_plot,
#     filter_fn=lambda x: filter_fn(x, True)
# )

# ^: list[Results]

all_results = [
    *dl_urban_results,
    # *dl_suburban_results,
    *ul_urban_results,
    # *ul_suburban_results
]
# ^: list[Results]

# transforming dBm / MHz to dB / kHz
# dBm -> dB means -30
# /MHz -> /kHz means -30
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
        plt = post_processor.get_plot_by_results_attribute_name(prop_name, plot_type=plot_type)
        if plt:
            plt.add_vline(
                interf_protection_criteria, line_dash="dash",
                name="0.1% criteria"
            )

system_dl_interf_power_plot = post_processor.get_plot_by_results_attribute_name("system_dl_interf_power_per_mhz")
system_ul_interf_power_plot = post_processor.get_plot_by_results_attribute_name("system_ul_interf_power_per_mhz")
aggregated_plot = None

if system_ul_interf_power_plot and system_dl_interf_power_plot:
    aggregated_plot = go.Figure()
    cutoff_percentage = 0.001
    next_tick = 1
    ticks_at = []
    while next_tick > cutoff_percentage:
        ticks_at.append(next_tick)
        next_tick /= 10
    ticks_at.append(cutoff_percentage)
    ticks_at.reverse()
    aggregated_plot.update_layout(
        title='CDF Plot for EESS Space Station receveid interference from Micro IMT in 7200 MHz',
        title_font_size=26,
        yaxis=dict(title="$\\text{P } (X < x)$",title_font_size=18, tickmode="array", tickvals=ticks_at, type="log", range=[np.log10(cutoff_percentage - cutoff_percentage/2), 0],tickfont=dict(size=18)),
        xaxis=dict(title="Interference [dBW/kHz]",title_font_size=18,tickmode="linear", dtick=5 , range=[-195,-159],tickfont=dict(size=18)), 
        legend=dict(
            x=0.5,  # Posição horizontal
            y=0.98,  # Posição vertical (98% do topo)
            bgcolor="rgba(255, 255, 255, 0.5)",  # Fundo semi-transparente
            bordercolor="rgba(0, 0, 0, 0.5)",  # Borda semi-transparente
            borderwidth=1,  # Largura da borda
            font=dict(size=18)  # Tamanho da legenda aumentado
        ),
        meta={"plot_type": "cdf"},
    )
    aggregated_ccdf_plot = go.Figure()
    aggregated_ccdf_plot.update_layout(
        title='CCDF Plot for EESS Space Station receveid interference from Micro IMT in 7200 MHz',
        title_font_size=26,
        yaxis=dict(title="$\\text{P } (X > x)$",title_font_size=18, tickmode="array", tickvals=ticks_at, type="log", range=[np.log10(cutoff_percentage - cutoff_percentage/2), 0],tickfont=dict(size=18)),
        xaxis=dict(title="Interference [dBW/kHz]",title_font_size=18,tickmode="linear", dtick=5 , range=[-195,-159],tickfont=dict(size=18)), 
        legend=dict(
            x=0.5,  # Posição horizontal 
            y=0.98,  # Posição vertical (98% do topo)
            bgcolor="rgba(255, 255, 255, 0.5)",  # Fundo semi-transparente
            bordercolor="rgba(0, 0, 0, 0.5)",  # Borda semi-transparente
            borderwidth=1,  # Largura da borda
            font=dict(size=20)  # Tamanho da legenda aumentado
        ),
        meta={"plot_type": "cdf"},
    )

    # we need to specify a common dir_name_contains substring so that we know which results we need to aggregate
    
    # Como há apenas um resultado de cada tipo, você pode acessá-los diretamente
    dl_urb_r = dl_urban_results[0]
    ul_urb_r = ul_urban_results[0]  # Pega o único resultado de UL urbano
    # ul_sub_r = ul_suburban_results[0]  # Pega o único resultado de UL suburban
    # dl_sub_r = dl_suburban_results[0]  # Pega o único resultado de DL suburban

    # Verifica se os resultados foram encontrados
    """print("Resultado de UL urbano:", ul_urb_r)
    print("Resultado de UL suburban:", ul_sub_r)
    print("Resultado de DL suburban:", dl_sub_r)"""

    # if None in [dl_sub_r, ul_sub_r, ul_urb_r]:
    #     raise Exception(f"Cannot aggregate {legend1} and {legend2}")

    n_bs_sim = 19 * 3 * 3 * 7

    rb = np.array([.01, .03])  # rb1, rb2
    ra_urban = np.array([.1, 0.45])  # ra1, ra2 urbano
    ra_suburban = np.array([0.05, 0.2])  # ra1, ra2 suburbano
    
    # area = 9867000  # US area (km²)
    area = 9867000

    ds_urb = 10
    ds_sub = 2.4

    for i in range(2):
        # Cálculo do número real de estações base macro
        # N_BS_ART ,  A * Ds*Ra*Rb
        
        # Ra1Rb1 (i=0)
        # Ra2Rb2 (i=1)
        
        n_bs_actual_urban = int(area * ds_urb * ra_urban[i] * rb[i])  #
        # n_bs_actual_suburban = int(area * ds_sub * ra_suburban[i] * rb[i])  #
      
        # print(f"Ra{i+1}Rb{i+1} : N_bs_urban = {n_bs_actual_urban} - N_bs_suburban = {n_bs_actual_suburban}")
        print(f"Ra{i+1}Rb{i+1} : N_bs_urban = {n_bs_actual_urban}")
        # Lista para armazenar os resultados
        aggregated_results_list = []

        aggregated_results = PostProcessor.aggregate_results(
            dl_samples=dl_urb_r.system_dl_interf_power_per_mhz,
            ul_samples=ul_urb_r.system_ul_interf_power_per_mhz,
            ul_tdd_factor=0.25,
            n_bs_sim=n_bs_sim,
            n_bs_actual=n_bs_actual_urban,
            n_drops=10000
        )
        # aggregated_results_sub = PostProcessor.aggregate_results(
        #     dl_samples=dl_sub_r.system_dl_interf_power_per_mhz,
        #     ul_samples=ul_sub_r.system_ul_interf_power_per_mhz,
        #     ul_tdd_factor=0.25,
        #     n_bs_sim=n_bs_sim,
        #     n_bs_actual=n_bs_actual_suburban,
        #     n_drops=10000
        # )
        # min_length = min(len(aggregated_results_sub), len(aggregated_results_urb))
   
   	#Somando os valores em escala linear 
        # aggregated_results = 10**(aggregated_results_sub[:min_length]/10) + 10**(aggregated_results_urb[:min_length]/10)
        
        #Retornando a escala Logaritma 
        # aggregated_results = 10*np.log10(aggregated_results_urb)
        #print(aggregated_results_sub[:3], aggregated_results_urb[:3], aggregated_results[:3])
        
        """
        #Testando apenas o Urbano e sub urbano por hora 
        x_urban, y_urban = PostProcessor.cdf_from(aggregated_results_urb)
        aggregated_plot.add_trace(
            go.Scatter(x=x_urban, y=y_urban, mode='lines', name=f'Aggregated Urban CDF Ra{i+1}Rb{i+1}',),
        )
        x_urban, y_urban = PostProcessor.ccdf_from(aggregated_results_urb)
        aggregated_ccdf_plot.add_trace(
            go.Scatter(x=x_urban, y=y_urban, mode='lines', name=f'Aggregated Urban CCDF Ra{i+1}Rb{i+1}',),
        )
        
        #Sub
        x_suburban, y_suburban = PostProcessor.cdf_from(aggregated_results_sub)
       # aggregated_plot.add_trace(
       #     go.Scatter(x=x_suburban, y=y_suburban, mode='lines', name=f'Aggregated Suburban CDF Ra{i+1}Rb{i+1}',),
       # )
        x_suburban, y_suburban = PostProcessor.ccdf_from(aggregated_results_sub)
       # aggregated_ccdf_plot.add_trace(
       #     go.Scatter(x=x_suburban, y=y_suburban, mode='lines', name=f'Aggregated Suburban CCDF Ra{i+1}Rb{i+1}',),
        #)
	
	"""
	#Agregado total
	
        x, y = PostProcessor.cdf_from(aggregated_results)

        aggregated_plot.add_trace(
            go.Scatter(x=x, y=y, mode='lines', name=f'Ra{i+1}Rb{i+1}',),
        )

        x, y = PostProcessor.ccdf_from(aggregated_results)
        aggregated_ccdf_plot.add_trace(
            go.Scatter(x=x, y=y, mode='lines', name=f'Ra{i+1}Rb{i+1}',),
        )
    
    """
    aggregated_plot.add_vline(
        interf_protection_criteria, line_dash="dash",
        name="0.1% criteria"
    )
    aggregated_ccdf_plot.add_vline(
        interf_protection_criteria, line_dash="dash",
        name="0.1% criteria"
    )
    """
    
    #Modificação para a linhar ser um traço padrão 
    # Adicionando a linha vertical com uma legenda
    aggregated_plot.add_trace(
        go.Scatter(
            x=[interf_protection_criteria, interf_protection_criteria],
            y=[0, 1],  # Ajuste os valores de y conforme necessário
            mode='lines',
            line=dict(dash='dash', color='gray'),  # Estilo da linha
            name='Protection criterion [-161 dBW/kHz, 0.1%]',  # Nome que aparecerá na legenda
            showlegend=True  # Garante que apareça na legenda
        )
    )

    # Adicionando a anotação (linha vertical)
    aggregated_plot.add_vline(
        x=interf_protection_criteria,
        line_dash="dash",
        line_color="gray",
        opacity=0.75  #  ajuste a opacidade
    )

    # Repetindo o mesmo para o gráfico CCDF
    aggregated_ccdf_plot.add_trace(
        go.Scatter(
            x=[interf_protection_criteria, interf_protection_criteria],
            y=[0, 1],  # Ajuste os valores de y conforme necessário
            mode='lines',
            line=dict(dash='dash', color='black'),  # Estilo da linha
            name='Protection criterion [-161 dBW/kHz, 0.1%]',  # Nome que aparecerá na legenda
            showlegend=True  # Garante que apareça na legenda
        )
    )

    # Adicionando a anotação (linha vertical)
    aggregated_ccdf_plot.add_vline(
        x=interf_protection_criteria,
        line_dash="dash",
        line_color="gray",
        opacity=0.75  # Ajuste a opacidade
    )
    #Linha horizontal no limite inferior

    aggregated_plot.add_hline(
        cutoff_percentage, line_dash="dash",
        name="limite inferior"
    )
    aggregated_ccdf_plot.add_hline(
        cutoff_percentage, line_dash="dash",
        name="limite inferior"
    )       
            
    
plots = [*post_processor.plots, aggregated_plot, aggregated_ccdf_plot]

PostProcessor.save_plots(
    os.path.join(campaign_base_dir, "output", "figs3"),
    plots,
    width = 1200,
    height= 800
)

plot_antenna_imt = PlotAntennaPattern("")

for plot in plots:
    plot.show()

"""
# Filtra e exibe apenas os gráficos Individuais coforme a escolha 
graf = "MHz from IMT DL"
for plot in post_processor.plots:
    title = plot.layout.title.text  # Obtém o título do gráfico
    if graf in title:  # Verifica se "graf" está no título
        plot.show()  # Exibe o gráfico no navegador
"""
        
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