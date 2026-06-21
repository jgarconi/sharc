import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt

from sharc.parameters.imt.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_antenna import ParametersAntenna
from antenna_beamforming_imt import AntennaBeamformingImt 
from antenna_f1245_fs import Atenna_f1245_fs

# Potência de Transmissão da BS IMT [dBm]
POT_BS_IMT = 30

# Configuração do IMT
bs_param = ParametersAntennaImt()
bs_param.adjacent_antenna_model = "BEAMFORMING"
bs_param.element_pattern = "M2101"
bs_param.element_max_g = 6.4
bs_param.element_phi_3db = 90
bs_param.element_theta_3db = 65
bs_param.element_am = 30
bs_param.element_sla_v = 30
bs_param.n_rows = 8
bs_param.n_columns = 16
bs_param.element_horiz_spacing = 0.5
bs_param.element_vert_spacing = 2.1
bs_param.multiplication_factor = 12
bs_param.downtilt = 6
bs_param.normalization = False
bs_param.minimum_array_gain = -200

bs_param.subarray.is_enabled = True
bs_param.subarray.n_rows = 3
bs_param.subarray.element_vert_spacing = 0.7
bs_param.subarray.eletrical_downtilt = 3

# Parâmetros físicos da antena
#bs_azimuth_fisico = 0.0 
#bs_antenna = AntennaBeamformingImt(bs_param.get_antenna_parameters(), bs_azimuth_fisico, -bs_param.downtilt)

# Configuração FS 
fs_par = ParametersAntenna()
fs_par.gain = 36.0
fs_par.frequency = 8000
fs_par.diameter = 2.0
fs_antenna = Atenna_f1245_fs(fs_par)
fs_antenna.add_beam(0, 0) 

# skiprows=1 pula a primeira linha (cabeçalho original)
df = pd.read_csv('/home/juliana/Downloads/ep_to_cells.csv')

# 2. LIMPEZA CRUCIAL: Remove espaços em branco dos nomes das colunas
df.columns = df.columns.str.strip()

resultados = []

for row in df.itertuples(index=False):
    # Pega a orentação real da BS
    bs_antenna = AntennaBeamformingImt(
        bs_param.get_antenna_parameters(), 
        azimuth=row.azimuth_bs_imt,          # Agora sim!
        elevation=-bs_param.downtilt         # Mantém o tilt mecânico de -6°
    )

    # 1. Sorteio do UE dentro da área de cobertura local da BS
    bs_antenna.reset_beams()
    
    # Geramos candidatos globais até encontrar um que caia no setor local correto
    while True:
        # Sorteia uma direção global qualquer ao redor da BS
        candidate_az_global = random.uniform(0, 360)
        candidate_el_global = random.uniform(0, 180) # Convenção zenital (0 a 180)
        
        # O próprio simulador calcula onde esse candidato cai no mundo LOCAL da antena
        lo_phi, lo_theta = bs_antenna.to_local_coord(candidate_az_global, candidate_el_global)
        
        # Verifica se o candidato atende estritamente às suas restrições locais
        if (-60 <= lo_phi <= 60) and (90 <= lo_theta <= 100):
            ue_azimuth = candidate_az_global
            ue_elevation = candidate_el_global
            break # Encontrou um UE válido, sai do loop
            
    # Adiciona o feixe apontado para o UE validado utilizando as coordenadas globais encontradas
    bs_antenna.add_beam(ue_azimuth, ue_elevation)

    # Para provar no terminal, vamos tirar um "print" textual dos ângulos locais do UE validado
    lo_phi_verif, lo_theta_verif = bs_antenna.to_local_coord(ue_azimuth, ue_elevation)
    print(f"--- VERIFICAÇÃO GEOMÉTRICA ---")
    print(f"BS Global -> Azim: {row.azimuth_bs_imt}°, Elev: {row.elevation_bs_imt}°")
    print(f"UE Global sorteado -> Azim: {ue_azimuth:.2f}°, Elev: {ue_elevation:.2f}°")
    print(f"Prova Real (ÂNGULOS LOCAIS DO UE) -> Azim Local: {lo_phi_verif[0]:.2f}°, Elev Local: {lo_theta_verif[0]:.2f}°")
    print(f"-------------------------------\n")

    # 2. Direção do Alvo (Vítima FS)
    target_az = row.azimuth_ep_fs % 360
    target_el = 90 - row.elevation_ep_fs

    # 3. CÁLCULO DE GANHO IMT (Interferente)
    imt_gain = bs_antenna.calculate_gain(
        phi_vec=np.array([target_az]),
        theta_vec=np.array([target_el]), 
        beams_l=np.array([0]), 
        co_channel=True
    )[0]

    # 4. CÁLCULO DE GANHO FS (Vítima)
    # Nota: Verifique se sua classe FS usa os ângulos absolutos ou o off-axis
    # Se precisar do ângulo off-axis em relação ao apontamento principal da FS:
    off_axis = fs_antenna.calculate_off_axis_angle(
        Az=row.azimuth_ep_fs, 
        b=row.elevation_ep_fs
    )
    
    fs_gain = fs_antenna.calculate_gain(
        off_axis_angle_vec=np.array([off_axis])
    )[0]

    # 5. Cálculo da Interferência Final (Link Budget)
    # POT_BS_IMT deve estar definida anteriormente (ex: 43 dBm)
    interferencia = POT_BS_IMT + imt_gain + fs_gain - row.path_loss_db
    
    resultados.append({
        'scenario': row.scenario,
        'i_dbm': interferencia,
        'imt_gain': imt_gain,
        'fs_gain': fs_gain
    })

# Converter resultados para um novo DataFrame se desejar analisar estatisticamente
df_resultados = pd.DataFrame(resultados)

# 1. SALVAR OS RESULTADOS EM CSV
df_resultados = pd.DataFrame(resultados)
nome_arquivo = 'resultados_interferencia_6G_FS.csv'
df_resultados.to_csv(nome_arquivo, index=False)
print(f"Resultados salvos com sucesso em: {nome_arquivo}")

# 2. CÁLCULO DA CCDF (Complementary Cumulative Distribution Function)
# A CCDF responde: "Qual a probabilidade da interferência ser maior que 'x'?"
interf_values = df_resultados['i_dbm'].sort_values().values
n = len(interf_values)
# Probabilidades de 1 até 1/n
ccdf_prob = 1.0 - np.arange(1, n + 1) / n

# 3. GERAÇÃO DO GRÁFICO
plt.figure(figsize=(10, 6))
plt.step(interf_values, ccdf_prob, where='post', color='blue', linewidth=2)

# Configurações estéticas e técnicas
plt.yscale('log') # Escala logarítmica é padrão para CCDF de interferência
plt.grid(True, which="both", ls="--", alpha=0.7)
plt.title('CCDF da Interferência Agregada (BS IMT -> FS)')
plt.xlabel('Interferência no Receptor Vítima (dBm)')
plt.ylabel('Probabilidade P(I > i)')

# Se você tiver um critério de proteção (ex: -100 dBm), pode plotar uma linha vertical
# plt.axvline(x=-100, color='red', linestyle=':', label='Critério de Proteção')
# plt.legend()

plt.tight_layout()
plt.savefig('ccdf_interferencia.png', dpi=300)
print("Gráfico CCDF gerado: ccdf_interferencia.png")