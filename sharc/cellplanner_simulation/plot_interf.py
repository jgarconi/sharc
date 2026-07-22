# -*- coding: utf-8 -*-
"""
Plota INR e interferência a partir de resultado_inr_por_fs_montecarlo.csv
(gerado por simulacao_montecarlo_interferencia_fs.py).
"""
import csv
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

CSV_RESULTADO = "/home/juliana/Documentos/projetos/sharc/sharc/cellplanner_simulation/output/resultado_inr_por_fs_montecarlo.csv"  # <<< AJUSTE >>>
INR_THRESHOLD = -6.0  # <<< AJUSTE >>> critério de proteção (dB), ex: -6 dB (I/N)

# ------------------------------------------------------------------
# 1) Ler resultados
# ------------------------------------------------------------------
data_inr = defaultdict(list)
data_interf = defaultdict(list)

with open(CSV_RESULTADO, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        fs = row["fs_id"]
        data_inr[fs].append(float(row["inr_db"]))
        data_interf[fs].append(float(row["interf_dbm"]))

fs_ids = sorted(data_inr.keys())
num_fs = len(fs_ids)
print(f"{num_fs} FS, {sum(len(v) for v in data_inr.values())} amostras no total")

# ------------------------------------------------------------------
# 2) Plots
# ------------------------------------------------------------------
# fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig, ax = plt.subplots(figsize=(13, 10))

# --- CDF do INR por FS ---
# ax = axes[0, 0]
# for fs in fs_ids:
#     vals = np.sort(data_inr[fs])
#     cdf = np.arange(1, len(vals) + 1) / len(vals)
#     ax.plot(vals, cdf, alpha=0.6, linewidth=1, label=fs)
# ax.axvline(INR_THRESHOLD, color='k', linestyle='--', linewidth=1, label=f"Critério {INR_THRESHOLD} dB")
# ax.set_xlabel("INR [dB]")
# ax.set_ylabel("Probabilidade (CDF)")
# ax.set_title("CDF do INR por estação FS")
# ax.grid(alpha=0.3)
# ax.legend(fontsize=6, ncol=2, loc="lower right")

# --- CDF agregada (todas as FS juntas) ---
# ax = axes[0, 1]
# all_inr = np.concatenate(list(data_inr.values()))
# vals = np.sort(all_inr)
# cdf = np.arange(1, len(vals) + 1) / len(vals)
# ax.plot(vals, cdf, color='crimson')
# ax.axvline(INR_THRESHOLD, color='k', linestyle='--', linewidth=1)
# pct_exceed = 100 * np.mean(all_inr > INR_THRESHOLD)
# ax.set_title(f"CDF agregada do INR (todas as FS)\nP(INR > {INR_THRESHOLD} dB) = {pct_exceed:.2f}%")
# ax.set_xlabel("INR [dB]")
# ax.set_ylabel("Probabilidade (CDF)")
# ax.grid(alpha=0.3)

# --- Boxplot do INR por FS ---
# ax = axes[1, 0]
# ax.boxplot([data_inr[fs] for fs in fs_ids], labels=fs_ids, showfliers=False)
# ax.axhline(INR_THRESHOLD, color='k', linestyle='--', linewidth=1)
# ax.set_ylabel("INR [dB]")
# ax.set_title("Distribuição do INR por estação FS")
# ax.tick_params(axis='x', rotation=90)
# ax.grid(alpha=0.3, axis='y')

# --- CDF da interferência agregada (dBm) ---
# ax = axes[1, 1]
all_interf = np.concatenate(list(data_interf.values()))
vals = np.sort(all_interf)
cdf = np.arange(1, len(vals) + 1) / len(vals)
ax.plot(vals, cdf, color='navy')
ax.set_xlabel("Interferência recebida [dBm]")
ax.set_ylabel("Probabilidade (CDF)")
ax.set_title("CDF da interferência agregada (todas as FS)")
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("inr_interferencia.png", dpi=150)
plt.show()

# ------------------------------------------------------------------
# 3) Resumo por FS
# ------------------------------------------------------------------
print(f"\n{'FS':8s} {'INR médio':>10s} {'INR P95':>10s} {'P(INR>thr)':>12s}")
for fs in fs_ids:
    arr = np.array(data_inr[fs])
    pct = 100 * np.mean(arr > INR_THRESHOLD)
    print(f"{fs:8s} {arr.mean():10.2f} {np.percentile(arr, 95):10.2f} {pct:11.2f}%")