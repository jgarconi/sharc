# -*- coding: utf-8 -*-
"""
Plot padrão (CDF) a partir de um arquivo de saída do sharc no formato
padrão (uma coluna, cabeçalho "samples"):

    samples
    45.917...
    42.509...
    ...
"""
import numpy as np
import matplotlib.pyplot as plt

ARQUIVO = "/home/juliana/Documentos/projetos/sharc/sharc/cellplanner_simulation/output/sharc_tradicional_2026-07-22_01/system_inr.csv"  # <<< AJUSTE >>> caminho do arquivo de saída do sharc
THRESHOLD = -10.0            # <<< AJUSTE >>> critério de proteção (dB), ou None pra omitir a linha
XLABEL = "INR [dB]"         # <<< AJUSTE >>> troque conforme a métrica do arquivo (interferência, path loss etc.)

# ------------------------------------------------------------------
# 1) Ler o arquivo (pula o cabeçalho "samples")
# ------------------------------------------------------------------
samples = np.loadtxt(ARQUIVO, skiprows=1)

# ------------------------------------------------------------------
# 2) CDF
# ------------------------------------------------------------------
vals = np.sort(samples)
cdf = np.arange(1, len(vals) + 1) / len(vals)

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(vals, cdf, color='crimson', linewidth=2)

titulo = f"CDF — {len(samples)} amostras"
if THRESHOLD is not None:
    pct_exceed = 100 * np.mean(samples > THRESHOLD)
    ax.axvline(THRESHOLD, color='k', linestyle='--', linewidth=1, label=f"Critério {THRESHOLD}")
    ax.legend(loc="lower right")
    titulo += f" | P(valor > {THRESHOLD}) = {pct_exceed:.2f}%"

ax.set_xlabel(XLABEL)
ax.set_ylabel("Probabilidade (CDF)")
ax.set_title(titulo)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("cdf_output.png", dpi=150)
plt.show()

print(f"{len(samples)} amostras | média = {samples.mean():.2f} | mediana = {np.median(samples):.2f}")