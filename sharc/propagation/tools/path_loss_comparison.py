# -*- coding: utf-8 -*-
"""Compare median path loss (PL50) of theoretical models vs ITU-R P.452 and P.1812.

Scenario
--------
Tx: 10 m height, 6 GHz, isotropic, EIRP 30 dBm, at (-22.931034, -47.096705)
Rx: isotropic, 20 m above ground, at (-22.971095, -47.143359)
"""
import os
import io
import sys
import json
import base64
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pyproj import Geod

# Imports do P.452
from sharc.parameters.parameters_p452 import ParametersP452
from sharc.propagation.propagation_clear_air_452 import PropagationClearAir

# Imports do P.1812
from sharc.parameters.parameters_p1812 import ParametersP1812
from sharc.propagation.propagation_p1812 import PropagationP1812

# --- Scenario Configuration --------------------------------------------------
TX_LAT, TX_LON = -22.931033961157787, -47.09670475303627
RX_LAT, RX_LON = -22.97109534026452, -47.143358503534564
FREQ_MHZ = 6000.0
FREQ_GHZ = FREQ_MHZ / 1000.0
HTX_M = 10.0
HRX_M = 20.0
EIRP_DBM = 30.0
STEP_M = 100.0
N_MC = 150               

# Atmospheric defaults
ATM_PRESSURE_HPA = 1017
AIR_TEMP_K = 293.15
N0 = 352.58
DELTA_N = 60
PROFILE_LEN = 100

def _out():
    out_dir = os.environ.get("CLAUDE_SCRATCH", os.path.join(os.path.dirname(__file__), "_campinas_out"))
    os.makedirs(out_dir, exist_ok=True)
    return out_dir

# --- Theoretical Models -----------------------------------------------------
def fspl(d_km):
    d_km = np.maximum(d_km, 1e-4)
    return 32.45 + 20 * np.log10(d_km) + 20 * np.log10(FREQ_MHZ)

def hata_cost231(d_km, hb=HTX_M, hm=HRX_M, f=FREQ_MHZ, big_city=False):
    d_km = np.maximum(d_km, 1e-3)
    logf = np.log10(f)
    a_hm = (1.1 * logf - 0.7) * hm - (1.56 * logf - 0.8)
    C = 3.0 if big_city else 0.0
    return (46.3 + 33.9 * logf - 13.82 * np.log10(hb) - a_hm
            + (44.9 - 6.55 * np.log10(hb)) * np.log10(d_km) + C)

# --- Factories para Modelos do SHARC ----------------------------------------
def create_p452_model(clutter_loss=False, clutter_type="one_end", is_terrain=False):
    """Fábrica para inicializar as variantes do ITU-R P.452."""
    par = ParametersP452()
    par.percentage_p = 50.0      # PL50
    par.Hte = HTX_M
    par.Hre = HRX_M
    par.tx_lat = TX_LAT
    par.rx_lat = RX_LAT
    par.atmospheric_pressure = ATM_PRESSURE_HPA
    par.air_temperature = AIR_TEMP_K
    par.N0 = N0
    par.delta_N = DELTA_N
    par.Dct = 100.0
    par.Dcr = 100.0
    
    par.clutter_loss = clutter_loss
    par.clutter_type = clutter_type
    par.is_terrain = is_terrain
    
    return PropagationClearAir(np.random.RandomState(1), par)

def create_p1812_model(terrain_profile="flat", clutter_mode="none", is_statistical=False):
    """Fábrica para inicializar as variantes do ITU-R P.1812."""
    par = ParametersP1812()
    par.percentage_p = 1.0  # Mediana no tempo
    par.percentage_l = 50.0  # Mediana no local
    par.Hte = HTX_M
    par.Hre = HRX_M
    par.tx_lat = TX_LAT
    par.rx_lat = RX_LAT
    
    # Atmosfera
    par.atmospheric_pressure = ATM_PRESSURE_HPA
    par.air_temperature = AIR_TEMP_K
    par.N0 = N0
    par.delta_N = DELTA_N
    
    # Configurações do P.1812
    par.terrain_profile = terrain_profile
    par.clutter_mode = clutter_mode
    par.clutter_statistical = is_statistical
    
    # Parâmetros operacionais estruturais
    par.profile_resolution = PROFILE_LEN
    par.srtm_directory = os.path.join(_out(), "srtm")
    par.srtm_auto_download = True
    
    return PropagationP1812(np.random.RandomState(1), par)

def _loss_low(prop, d_km):
    d = np.array([[d_km]])
    f = FREQ_GHZ * np.ones((1, 1))
    ind = np.zeros((1, 1), dtype=bool)
    el = np.zeros((1, 1))
    return float(np.ravel(prop.get_loss(d, f, ind, el, np.array([0.0]), np.array([0.0])))[0])

# --- Main Execution ---------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="PL50 model comparison along a path.")
    ap.add_argument("--tx-lat", type=float, default=TX_LAT)
    ap.add_argument("--tx-lon", type=float, default=TX_LON)
    ap.add_argument("--rx-lat", type=float, default=RX_LAT)
    ap.add_argument("--rx-lon", type=float, default=RX_LON)
    ap.add_argument("--tag", default="", help="suffix for output files and title")
    args = ap.parse_args(argv)
    
    tag = ("_" + args.tag) if args.tag else ""
    title_tag = f" [{args.tag}]" if args.tag else ""
    out = _out()

    # Setup distance array
    geod = Geod(ellps="WGS84")
    _, _, total_m = geod.inv(args.tx_lon, args.tx_lat, args.rx_lon, args.rx_lat)
    n = int(round(total_m / STEP_M)) + 1
    d_km = np.linspace(0.0, total_m / 1000.0, n)
    print(f"Path Tx->Rx: {total_m/1000:.2f} km, {n} samples @ {STEP_M:.0f} m, f={FREQ_MHZ:.0f} MHz")
    
    # Initialize P.452 variants
    p452_smooth = create_p452_model(clutter_loss=False)
    p452_2108_one = create_p452_model(clutter_loss=True, clutter_type="one_end")
    p452_2108_both = create_p452_model(clutter_loss=True, clutter_type="both_ends")
    p452_stat = create_p452_model(is_terrain=True)
    
    # Initialize P.1812 variants
    p1812_flat = create_p1812_model(terrain_profile="flat")
    p1812_stat = create_p1812_model(terrain_profile="statistical", is_statistical=True)
    p1812_statc = create_p1812_model(terrain_profile="statistical", clutter_mode="terrain", is_statistical=True)

    # Initialize results dictionary
    res_keys = ["fspl", "hata", "p452_smooth", "p452_2108_one", "p452_2108_both", "p452_stat", 
                "p1812_flat", "p1812_stat", "p1812_statc"]
    res = {k: np.full(n, np.nan) for k in res_keys}

    # Vectorized Theoretical Models
    res["fspl"] = fspl(d_km)
    res["hata"] = hata_cost231(d_km)

    # Point-by-point evaluation
    for k in range(1, n):
        dk = d_km[k]
        sys.stdout.write(f"\rComputing point {k}/{n-1} ({dk:.2f} km)...")
        sys.stdout.flush()

        # Deterministic / Flat Earth Models
        res["p452_smooth"][k]    = _loss_low(p452_smooth, dk)
        res["p452_2108_one"][k]  = _loss_low(p452_2108_one, dk)
        res["p452_2108_both"][k] = _loss_low(p452_2108_both, dk)
        res["p1812_flat"][k]     = _loss_low(p1812_flat, dk)

        # Statistical Models (Monte Carlo)
        mc_452, mc_1812, mc_1812c = [], [], []
        
        for s in range(N_MC):
            # Sequências aleatórias independentes para cada modelo
            p452_stat.random_number_gen   = np.random.RandomState(1000 * k + s)
            p1812_stat.random_number_gen  = np.random.RandomState(2000 * k + s)
            p1812_statc.random_number_gen = np.random.RandomState(3000 * k + s)
            
            # Workaround: Evita loop infinito para distâncias muito curtas (< 1km)
            if dk < 1.0:
                mc_452.append(res["p452_smooth"][k])
                mc_1812.append(res["p1812_flat"][k])
                mc_1812c.append(res["p1812_flat"][k])
            else:
                mc_452.append(_loss_low(p452_stat, dk))
                mc_1812.append(_loss_low(p1812_stat, dk))
                mc_1812c.append(_loss_low(p1812_statc, dk))
                
        # Extrai a mediana (PL50) das realizações de Monte Carlo
        res["p452_stat"][k]   = np.median(mc_452)
        res["p1812_stat"][k]  = np.median(mc_1812)
        res["p1812_statc"][k] = np.median(mc_1812c)

    print("\nComputation complete.")

    # --- Plotting -----------------------------------------------------------
    styles = [
        ("fspl",           "FSPL (Espaço livre)",                    "#7f8c8d", "--"),
        ("hata",           "Okumura-Hata / COST-231 (extrap.)",      "#9b59b6", "--"),
        ("p452_smooth",    "ITU-R P.452 (Smooth Earth)",             "#16a085", "-"),
        ("p452_2108_one",  "ITU-R P.452 (Clutter 2108 - One End)",   "#29b930", "-"),
        ("p452_2108_both", "ITU-R P.452 (Clutter 2108 - Both Ends)", "#b92965", "-"),
        ("p452_stat",      "ITU-R P.452 (Terreno Estatístico)",      "#f1c40f", "-"),
        ("p1812_flat",     "ITU-R P.1812 (Smooth Earth)",            "#2980b9", "-"),
        ("p1812_stat",     "ITU-R P.1812 (Terreno Estatístico)",     "#e67e22", "-"),
        ("p1812_statc",    "ITU-R P.1812 (Terr. Estat. + Clutter)",  "#c0392b", "-"),
    ]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    for key, lbl, col, ls in styles:
        ax.plot(d_km, res[key], ls, color=col, lw=1.8, label=lbl)
        
    ax.set_xlabel("Distância Tx-Rx (km)")
    ax.set_ylabel("Perda de Propagação Mediana PL50 (dB)")
    ax.set_title(f"Comparação PL50 — Campinas-SP{title_tag}\n"
                 f"{FREQ_MHZ:.0f} MHz, hTx={HTX_M:.0f} m, hRx={HRX_M:.0f} m")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower right", bbox_to_anchor=(1, 0))
    fig.tight_layout()

    # Save outputs
    png = os.path.join(out, f"path_loss_comparison{tag}.png")
    fig.savefig(png, dpi=120)
    print(f"Gráfico salvo em -> {png}")

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    
    with open(os.path.join(out, f"path_loss_comparison{tag}_b64.json"), "w") as fh:
        json.dump({"img": base64.b64encode(buf.getvalue()).decode("ascii")}, fh)

    # Endpoint summary
    print(f"\nResumo no Rx (d={d_km[-1]:.2f} km):")
    for key, lbl, _, _ in styles:
        print(f"  {lbl:<40} {res[key][-1]:7.1f} dB   (Prx = {EIRP_DBM-res[key][-1]:7.1f} dBm)")

    np.savez(os.path.join(out, f"path_loss_comparison{tag}.npz"), d_km=d_km, **res)

if __name__ == "__main__":
    main(sys.argv[1:])