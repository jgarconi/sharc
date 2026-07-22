# -*- coding: utf-8 -*-
"""
Interferência agregada da rede IMT (57 células fixas, 19 sites x 3 setores)
sobre 20 estações FS (posições fixas), via Monte Carlo real do
SimulationDownlink (UEs remanejados a cada snapshot, beams reais, power
control real) -- SUBSTITUINDO apenas o path loss (propagation_system) pelo
valor isotrópico pronto (calculado externamente, sem ganho de antena).

Entradas:
  - parameters_FS_8000.yaml         -> parâmetros reais (antenas, banda, etc.)
  - positions_imt_fs.csv            -> posições fixas (57 células IMT + 20 FS)
  - Dados_EP_-_Cell_-_aux_sharc__1_.csv -> path loss isotrópico (EPxx-CellYYYY -> dB)

Saída:
  - resultado_inr_por_fs_montecarlo.csv  -> INR/interferência por snapshot, por FS
"""
import re
import csv
import numpy as np

from sharc.parameters.parameters import Parameters
from sharc.simulation_downlink import SimulationDownlink
from sharc.propagation.propagation import Propagation

YAML_PATH = "/home/juliana/Documentos/projetos/sharc/sharc/cellplanner_simulation/input/parameters_FS_8000.yaml"          # <<< AJUSTE o caminho >>>
CSV_POSICOES = "/home/juliana/Documentos/projetos/sharc/sharc/cellplanner_simulation/input/positions_imt_fs.csv"          # <<< AJUSTE o caminho >>>
CSV_PATHLOSS = "/home/juliana/Documentos/projetos/sharc/sharc/cellplanner_simulation/input/pathloss_azim_elev.csv"        # <<< AJUSTE o caminho >>>
NUM_SNAPSHOTS = 100                              # <<< AJUSTE (ou leia do YAML) >>_


def cell_number(cell_id: str) -> int:
    """Extrai o número da célula de qualquer formatação (Cell0001, Cell00010, ...)."""
    return int(re.search(r"\d+", cell_id).group())


# ------------------------------------------------------------------
# Propagação substituída: devolve o path loss isotrópico já pronto
# (mesma ordem de colunas que a topologia das 57 células)
# ------------------------------------------------------------------
class PropagationFixedIsotropic(Propagation):
    """Devolve um path loss pré-calculado (antenas isotrópicas) em vez de
    rodar um modelo de propagação real (P.452/P.1812/etc)."""

    def __init__(self, random_number_gen, pl_vector_bs_order):
        super().__init__(random_number_gen)
        self.pl_vector = np.asarray(pl_vector_bs_order)

    def get_loss(
        self, params, frequency, station_a, station_b,
        station_a_gains=None, station_b_gains=None,
    ) -> np.array:
        # station_a = sistema (FS, 1 estação) | station_b = IMT (57 células)
        num_a = station_a.num_stations
        num_b = station_b.num_stations
        assert num_b == len(self.pl_vector), (
            f"path loss tem {len(self.pl_vector)} valores, mas topologia "
            f"tem {num_b} células -- confira a ordem/tamanho"
        )
        return np.tile(self.pl_vector, (num_a, 1))


# ------------------------------------------------------------------
# 1) Ler posições fixas (57 células IMT + 20 FS), normalizando os IDs
# ------------------------------------------------------------------
bs_entries, fs_x, fs_y = [], {}, {}
with open(CSV_POSICOES, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["station_type"] == "IMT":
            bs_entries.append({
                "num": cell_number(row["id"]),
                "x": float(row["x_m"]),
                "y": float(row["y_m"]),
            })
        elif row["station_type"] == "FS":
            fs_x[row["id"]] = float(row["x_m"])
            fs_y[row["id"]] = float(row["y_m"])

# ordena pelas células 1..57 -> garante a mesma ordem usada pela topologia
bs_entries.sort(key=lambda e: e["num"])
bs_cell_numbers = [e["num"] for e in bs_entries]
bs_x = np.array([e["x"] for e in bs_entries])
bs_y = np.array([e["y"] for e in bs_entries])
# boresight físico do setor: 30 / 150 / 270, repetido a cada site (3 células)
bs_azimuth = np.tile([30.0, 150.0, 270.0], len(bs_entries) // 3)

fs_ids = sorted(fs_x.keys())
num_bs = len(bs_entries)
num_fs = len(fs_ids)
print(f"{num_bs} células IMT, {num_fs} FS lidas de {CSV_POSICOES}")

# ------------------------------------------------------------------
# 2) Ler path loss isotrópico (EPxx-CellYYYY -> dB), indexado por número
# ------------------------------------------------------------------
pl_lookup = {}  # (fs_id, cell_number) -> path_loss_db
with open(CSV_PATHLOSS, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        link = row["link"]
        ep_id, cell_str = link.split("-")
        pl_lookup[(ep_id, cell_number(cell_str))] = float(row["path_loss_db"])

# ------------------------------------------------------------------
# 3) Parâmetros reais do YAML
# ------------------------------------------------------------------
parameters = Parameters()
parameters.set_file_name(YAML_PATH)
parameters.read_params()

num_snapshots = NUM_SNAPSHOTS  # ou: parameters.general.num_snapshots

resultados = []

for fs_idx, fs_id in enumerate(fs_ids):
    print(f"\n=== {fs_id} ({fs_idx + 1}/{num_fs}) ===")

    # posição fixa desta FS
    parameters.single_earth_station.geometry.location.type = "FIXED"
    parameters.single_earth_station.geometry.location.fixed.x = fs_x[fs_id]
    parameters.single_earth_station.geometry.location.fixed.y = fs_y[fs_id]

    simulation = SimulationDownlink(parameters, YAML_PATH)

    # -- posições fixas das 57 células (sobrescreve o hexágono padrão) --
    simulation.topology.calculate_coordinates()
    simulation.topology.x = bs_x.copy()
    simulation.topology.y = bs_y.copy()
    simulation.topology.azimuth = bs_azimuth.copy()

    # -- path loss isotrópico pronto, na MESMA ordem das 57 células --
    pl_vector = np.array([
        pl_lookup[(fs_id, n)] for n in bs_cell_numbers
    ])
    rand_gen = np.random.RandomState(parameters.general.seed)
    simulation.propagation_system = PropagationFixedIsotropic(rand_gen, pl_vector)

    simulation.initialize()

    seed_gen = np.random.RandomState(parameters.general.seed)
    for snap in range(1, num_snapshots + 1):
        seed = int(seed_gen.randint(1, 2**32 - 1))
        simulation.snapshot(write_to_file=False, snapshot_number=snap, seed=seed)

    inr = np.array(simulation.results.system_inr)
    interf_dbm = np.array(simulation.results.system_dl_interf_power)

    for i in range(len(inr)):
        resultados.append({
            "fs_id": fs_id, "snapshot": i + 1,
            "inr_db": inr[i], "interf_dbm": interf_dbm[i],
        })

    print(
        f"  {len(inr)} amostras | INR médio = {inr.mean():.2f} dB | "
        f"P(INR > -6 dB) = {100 * np.mean(inr > -6):.2f}%",
    )

# ------------------------------------------------------------------
# 4) Exportar
# ------------------------------------------------------------------
with open("resultado_inr_por_fs_montecarlo.csv", "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=["fs_id", "snapshot", "inr_db", "interf_dbm"])
    writer.writeheader()
    writer.writerows(resultados)

print("\nResultado salvo em resultado_inr_por_fs_montecarlo.csv")