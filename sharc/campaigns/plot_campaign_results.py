import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import matplotlib
print(matplotlib.get_backend())

class CCDFPlotter:

    def __init__(self, cutoff=0.0002):
        self.cutoff = cutoff
        self.datasets = []
        self.criteria = {}

    # ==============================
    # ADD DATA
    # ==============================
    def add_csv(self, path, label, linestyle="-", color=None):
        df = pd.read_csv(path)
        col = df.columns[0]
        self.datasets.append({
            "data": df[col].values,
            "label": label,
            "linestyle": linestyle,
            "color": color
        })

    def add_array(self, data, label, linestyle="-", color=None):
        self.datasets.append({
            "data": data,
            "label": label,
            "linestyle": linestyle,
            "color": color
        })

    # ==============================
    # CRITÉRIOS
    # ==============================
    def add_criterion(self, name, percentage, db_value, linestyle="dashed", color="black"):
        self.criteria[name] = (percentage, db_value, linestyle, color)

    # ==============================
    # CCDF
    # ==============================
    @staticmethod
    def ccdf_from(data):
        data_sorted = np.sort(data)
        n = len(data_sorted)
        y = 1.0 - np.arange(1, n + 1) / n
        return data_sorted, y

    # ==============================
    # CROSSING
    # ==============================
    @staticmethod
    def find_crossing(x, y, y_target):
        for i in range(len(y) - 1):
            if y[i] >= y_target and y[i+1] <= y_target:
                x1, x2 = x[i], x[i+1]
                y1, y2 = y[i], y[i+1]

                if y1 == y2:
                    return x1

                return x1 + (y_target - y1) * (x2 - x1) / (y2 - y1)

        return None

    # ==============================
    # PLOT
    # ==============================
    def plot(self, save_path=None):
        plt.figure(figsize=(8, 6))

        lim_x = [np.inf, -np.inf]

        for ds in self.datasets:
            x, y = self.ccdf_from(ds["data"])

            plt.plot(
                x, y,
                label=ds["label"],
                linestyle=ds["linestyle"],
                color=ds["color"]
            )

            lim_x[0] = min(lim_x[0], np.min(x))
            lim_x[1] = max(lim_x[1], np.max(x))

            # CROSSINGS
            for name, (perc, db_val, _, _) in self.criteria.items():
                x_cross = self.find_crossing(x, y, perc)

                if x_cross is not None:
                    margin = db_val - x_cross
                    print(f"[{ds['label']}] {name} → X ≈ {x_cross:.2f} dB → Margin = {margin:.2f} dB")
                else:
                    print(f"[{ds['label']}] {name} → NÃO cruza")

        # CRITÉRIOS VISUAIS
        for name, (perc, db_val, linestyle, color) in self.criteria.items():
            plt.axhline(y=perc, linestyle=linestyle, color=color, label=name)
            plt.axvline(x=db_val, linestyle=linestyle, color=color)

        plt.yscale("log")
        plt.xlabel("INR [dB]")
        plt.ylabel("P(I > X)")
        plt.grid(True, which="both", linestyle="--", alpha=0.6)

        plt.ylim(self.cutoff, 1)
        plt.legend(fontsize=8)

        if save_path:
            plt.savefig(save_path, dpi=300)
            print(f"[✓] Figura salva em {save_path}")
        else:
            plt.show()

if __name__ == '__main__':
    plotter = CCDFPlotter()

    plotter.add_csv(
        "/home/juliana/Documentos/projetos/sharc/sharc/campaigns/resim_study_A.1_FSS/2MC_files/resim_study_A.1_FSS_DL_ra1rb1.csv",
        "Study A.1 FSS - Ra1Rb1",
        linestyle="-"
        )
        
    plotter.add_csv(
        "/home/juliana/Documentos/projetos/sharc/sharc/campaigns/resim_study_A.1_FSS/2MC_files/resim_study_A.1_FSS_DL_ra2rb1.csv",
        "Study A.1 FSS - Ra2Rb1",
        linestyle="--"
        )

    plotter.add_criterion("[-6 dB, 0.03%]", 0.0003, -6, "--", "gray")
    plotter.add_criterion("[-7 dB, 0.1%]", 0.001, -7, "-.", "gray")
    plotter.add_criterion("[-10.5 dB, 20%]", 0.2, -10.5, ":", "gray")

    plotter.plot()