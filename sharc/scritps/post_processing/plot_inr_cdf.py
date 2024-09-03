"""Plots INR CDF curves from INR samples in an input file.
"""
import argparse
import math
import pandas as pd
import matplotlib.pyplot as plt
from sharc.parameters.constants import BOLTZMANN_CONSTANT


def plot_inr_cdf(samples_file: str, system_bw_mhz: float, system_noise_temp_k: float):
    """Plots INR CDF curves from INR samples in an input file.

    Parameters
    ----------
    samples_file : str
        input file name with INR samples
    system_bw_mhz : float
        system bandwidth in MHz
    system_noise_temp_k : float
        system noise temperature in Kelvin
    """
    system_noise_temp_k = 500

    thermal_noise_db = \
        10 * math.log10(BOLTZMANN_CONSTANT * system_noise_temp_k) + \
        10 * math.log10(system_bw_mhz * 1e6)

    dl_inr_samples = pd.read_csv(samples_file,
                                 skiprows=1,
                                 sep='	',
                                 header=None)

    plt.figure(figsize=(8, 6))

    y = dl_inr_samples[1] + thermal_noise_db

    axs = plt.subplot(1, 1, 1)

    axs.ecdf(y)
    axs.grid()
    plt.title(f"CDF of input file \"{samples_file}\"")
    plt.ylabel('Probability of P < X', fontsize=10, color='black')
    plt.xlabel('Interference [dBW/10MHz]', fontsize=10, color='black')
    plt.show()

if __name__ == "__main__":
    cli_parser = argparse.ArgumentParser(prog="INR CDF plotter",
                                         description="Plots INR CDF from INR samples. \
                                            Usage example: \
                                            python plot_inr_cdf.py -b 100 -t 500 \
                                                ../../output/\[SYS\]\ INR\ samples.txt")

    cli_parser.add_argument("samples_file", help="Input file with INR samples")
    cli_parser.add_argument("-b", "--system_bw", type=float,
                            help="System bandwidth in MHz")
    cli_parser.add_argument("-t", "--noise_temperature", type=float,
                        help="System noise temperature in K")

    args = cli_parser.parse_args()
    plot_inr_cdf(args.samples_file, args.system_bw, args.noise_temperature)
