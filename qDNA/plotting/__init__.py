import matplotlib.pyplot as plt

plt.style.use("seaborn-v0_8")
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "DejaVu Sans"
plt.rcParams["axes.labelsize"] = 15
plt.rcParams["legend.fontsize"] = 12
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12

import seaborn as sns
from qDNA.tools import get_config

DNA_BASES = get_config()["DNA_BASES"]
PARTICLES = get_config()["PARTICLES"]
COLORS_DNA_BASES = dict(
    zip(DNA_BASES, sns.color_palette("colorblind", n_colors=len(DNA_BASES)))
)
COLORS_PARTICLES = dict(zip(PARTICLES, sns.color_palette()[:3]))

from .plot_eigenspectrum import *
from .plot_pop import *
from .plot_fourier import *
from .plot_dna_base_frequency import *
