import matplotlib.pyplot as plt
import seaborn as sns

from ..tools import CONFIG

plt.style.use("seaborn-v0_8")

plt.rcParams.update(
    {
        # "text.usetex": True,
        # "font.family": "serif",
        # "font.serif": ["Times"],  # used in elsarticle (Elsevier)
        # "font.serif": ["Computer Modern Roman"], # used in revtex4-2 (APS)
        # "text.latex.preamble": r"\usepackage{amsmath}\usepackage{siunitx}",
        "legend.fontsize": 20,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "axes.labelsize": 20,
        "axes.titlesize": 20,
        "font.size": 20,
    }
)

DNA_BASES = CONFIG["DNA_BASES"]
PARTICLES = CONFIG["PARTICLES"]
COLORS_DNA_BASES = dict(
    zip(DNA_BASES, sns.color_palette("colorblind", n_colors=len(DNA_BASES)))
)
COLORS_PARTICLES = dict(zip(PARTICLES, sns.color_palette()[:3]))

from .plot_dna_base_frequency import *
from .plot_eigenspectrum import *
from .plot_fourier import *
from .plot_pop import *
