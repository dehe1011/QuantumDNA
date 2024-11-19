import matplotlib.pyplot as plt
import seaborn as sns

from ..tools import CONFIG

plt.style.use("seaborn-v0_8-paper")

# Update matplotlib default parameters
plt.rcParams.update(
    {
        # "text.usetex": True,
        # "text.latex.preamble": r"\usepackage{amsmath}\usepackage{siunitx}",
        "font.family": "serif",
        # "font.serif": ["Times"],  # used in elsarticle (Elsevier)
        "font.serif": ["Computer Modern Roman"],  # used in revtex4-2 (APS)
        "font.size": 15,
        "axes.grid": True,  # Always show grid
        "grid.linestyle": "--",  # Grid line style
        "grid.alpha": 0.7,  # Grid transparency
        "axes.labelsize": 15,  # Font size for axis labels
        "legend.fontsize": 12,  # Font size for legend
        "xtick.labelsize": 12,  # Font size for x-tick labels
        "ytick.labelsize": 12,  # Font size for y-tick labels
        "figure.autolayout": True,  # Automatically apply tight_layout
        "lines.linewidth": 2.0,  # Default line width
        "lines.markersize": 8.0,  # Default marker size
    }
)

DNA_BASES = CONFIG["DNA_BASES"]
PARTICLES = CONFIG["PARTICLES"]
COLORS_DNA_BASES = dict(
    zip(DNA_BASES, sns.color_palette("Paired", n_colors=len(DNA_BASES)))
)
COLORS_PARTICLES = dict(zip(PARTICLES, sns.color_palette()[:3]))

from .plot_dna_base_frequency import *
from .plot_eigenspectrum import *
from .plot_fourier import *
from .plot_pop import *
