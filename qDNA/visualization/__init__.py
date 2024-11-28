import matplotlib.pyplot as plt
import seaborn as sns

from ..tools import CONFIG

plt.style.use("seaborn-v0_8-paper")

# Update matplotlib default parameters in the mpl-data/matplotlibrc file
plt.rcParams.update(
    {
        # latex rendering
        # "text.usetex": True,
        # "text.latex.preamble": r"\usepackage{amsmath}\usepackage{siunitx}",
        # font
        "font.family": "serif",
        # "font.serif": ["Times New Roman"],  # used in elsarticle (Elsevier)
        "font.serif": ["Computer Modern Roman"],  # used in revtex4-2 (APS)
        "font.size": 15,
        # figure
        # revtex4-2 (APS) one-column: [8.64468/2.54, 8.64468/2.54 * 3/4]
        # elsarticle (Elsevier) one-column: [8.85553/2.54, 8.85553/2.54 * 3/4]
        "figure.figsize": [6.4, 4.8],  # Default: [6.4, 4.8], ratio 4:3
        # axes, labels and legend
        "axes.titlesize": "large",  # Font size for title, default: large
        "axes.labelsize": "medium",  # Font size for axis labels, default: medium
        "legend.fontsize": "small",  # Font size for legend, default: medium
        "xtick.labelsize": "small",  # Font size for x-tick labels, default: medium
        "ytick.labelsize": "small",  # Font size for y-tick labels, default: medium
        # lines and markers
        "lines.linewidth": 2.0,  # Default line width
        "lines.markersize": 8.0,  # Default marker size
        # grid and layout
        "axes.grid": True,  # Always show grid
        "grid.linestyle": "--",  # Grid line style
        "grid.alpha": 0.7,  # Grid transparency
        "figure.autolayout": True,  # Automatically apply tight_layout
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
