import matplotlib as mpl

from ..tools import DNA_BASES, PARTICLES

mpl.rcdefaults()
try:
    mpl.style.use("qDNA-default")
except OSError:
    print("Could not load qDNA-default style. Using seaborn-v0_8-paper style instead.")
    mpl.style.use("seaborn-v0_8-paper")

# A, T, G, C, F
COLORS_DNA_BASES = dict(
    zip(DNA_BASES, ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"])
)  # seaborn color palette "Paired"
COLORS_DNA_BASES = dict(
    zip(DNA_BASES, ["#b93540", "#ed7e40", "#4167c7", "#60abcd", "#1f1e1e"])
)  # seaborn color palette "icefire" with 7 colors

COLORS_PARTICLES = dict(
    zip(PARTICLES, ["#4c72b0", "#55a868", "#c44e52"])
)  # seaborn color palette "deep"
COLORS_PARTICLES = dict(
    zip(PARTICLES, ["#4798ce", "#e25f33", "#1f1e1e"])
)  # seaborn color palette "icefire" with 5 colors

from .plot_dna_base_frequency import *
from .plot_eigenspectrum import *
from .plot_fourier import *
from .plot_pop import *
