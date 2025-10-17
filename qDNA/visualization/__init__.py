import matplotlib as mpl
from .visualization import *
from .plot_base_frequency import *

mpl.rcdefaults()
try:
    mpl.style.use("qDNA-default")
except OSError:
    print("Could not load qDNA-default style. Using seaborn-v0_8-paper style instead.")
    mpl.style.use("seaborn-v0_8-paper")
