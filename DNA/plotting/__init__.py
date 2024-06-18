import matplotlib.pyplot as plt

plt.style.use('seaborn-v0_8')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

from .plot_eigenspectrum import *
from .plot_pop import *
from .plot_fourier import *

# particle colors
# DNA base colors 