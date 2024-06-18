# import os
# import sys
# ROOT_DIR = __file__[:__file__.rfind('Quantum_DNA_1.0')]+ 'Quantum_DNA_1.0'
# if ROOT_DIR not in sys.path:
#     del sys.path[0]
#     sys.path.insert(0, ROOT_DIR)
# print(ROOT_DIR)

# Latin Modern Sans
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

def my_func():
    fig, ax = plt.subplots(2,1)
    return fig, ax