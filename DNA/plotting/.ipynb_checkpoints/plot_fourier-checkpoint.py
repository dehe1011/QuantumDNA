from DNA.model import TBHamType, fourier_analysis_1P, fourier_analysis_2P, PARTICLES
from utils import get_conversion
from typing import Any
from DNA.dynamics import MESolverType
import numpy as np

def plot_fourier_1P(ax, tb_ham: TBHamType, init_state: Any, end_state: Any, x_axis: str):
    average_pop, amplitudes, frequencies = fourier_analysis_1P(tb_ham, init_state, end_state) 
    frequencies = np.array(frequencies) * get_conversion(tb_ham.unit, 'rad/ps')
    if x_axis.lower() == 'frequency':
        ax.plot(amplitudes, frequencies, '.', markersize=12)
        ax.set_xlabel(f"Frequency in rad/ps")
    if x_axis.lower() == 'period':
        periods = 2*np.pi/frequencies
        ax.plot(amplitudes, periods, '.', markersize=12)
        ax.set_xlabel(f"Period in ps")
    ax.set_ylabel("Amplitude")

def plot_fourier_2P(ax, tb_ham: TBHamType, init_state: Any, end_state: Any):
    for particle in PARTICLES:
        average_pop, amplitudes, frequencies = fourier_analysis_2P(tb_ham, init_state, end_state, particle)
        frequencies = np.array(frequencies) * get_conversion(tb_ham.unit, 'rad/ps')
        if x_axis.lower() == 'frequency':
            ax.plot(amplitudes, frequencies, '.', markersize=12, label = particle)
        if x_axis.lower() == 'period':
            periods = 2*np.pi/frequencies
            ax.plot(amplitudes, periods, '.', markersize=12)
    ax.legend()
    if x_axis.lower() == 'frequency':
        ax.set_xlabel(f"Frequency in rad/ps")
    if x_axis.lower() == 'period':
        ax.set_xlabel(f"Frequency in ps")
    ax.set_ylabel("Amplitude")

def test_fourier(tb_site: str, me_solver: MESolverType):
    compare = []
    end_state = tb_site
    tb_site = basis_converter(tb_site, me_solver.tb_model.tb_basis)
    for particle in PARTICLES:
        average_pop, amplitudes, frequencies = fourier_analysis_2P(me_solver.tb_ham, me_solver.init_state, end_state, particle)
        pop_fourier = [get_pop_fourier(t, average_pop, amplitudes, frequencies) for t in me_solver.times]
        pop_me_solver = [dm[tb_site, tb_site].real for dm in me_solver.get_result_particle(particle)]
        compare.append( np.allclose(pop_me_solver, pop_fourier, atol = 1e-6) )
    return compare