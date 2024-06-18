from DNA.model import TBHamType, fourier_analysis_1P, fourier_analysis_2P, PARTICLES
from utils import get_conversion
from typing import Any
from DNA.dynamics import MESolverType
import numpy as np

def plot_fourier_1P(ax, tb_ham: TBHamType, init_state: Any, end_state: Any, x_axis: str):
    get_frame_fourier(ax, x_axis)
    average_pop, amplitudes, frequencies = fourier_analysis_1P(tb_ham, init_state, end_state) 
    frequencies = np.array(frequencies) * get_conversion(tb_ham.unit, 'rad/ps')
    if x_axis.lower() == 'frequency':
        ax.plot(amplitudes, frequencies, '.', markersize=12)
    if x_axis.lower() == 'period':
        periods = 2*np.pi/frequencies
        ax.plot(amplitudes, periods, '.', markersize=12)

def plot_fourier_2P(ax, tb_ham: TBHamType, init_state: Any, end_state: Any):
    get_frame_fourier(ax, x_axis)
    for particle in PARTICLES:
        average_pop, amplitudes, frequencies = fourier_analysis_2P(tb_ham, init_state, end_state, particle)
        frequencies = np.array(frequencies) * get_conversion(tb_ham.unit, 'rad/ps')
        if x_axis.lower() == 'frequency':
            ax.plot(amplitudes, frequencies, '.', markersize=12, label = particle)
        if x_axis.lower() == 'period':
            periods = 2*np.pi/frequencies
            ax.plot(amplitudes, periods, '.', markersize=12)
    ax.legend()

def get_frame_fourier(ax, x_axis: str):
    if x_axis.lower() == 'frequency':
        ax.set_xlabel(f"Frequency in rad/ps")
    if x_axis.lower() == 'period':
        ax.set_xlabel(f"Frequency in ps")
    ax.set_ylabel("Amplitude")


# ---------------------------------------------------------------------------------------------------------

def get_frame(ax, dna_sequence):
    ax[0].set_ylabel(r'$P_k^{(\mathrm{acc})}$', fontsize=20)
    for particle_idx, particle in enumerate(PARTICLES):
        ax[particle_idx].plot(J_list, [0]*len(J_list), color='k') # plot the bottom line
        ax[particle_idx].set_xlabel(r'$J\,(\mathrm{meV})$')
        ax[particle_idx].set_title(particle, fontsize=20)
        ax[particle_idx].tick_params(axis='both', which='both')
        ax[particle_idx].set_ylim(0, 1.05)
    
    for dna_base_idx, dna_base in enumerate(dna_sequence):
        ax[0].text(0, dna_base_idx/len(dna_sequence), dna_base, fontsize=25, color=COLORS_DNA_BASES[dna_base], alpha=0.6)

# --------------------------------------------------------------------------------------------------

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