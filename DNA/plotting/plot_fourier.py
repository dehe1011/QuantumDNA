import numpy as np
from typing import Any, List
import seaborn as sns

from utils import get_conversion, get_config
from DNA.model import TBHamType
from DNA.dynamics import MESolverType

DNA_BASES = get_config()['DNA_BASES']
COLORS_DNA_BASES = dict(zip( DNA_BASES, sns.color_palette('colorblind', n_colors=len(DNA_BASES) )))
PARTICLES = get_config()['PARTICLES']
COLORS_PARTICLES = dict( zip(PARTICLES, sns.color_palette()[:3]) ) 

# -----------------------------------------------------------------------------------------------------

def plot_fourier(ax, tb_ham: TBHamType, init_state: Any, end_state: Any, x_axis: str):
    get_frame_fourier(ax, x_axis)
    amplitudes_dict = tb_ham.get_amplitudes(init_state, end_state) 
    frequencies_dict = tb_ham.get_frequencies(init_state, end_state) 
    for particle in tb_ham.particles: 
        amplitudes = amplitudes_dict[particle] 
        frequencies = np.array(frequencies_dict[particle]) * get_conversion(tb_ham.unit, 'rad/ps')
        if x_axis.lower() == 'frequency':
            ax.plot(frequencies, amplitudes, '.',label=particle, color=COLORS_PARTICLES[particle], markersize=12)
        if x_axis.lower() == 'period':
            periods = 1e3/frequencies
            ax.plot(periods, amplitudes, '.', label=particle, color=COLORS_PARTICLES[particle], markersize=12)
    ax.legend()

def get_frame_fourier(ax, x_axis: str):
    if x_axis.lower() == 'frequency':
        ax.set_xlabel(f"Frequency in rad/ps")
    if x_axis.lower() == 'period':
        ax.set_xlabel(f"Period in fs")
    ax.set_ylabel("Amplitude")

# ---------------------------------------------------------------------------------------------------------

def get_cumulative_average_pop(tb_ham: TBHamType, J_list: List): 
    """
    Returns:
        cumulative_average_pop[dna_base_idx][J_idx][particle_idx]
    """
    assert tb_ham.description == '2P', "this analysis is only available in the 2P description" 
    init_state = (get_config()['me_kwargs_default'].get('init_e_state'), get_config()['me_kwargs_default'].get('init_h_state') )
    num_sites = len(tb_ham.tb_basis)
    
    pop_list = np.zeros( (len(tb_ham.particles), len(J_list), num_sites ))
    for particle_idx, particle in enumerate(tb_ham.particles):
        for J_idx, J in enumerate(J_list):
            tb_ham.interaction_param = J
            for tb_site_idx, tb_site in enumerate(tb_ham.tb_basis):
                pop_list[particle_idx][J_idx][tb_site_idx] = tb_ham.get_average_pop(init_state, tb_site)[particle]
    cumulative_pop_list = [0] * ( num_sites+1 )
    running_pop_list = np.zeros( (len(tb_ham.particles), len(J_list)) )
    cumulative_pop_list[-1] = np.array(running_pop_list)
    for tb_basis_idx in range(num_sites):
        running_pop_list += pop_list[:,:,tb_basis_idx]
        cumulative_pop_list[tb_basis_idx] = np.array(running_pop_list)
    return np.array(cumulative_pop_list)

def get_frame_average_pop(ax, dna_seq, J_list, J_unit):
    ax[0].set_ylabel(r'$P_k^{(\mathrm{acc})}$', fontsize=20)
    PARTICLES = get_config()['PARTICLES']
    for particle_idx, particle in enumerate(PARTICLES):
        ax[particle_idx].plot(J_list, [0]*len(J_list), color='k') # plot the bottom line
        ax[particle_idx].set_xlabel('J in '+J_unit)
        ax[particle_idx].set_title(particle, fontsize=20)
        ax[particle_idx].tick_params(axis='both', which='both')
        ax[particle_idx].set_ylim(0, 1.05)
    
    for dna_base_idx, dna_base in enumerate(dna_seq):
        ax[0].text(0, dna_base_idx/len(dna_seq), dna_base, fontsize=25, color=COLORS_DNA_BASES[dna_base], alpha=0.6)

def plot_average_pop(ax, tb_ham: TBHamType, J_list, J_unit):
    J_list *= get_conversion(J_unit, tb_ham.unit)
    J_unit = tb_ham.unit
    dna_seq = tb_ham.tb_sites_flattened
    cumulative_average_pop = get_cumulative_average_pop(tb_ham, J_list)
    PARTICLES = get_config()['PARTICLES']
    for particle_idx, particle in enumerate(PARTICLES):
        for dna_base_idx, dna_base in enumerate(dna_seq):
            ax[particle_idx].plot( J_list, cumulative_average_pop[dna_base_idx][:][particle_idx], color='k')
            ax[particle_idx].fill_between(J_list, cumulative_average_pop[dna_base_idx-1][:][particle_idx], cumulative_average_pop[dna_base_idx][:][particle_idx], color=COLORS_DNA_BASES[dna_base], alpha=0.15)
    get_frame_average_pop(ax, dna_seq, J_list, J_unit)
