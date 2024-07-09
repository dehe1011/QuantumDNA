import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List 

from utils import get_config
from DNA import MESolverType, TBHamType, calc_coherence, get_pop_fourier

PARTICLES = get_config()['PARTICLES']
COLORS_PARTICLES = dict( zip(PARTICLES, sns.color_palette()[:3]) ) 

def plot_pop_fourier(ax, tb_ham: TBHamType, init_state: str, end_state: str, times: List[float], t_unit: str):
    """
    Plots population using the Fourier decomposition.
    """
    tb_ham.unit = 'rad/'+t_unit
    amplitudes_dict, frequencies_dict, average_pop_dict = tb_ham.get_fourier(init_state, end_state, ['amplitude', 'frequency', 'average_pop'])
    for particle in tb_ham.particles:
        amplitudes = amplitudes_dict[particle]
        frequencies = frequencies_dict[particle]
        average_pop = average_pop_dict[particle]
        pop_list = [get_pop_fourier(t, average_pop, amplitudes, frequencies) for t in times]
        ax.plot(times, pop_list, label = particle, color = COLORS_PARTICLES[particle])

    ax.set_ylabel("Population")
    ax.set_xlabel("Time in "+t_unit)
    ax.set_ylim(0, 1.02)
    ax.legend()

def plot_pop(ax, tb_site: str, me_solver: MESolverType, add_legend: bool = True):
    """ 
    Plots population of one base.
    """
    tb_basis = me_solver.tb_model.tb_basis
    tb_site_idx = tb_basis.index(tb_site)
    description = me_solver.tb_ham.description
    particles = me_solver.tb_ham.particles
    for particle in particles:
        if description == '2P':
            ax.plot( me_solver.times, me_solver.get_pop()[particle+'_'+tb_site], label=particle, color=COLORS_PARTICLES[particle])
        if description == '1P':
            ax.plot( me_solver.times, [dm[tb_site_idx, tb_site_idx].real for dm in me_solver.get_result()], label=particle, color=COLORS_PARTICLES[particle])
    ax.set_ylabel("Population")
    ax.set_xlabel("Time in "+me_solver.t_unit)
    ax.set_ylim(0, 1.02)
    if add_legend:
        ax.legend()

def plot_pops(me_solver: MESolverType):
    """ 
    Plots populations of all bases.
    """
    num_strands = me_solver.tb_model.num_strands
    num_sites_per_strand = me_solver.tb_model.num_sites_per_strand
    tb_basis = me_solver.tb_model.tb_basis
    fig, axes = plt.subplots(num_strands, num_sites_per_strand, figsize=(5*num_sites_per_strand, 5*num_strands) )
    for i, (ax, tb_site) in enumerate( zip(axes.flatten(), tb_basis) ):
        plot_pop(ax, tb_site, me_solver, add_legend = (i==0))
    # return fig, axes
    
def plot_coh(ax, me_solver: MESolverType):
    """
    Plots a measure of coherence as the sum of absolute values of the off-diagonal elements 
    """
    description = me_solver.tb_ham.description
    particles = me_solver.tb_ham.particles
    for particle in particles:
        if description == '2P':
            ax.plot( me_solver.times, me_solver.get_coh()[particle], label = particle)
        if description == '1P':
            ax.plot( me_solver.times, [calc_coherence(dm.full()) for dm in me_solver.get_result()], label = particle)
    ax.set_ylabel("Coherence")
    ax.set_xlabel("Time in "+me_solver.t_unit)
    ax.legend()

def plot_test_fourier(ax, tb_site: str, me_solver: MESolverType):
    """
    Plots a comparision of the result obtained by Fourier decomposition and the result obtained solving the Schrödinger equation with the qutip mesolver.
    """
    plot_pop_fourier(ax, me_solver.tb_ham, me_solver.init_state, tb_site, me_solver.times, me_solver.t_unit)
    plot_pop(ax, tb_site, me_solver)
    