import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List 
from DNA.dynamics import MESolverType
from DNA.model import get_pop_fourier, basis_converter, PARTICLES, fourier_analysis_2P, fourier_analysis_1P, TBHamType
from DNA.observables import calc_coherence

COLORS = dict( zip(PARTICLES, sns.color_palette()[:3]) ) 

def plot_pop_fourier(ax, tb_ham: TBHamType, init_state: str, end_state: str, times: List[float]):
    description = tb_ham.particle
    if description == 'exciton':
        for i, particle in enumerate(PARTICLES):
            average_pop, amplitudes, frequencies = fourier_analysis_2P(tb_ham, init_state, end_state, particle)
            pop_list = [get_pop_fourier(t, average_pop, amplitudes, frequencies) for t in times]
            ax.plot(times, pop_list, label = particle, color = COLORS[particle])
        
    if description in ('electron', 'hole'):
        particle = description
        average_pop, amplitudes, frequencies = fourier_analysis_1P(tb_ham, init_state, end_state)
        pop_list = [get_pop_fourier(t, average_pop, amplitudes, frequencies) for t in times]
        ax.plot(times, pop_list, label = particle, color = COLORS[particle])
    ax.set_ylabel("Population")
    ax.set_xlabel("Time in ps")
    ax.set_ylim(0, 1.02)
    ax.legend()

def plot_pop(ax, tb_site: str, me_solver: MESolverType):
    tb_site = basis_converter(tb_site, me_solver.tb_model.tb_basis)
    description = me_solver.tb_ham.particle
    if description == 'exciton':
        for particle in PARTICLES:
            ax.plot( me_solver.times, [dm[tb_site, tb_site].real for dm in me_solver.get_result_particle(particle)], label = particle, color=COLORS[particle])
    if description in ('electron', 'hole'):
        particle = description
        ax.plot( me_solver.times, [dm[tb_site, tb_site].real for dm in me_solver.result], label = particle, color=COLORS[particle])
    ax.set_ylabel("Population")
    ax.set_xlabel("Time in ps")
    ax.set_ylim(0, 1.02)
    ax.legend()
    
def plot_coh(ax, me_solver: MESolverType):
    description = me_solver.tb_ham.particle
    if description == 'exciton':
        for particle in PARTICLES:
            ax.plot( me_solver.times, [calc_coherence(dm) for dm in me_solver.get_result_particle(particle)], label = particle)
    if description in ('electron', 'hole'):
        particle = description
        ax.plot( me_solver.times, [calc_coherence(dm.full()) for dm in me_solver.result], label = particle)
    ax.set_ylabel("Coherence")
    ax.set_xlabel("Time in ps")
    ax.legend()