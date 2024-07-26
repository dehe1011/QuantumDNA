"""
This module provides functions to plot population and coherence for quantum DNA models.
"""

import matplotlib.pyplot as plt
import seaborn as sns

from qDNA.tools import get_config
from qDNA.utils import calc_coherence, get_pop_fourier

__all__ = ["plot_pop_fourier", "plot_pop", "plot_pops", "plot_coh", "plot_test_fourier"]

PARTICLES = get_config()["PARTICLES"]
COLORS_PARTICLES = dict(zip(PARTICLES, sns.color_palette()[:3]))


def plot_pop_fourier(ax, tb_ham, init_state, end_state, times, t_unit):
    """
    Plots population using the Fourier decomposition.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    tb_ham : TBHamType
        The tight-binding Hamiltonian.
    init_state : str
        The initial state.
    end_state : str
        The end state.
    times : list of float
        The time points for plotting.
    t_unit : str
        The unit of time.
    """
    tb_ham.unit = "rad/" + t_unit
    amplitudes_dict, frequencies_dict, average_pop_dict = tb_ham.get_fourier(
        init_state, end_state, ["amplitude", "frequency", "average_pop"]
    )
    for particle in tb_ham.particles:
        amplitudes = amplitudes_dict[particle]
        frequencies = frequencies_dict[particle]
        average_pop = average_pop_dict[particle]
        pop_list = [
            get_pop_fourier(t, average_pop, amplitudes, frequencies) for t in times
        ]
        ax.plot(times, pop_list, label=particle, color=COLORS_PARTICLES[particle])

    ax.set_ylabel("Population")
    ax.set_xlabel("Time in " + t_unit)
    ax.legend()


def plot_pop(ax, tb_site, me_solver, add_legend=True):
    """
    Plots population of one base.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    tb_site : str
        The tight-binding site.
    me_solver : MESolverType
        The master equation solver instance.
    add_legend : bool, optional
        Whether to add a legend, by default True.
    """
    tb_basis = me_solver.tb_model.tb_basis
    tb_site_idx = tb_basis.index(tb_site)
    description = me_solver.tb_ham.description
    particles = me_solver.tb_ham.particles
    for particle in particles:
        if description == "2P":
            ax.plot(
                me_solver.times,
                me_solver.get_pop()[particle + "_" + tb_site],
                label=particle,
                color=COLORS_PARTICLES[particle],
            )
        elif description == "1P":
            ax.plot(
                me_solver.times,
                [dm[tb_site_idx, tb_site_idx].real for dm in me_solver.get_result()],
                label=particle,
                color=COLORS_PARTICLES[particle],
            )
    ax.set_ylabel("Population")
    ax.set_xlabel("Time in " + me_solver.t_unit)
    ax.set_ylim(0, 1.02)
    if add_legend:
        ax.legend()


def plot_pops(me_solver):
    """
    Plots populations of all bases.

    Parameters
    ----------
    me_solver : MESolverType
        The master equation solver instance.

    Returns
    -------
    tuple
        A tuple containing the figure and axes.
    """
    num_strands = me_solver.tb_model.num_strands
    num_sites_per_strand = me_solver.tb_model.num_sites_per_strand
    tb_basis = me_solver.tb_model.tb_basis
    fig, axes = plt.subplots(
        num_strands,
        num_sites_per_strand,
        figsize=(5 * num_sites_per_strand, 5 * num_strands),
    )
    for i, (ax, tb_site) in enumerate(zip(axes.flatten(), tb_basis)):
        add_legend = i == 0
        plot_pop(ax, tb_site, me_solver, add_legend=add_legend)
    return fig, axes


def plot_coh(ax, me_solver):
    """
    Plots a measure of coherence as the sum of absolute values of the off-diagonal elements.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    me_solver : MESolverType
        The master equation solver instance.
    """
    description = me_solver.tb_ham.description
    particles = me_solver.tb_ham.particles
    for particle in particles:
        if description == "2P":
            ax.plot(me_solver.times, me_solver.get_coh()[particle], label=particle)
        elif description == "1P":
            ax.plot(
                me_solver.times,
                [calc_coherence(dm.full()) for dm in me_solver.get_result()],
                label=particle,
            )
    ax.set_ylabel("Coherence")
    ax.set_xlabel("Time in " + me_solver.t_unit)
    ax.legend()


def plot_test_fourier(ax, tb_site, me_solver):
    """
    Plots a comparison of the result obtained by Fourier decomposition and the result obtained solving the Schrödinger equation with the qutip mesolver.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    tb_site : str
        The tight-binding site.
    me_solver : MESolverType
        The master equation solver instance.
    """
    plot_pop_fourier(
        ax,
        me_solver.tb_ham,
        me_solver.init_state,
        tb_site,
        me_solver.times,
        me_solver.t_unit,
    )
    plot_pop(ax, tb_site, me_solver)
