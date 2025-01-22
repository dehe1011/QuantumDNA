"""This module provides functions to plot population and coherence for quantum DNA
models."""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from qDNA.utils import calc_coherence, get_pop_fourier
from . import COLORS_PARTICLES

__all__ = [
    "plot_pops_heatmap",
    "plot_pop_fourier",
    "plot_pop",
    "plot_pops",
    "plot_coh",
    "plot_test_fourier",
]

# -------------------------------------------------------------------------------------


def plot_pops_heatmap(me_solver, vmax_list=[1, 1, 1], heatmap_type="seaborn"):
    """
    Plots a heatmap of populations for each particle in the system.

    Parameters
    ----------
    me_solver : object
        An instance of a solver that contains the time evolution data and Hamiltonian information.
    vmax_list : list of float, optional
        A list of maximum values for the color scale of each particle's heatmap. Default is [1, 1, 1].
    heatmap_type : str, optional
        The type of heatmap to plot. Options are "seaborn" or "matplotlib". Default is "seaborn".

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the heatmaps.
    ax : numpy.ndarray of matplotlib.axes._subplots.AxesSubplot
        An array of axes objects for each particle's heatmap.

    Notes
    -----
    .. note::

        - The function uses seaborn or matplotlib to generate heatmaps.
        - The color maps used are "Blues", "Reds", and "Greens" for different particles.
    """

    num_particles = len(me_solver.tb_ham.particles)
    fig, ax = plt.subplots(
        1, num_particles, figsize=(6.4 * num_particles, 4.8), sharey=False
    )

    pop_dict = me_solver.get_pop()
    cmaps = ["Blues", "Reds", "Greys"]  # grey scale for exciton
    cmaps = ["Blues", "Reds", "Greens"]  # green scale for exciton
    for i, particle in enumerate(me_solver.tb_ham.particles):
        particle_pop = np.array(
            [value for key, value in pop_dict.items() if key.startswith(particle)]
        )

        # seaborn heatmap (looks prettier in my opinion)
        if heatmap_type == "seaborn":
            heatmap = sns.heatmap(
                particle_pop,
                xticklabels=[],
                yticklabels=[],
                cmap=cmaps[i],
                ax=ax[i],
                cbar=False,
                vmax=vmax_list[i],
            )
            heatmap.figure.colorbar(heatmap.collections[0], ax=ax[i])

        # matplotlib heatmap
        if heatmap_type == "matplotlib":
            im = ax[i].imshow(
                particle_pop, cmap=cmaps[i], aspect="auto", vmax=vmax_list[i]
            )
            fig.colorbar(im, ax=ax[i])

        ax[i].set_xticks([])
        ax[i].set_yticks([])

        y_len, x_len = particle_pop.shape
        ax[i].set_xlabel(f"Time [{me_solver.t_unit}]")
        ax[i].set_title(particle.capitalize())
        ax[i].set_xticks(
            np.linspace(0, x_len, 7), np.round(np.linspace(0, me_solver.t_end, 7), 0)
        )
        ax[i].set_yticks(
            np.arange(y_len) + 0.5, labels=me_solver.tb_ham.tb_sites_flattened
        )

    ax[0].set_ylabel("DNA Bases")
    return fig, ax


def plot_pop_fourier(ax, tb_ham, init_state, end_state, times, t_unit):
    """Plots population using the Fourier decomposition.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    tb_ham : TB_Ham
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

    # calculation
    tb_ham.unit = "rad/" + t_unit
    amplitudes_dict, frequencies_dict, average_pop_dict = tb_ham.get_fourier(
        init_state, end_state, ["amplitude", "frequency", "average_pop"]
    )
    pop_dict = {}
    for particle in tb_ham.particles:
        amplitudes = amplitudes_dict[particle]
        frequencies = frequencies_dict[particle]
        average_pop = average_pop_dict[particle]
        pop_dict[particle] = [
            get_pop_fourier(t, average_pop, amplitudes, frequencies) for t in times
        ]

    # plotting
    for particle in tb_ham.particles:
        ax.plot(
            times, pop_dict[particle], label=particle, color=COLORS_PARTICLES[particle]
        )

    # plot settings
    ax.set_ylabel("Population")
    ax.set_xlabel("Time [" + t_unit + "]")
    ax.legend()
    return ax


def plot_pop(ax, tb_site, me_solver, add_legend=True):
    """Plots population of one base.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    tb_site : str
        The tight-binding site.
    me_solver : ME_Solver
        The master equation solver instance.
    add_legend : bool, optional
        Whether to add a legend, by default True.
    """

    # calculation
    description = me_solver.tb_ham.description
    if description == "2P":
        me_solver.get_pop()
    elif description == "1P":
        me_solver.get_result()

    # plotting
    particles = me_solver.tb_ham.particles
    for particle in particles:
        if description == "2P":
            ax.plot(
                me_solver.times,
                me_solver.pop[particle + "_" + tb_site],
                color=COLORS_PARTICLES[particle],
            )
        elif description == "1P":
            tb_basis = me_solver.tb_model.tb_basis
            tb_site_idx = tb_basis.index(tb_site)
            ax.plot(
                me_solver.times,
                [dm[tb_site_idx, tb_site_idx].real for dm in me_solver.result],
                color=COLORS_PARTICLES[particle],
            )

    # plot settings
    ax.set_ylim(0, 1.02)
    if add_legend:
        ax.set_ylabel("Population")
        ax.set_xlabel("Time [" + me_solver.t_unit + "]")
        ax.legend(particles)
    return ax


def plot_pops(me_solver):
    """Plots populations of all bases.

    Parameters
    ----------
    me_solver : ME_Solver
        The master equation solver instance.

    Returns
    -------
    tuple
        A tuple containing the figure and axes.
    """

    num_strands = me_solver.tb_model.num_strands
    num_sites_per_strand = me_solver.tb_model.num_sites_per_strand
    tb_basis = me_solver.tb_model.tb_basis

    # plotting
    fig, axes = plt.subplots(
        num_strands,
        num_sites_per_strand,
        figsize=(5 * num_sites_per_strand, 5 * num_strands),
        sharex=True,
        sharey=True,
    )
    axes = axes.flatten()
    for i, tb_site in enumerate(tb_basis):
        add_legend = i == 0
        axes[i] = plot_pop(axes[i], tb_site, me_solver, add_legend=add_legend)
    return fig, axes


def plot_coh(ax, me_solver):
    """Plots a measure of coherence as the sum of absolute values of the off- diagonal
    elements.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    me_solver : ME_Solver
        The master equation solver instance.
    """

    # calculation
    description = me_solver.tb_ham.description
    if description == "2P":
        me_solver.get_coh()
    elif description == "1P":
        me_solver.get_result()

    # plotting
    particles = me_solver.tb_ham.particles
    for particle in particles:
        if description == "2P":
            ax.plot(
                me_solver.times,
                me_solver.coh[particle],
                label=particle,
                color=COLORS_PARTICLES[particle],
            )
        elif description == "1P":
            ax.plot(
                me_solver.times,
                [calc_coherence(dm.full()) for dm in me_solver.result],
                label=particle,
                color=COLORS_PARTICLES[particle],
            )

    # plot settings
    ax.set_ylabel("Coherence")
    ax.set_xlabel("Time [" + me_solver.t_unit + "]")
    ax.legend()
    return ax


def plot_test_fourier(ax, tb_site, me_solver):
    """Plots a comparison of the result obtained by Fourier decomposition and the result
    obtained solving the Schrödinger equation with the qutip mesolver.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    tb_site : str
        The tight-binding site.
    me_solver : ME_Solver
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
    return ax
