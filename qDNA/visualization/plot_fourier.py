"""This module provides functions to plot the Fourier transform and the average
population for quantum DNA models."""

import numpy as np

from ..utils import get_conversion
from ..tools import DEFAULTS
from . import COLORS_DNA_BASES, PARTICLES, COLORS_PARTICLES

__all__ = [
    "plot_fourier",
    "plot_average_pop",
    "get_frame_average_pop",
]

# -----------------------------------------------------------------------------------------------------


def plot_fourier(ax, tb_ham, init_state, end_state, x_axis):
    """Plots the Fourier transform of the transition amplitudes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    tb_ham : TB_Ham
        The tight-binding Hamiltonian.
    init_state : Any
        The initial state.
    end_state : Any
        The end state.
    x_axis : str
        The x-axis type ('frequency' or 'period').
    """

    # calculation
    amplitudes_dict = tb_ham.get_amplitudes(init_state, end_state)
    frequencies_dict = tb_ham.get_frequencies(init_state, end_state)
    # transform frequencies to rad/ps
    for particle in tb_ham.particles:
        conversion = get_conversion(tb_ham.unit, "rad/ps") / (2 * np.pi)
        frequencies_dict[particle] = np.array(frequencies_dict[particle]) * conversion

    # plotting
    for particle in tb_ham.particles:
        amplitudes = amplitudes_dict[particle]
        frequencies = frequencies_dict[particle]
        # transform frequency to period (in fs)
        periods = 1e3 / frequencies

        # frequency as x-axis
        if x_axis.lower() == "frequency":
            ax.plot(
                frequencies,
                np.abs(amplitudes),
                ".",
                label=particle,
                color=COLORS_PARTICLES[particle],
                markersize=12,
            )

        # period as x-axis
        elif x_axis.lower() == "period":
            ax.plot(
                periods,
                np.abs(amplitudes),
                ".",
                label=particle,
                color=COLORS_PARTICLES[particle],
                markersize=12,
            )

    # plot settings
    if x_axis.lower() == "frequency":
        ax.set_xlabel("Frequency in rad/ps")
    elif x_axis.lower() == "period":
        ax.set_xlabel("Period in fs")
    ax.set_ylabel("Amplitude")
    ax.set_ylim(0.01)
    ax.legend()


# ---------------------------------------------------------------------------------------------------------


def get_cumulative_average_pop(tb_ham, J_list):
    """Computes the cumulative average population.

    Parameters
    ----------
    tb_ham : TB_Ham
        The tight-binding Hamiltonian.
    J_list : list
        List of interaction parameters.

    Returns
    -------
    numpy.ndarray
        Cumulative average population.
    """

    # since we investigate the effect of Coulomb interaction, we are limited to the 2P description
    assert (
        tb_ham.description == "2P"
    ), "this analysis is only available in the 2P description"

    # TODO: make this work for other initial states
    init_state = (
        DEFAULTS["me_kwargs_default"].get("init_e_state"),
        DEFAULTS["me_kwargs_default"].get("init_h_state"),
    )
    num_sites = len(tb_ham.tb_basis)

    # pop_list contains the average population for each particle, J, and tb_site
    pop_list = np.zeros((len(tb_ham.particles), len(J_list), num_sites))

    # calculate the average population for each particle, J, and tb_site using tb_ham.get_average_pop
    for J_idx, J in enumerate(J_list):
        tb_ham.interaction_param = J
        for tb_site_idx, tb_site in enumerate(tb_ham.tb_basis):
            average_pop = tb_ham.get_average_pop(init_state, tb_site)
            for particle_idx, particle in enumerate(tb_ham.particles):
                pop_list[particle_idx][J_idx][tb_site_idx] = average_pop[particle]

    # calculate the cumulative average population
    cumulative_pop_list = [0] * (num_sites + 1)

    # add zero population
    running_pop_list = np.zeros((len(tb_ham.particles), len(J_list)))
    cumulative_pop_list[-1] = np.array(running_pop_list)

    # add cumulative population
    for tb_basis_idx in range(num_sites):
        running_pop_list += pop_list[:, :, tb_basis_idx]
        cumulative_pop_list[tb_basis_idx] = np.array(running_pop_list)

    return np.array(cumulative_pop_list)


def get_frame_average_pop(ax, dna_seq, J_list, J_unit):
    """Sets up the frame for the average population plot.

    Parameters
    ----------
    ax : list of matplotlib.axes.Axes
        The axes to plot on.
    dna_seq : list
        List of DNA bases in the sequence.
    J_list : list
        List of interaction parameters.
    J_unit : str
        The unit of the interaction parameters.
    """

    ax[0].set_ylabel("Acc. Population")
    for particle_idx, particle in enumerate(PARTICLES):
        # plot the bottom line
        ax[particle_idx].plot(J_list, [0] * len(J_list), color="k")
        ax[particle_idx].set_xlabel("J [" + J_unit + "]")
        ax[particle_idx].set_title(particle.capitalize(), fontsize=18)
        # ax[particle_idx].set_ylim(0, 1.05)

    # plot the DNA bases as letters
    for dna_base_idx, dna_base in enumerate(dna_seq):
        ax[0].text(
            0,
            dna_base_idx / len(dna_seq),
            dna_base,
            fontsize=20,
            color="k",
            alpha=0.8,
        )


def plot_average_pop(ax, tb_ham, J_list, J_unit):
    """Plots the average population for each interaction parameter.

    Parameters
    ----------
    ax : list of matplotlib.axes.Axes
        The axes to plot on.
    tb_ham : TB_Ham
        The tight-binding Hamiltonian.
    J_list : list
        List of interaction parameters.
    J_unit : str
        The unit of the interaction parameters.
    """

    # calculation
    J_list = np.array(J_list) * get_conversion(J_unit, tb_ham.unit)
    cumulative_average_pop = get_cumulative_average_pop(tb_ham, J_list)

    # plotting
    J_unit = tb_ham.unit
    dna_seq = tb_ham.tb_sites_flattened

    for particle_idx, _ in enumerate(PARTICLES):
        for dna_base_idx, dna_base in enumerate(dna_seq):
            # black lines
            ax[particle_idx].plot(
                J_list,
                cumulative_average_pop[dna_base_idx][:][particle_idx],
                color="k",
            )

            # fill between the black lines
            ax[particle_idx].fill_between(
                J_list,
                cumulative_average_pop[dna_base_idx - 1][:][particle_idx],
                cumulative_average_pop[dna_base_idx][:][particle_idx],
                color=COLORS_DNA_BASES[dna_base],
                alpha=0.3,
            )

    # plot settings
    get_frame_average_pop(ax, dna_seq, J_list, J_unit)
