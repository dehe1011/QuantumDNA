"""
This module provides functions to plot the Fourier transform and the average population for quantum DNA models.
"""

import seaborn as sns
import numpy as np

from qDNA.tools import get_config
from qDNA.utils import get_conversion

__all__ = [
    "plot_fourier",
    "plot_average_pop",
    "get_frame_fourier",
    "get_frame_average_pop",
]

DNA_BASES = get_config()["DNA_BASES"]
COLORS_DNA_BASES = dict(
    zip(DNA_BASES, sns.color_palette("colorblind", n_colors=len(DNA_BASES)))
)
PARTICLES = get_config()["PARTICLES"]
COLORS_PARTICLES = dict(zip(PARTICLES, sns.color_palette()[:3]))

# -----------------------------------------------------------------------------------------------------


def plot_fourier(ax, tb_ham, init_state, end_state, x_axis):
    """
    Plots the Fourier transform of the transition amplitudes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    tb_ham : TBHamType
        The tight-binding Hamiltonian.
    init_state : Any
        The initial state.
    end_state : Any
        The end state.
    x_axis : str
        The x-axis type ('frequency' or 'period').
    """
    get_frame_fourier(ax, x_axis)
    amplitudes_dict = tb_ham.get_amplitudes(init_state, end_state)
    frequencies_dict = tb_ham.get_frequencies(init_state, end_state)
    for particle in tb_ham.particles:
        amplitudes = amplitudes_dict[particle]
        frequencies = np.array(frequencies_dict[particle]) * get_conversion(
            tb_ham.unit, "rad/ps"
        )
        if x_axis.lower() == "frequency":
            ax.plot(
                frequencies,
                amplitudes,
                ".",
                label=particle,
                color=COLORS_PARTICLES[particle],
                markersize=12,
            )
        elif x_axis.lower() == "period":
            periods = 1e3 / frequencies
            ax.plot(
                periods,
                amplitudes,
                ".",
                label=particle,
                color=COLORS_PARTICLES[particle],
                markersize=12,
            )
    ax.legend()


def get_frame_fourier(ax, x_axis):
    """
    Sets up the frame for the Fourier plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    x_axis : str
        The x-axis type ('frequency' or 'period').
    """
    if x_axis.lower() == "frequency":
        ax.set_xlabel("Frequency in rad/ps")
    elif x_axis.lower() == "period":
        ax.set_xlabel("Period in fs")
    ax.set_ylabel("Amplitude")


# ---------------------------------------------------------------------------------------------------------


def get_cumulative_average_pop(tb_ham, J_list):
    """
    Computes the cumulative average population.

    Parameters
    ----------
    tb_ham : TBHamType
        The tight-binding Hamiltonian.
    J_list : list
        List of interaction parameters.

    Returns
    -------
    np.ndarray
        Cumulative average population.
    """
    assert (
        tb_ham.description == "2P"
    ), "this analysis is only available in the 2P description"
    init_state = (
        get_config()["me_kwargs_default"].get("init_e_state"),
        get_config()["me_kwargs_default"].get("init_h_state"),
    )
    num_sites = len(tb_ham.tb_basis)

    pop_list = np.zeros((len(tb_ham.particles), len(J_list), num_sites))
    for particle_idx, particle in enumerate(tb_ham.particles):
        for J_idx, J in enumerate(J_list):
            tb_ham.interaction_param = J
            for tb_site_idx, tb_site in enumerate(tb_ham.tb_basis):
                pop_list[particle_idx][J_idx][tb_site_idx] = tb_ham.get_average_pop(
                    init_state, tb_site
                )[particle]
    cumulative_pop_list = [0] * (num_sites + 1)
    running_pop_list = np.zeros((len(tb_ham.particles), len(J_list)))
    cumulative_pop_list[-1] = np.array(running_pop_list)
    for tb_basis_idx in range(num_sites):
        running_pop_list += pop_list[:, :, tb_basis_idx]
        cumulative_pop_list[tb_basis_idx] = np.array(running_pop_list)
    return np.array(cumulative_pop_list)


def get_frame_average_pop(ax, dna_seq, J_list, J_unit):
    """
    Sets up the frame for the average population plot.

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
    ax[0].set_ylabel(r"$P_k^{(\mathrm{acc})}$", fontsize=20)
    for particle_idx, particle in enumerate(PARTICLES):
        ax[particle_idx].plot(
            J_list, [0] * len(J_list), color="k"
        )  # plot the bottom line
        ax[particle_idx].set_xlabel("J in " + J_unit)
        ax[particle_idx].set_title(particle, fontsize=20)
        ax[particle_idx].tick_params(axis="both", which="both")
        ax[particle_idx].set_ylim(0, 1.05)

    for dna_base_idx, dna_base in enumerate(dna_seq):
        ax[0].text(
            0,
            dna_base_idx / len(dna_seq),
            dna_base,
            fontsize=25,
            color=COLORS_DNA_BASES[dna_base],
            alpha=0.6,
        )


def plot_average_pop(ax, tb_ham, J_list, J_unit):
    """
    Plots the average population for each interaction parameter.

    Parameters
    ----------
    ax : list of matplotlib.axes.Axes
        The axes to plot on.
    tb_ham : TBHamType
        The tight-binding Hamiltonian.
    J_list : list
        List of interaction parameters.
    J_unit : str
        The unit of the interaction parameters.
    """
    J_list = np.array(J_list) * get_conversion(J_unit, tb_ham.unit)
    J_unit = tb_ham.unit
    dna_seq = tb_ham.tb_sites_flattened
    cumulative_average_pop = get_cumulative_average_pop(tb_ham, J_list)
    for particle_idx, _ in enumerate(PARTICLES):
        for dna_base_idx, dna_base in enumerate(dna_seq):
            ax[particle_idx].plot(
                J_list, cumulative_average_pop[dna_base_idx][:][particle_idx], color="k"
            )
            ax[particle_idx].fill_between(
                J_list,
                cumulative_average_pop[dna_base_idx - 1][:][particle_idx],
                cumulative_average_pop[dna_base_idx][:][particle_idx],
                color=COLORS_DNA_BASES[dna_base],
                alpha=0.15,
            )
    get_frame_average_pop(ax, dna_seq, J_list, J_unit)
