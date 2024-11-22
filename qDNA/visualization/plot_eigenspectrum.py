"""
This module provides functions to plot eigenenergies and eigenstates for quantum DNA models.
"""

import numpy as np

from ..utils import get_conversion
from ..dynamics import get_reduced_dm_eigs

from . import PARTICLES, COLORS_PARTICLES

__all__ = ["PARTICLES", "COLORS_PARTICLES", "plot_eigv", "plot_eigs"]

# -------------------------------------------------------------------------------------


def plot_eigv(ax, tb_ham, energy_unit="100meV"):
    """
    Plots the eigenenergies.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    tb_ham : TB_Ham
        The tight-binding Hamiltonian.
    energy_unit : str, optional
        The unit of energy to plot, by default '100meV'.
    """

    # calculation
    eigv, _ = tb_ham.get_eigensystem()
    eigv *= get_conversion(tb_ham.unit, energy_unit)

    # plotting
    ax.set_title("Eigenvalues")
    ax.plot(eigv, "o")
    ax.set_ylabel("Energy in " + energy_unit)


def plot_eigs(ax, tb_ham, eigenstate_idx):
    """
    Plots the distribution of eigenstates over tight-binding sites.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    tb_ham : TB_Ham
        The tight-binding Hamiltonian.
    eigenstate_idx : int
        The index of the eigenstate to plot.
    """

    # calculation
    _, eigs = tb_ham.get_eigensystem()
    eigs_prob_distribution = {}

    # 1P description
    if tb_ham.description == "2P":
        for particle in tb_ham.particles:
            # extract eigenstate population for each local state
            reduced_dm_eigs = get_reduced_dm_eigs(tb_ham, particle, eigenstate_idx)
            eigs_prob_distribution[particle] = np.diag(reduced_dm_eigs).real

    # 2P description
    elif tb_ham.description == "1P":
        for particle in tb_ham.particles:
            # extract eigenstate population for each local state
            eigs_distribution = np.diag(eigs)
            eigs_prob_distribution[particle] = np.multiply(
                eigs_distribution, eigs_distribution.conj()
            ).real

    # plotting
    for particle in tb_ham.particles:
        ax.plot(
            tb_ham.tb_basis,
            eigs_prob_distribution[particle],
            label=particle,
            color=COLORS_PARTICLES[particle],
        )
    ax.set_ylim(0, 1.02)
    ax.set_title(f"Distribution of Eigenstate {eigenstate_idx}")
    ax.legend()
