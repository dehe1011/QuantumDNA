"""
This module provides functions to plot eigenenergies and eigenstates for quantum DNA models.
"""

import numpy as np
import seaborn as sns

from qDNA.tools import get_config
from qDNA import get_reduced_dm_eigs, get_conversion

__all__ = ["PARTICLES", "COLORS_PARTICLES", "plot_eigv", "plot_eigs"]

PARTICLES = get_config()["PARTICLES"]
COLORS_PARTICLES = dict(zip(PARTICLES, sns.color_palette()[:3]))

# -------------------------------------------------------------------------------------


def plot_eigv(ax, tb_ham, energy_unit="100meV"):
    """
    Plots the eigenenergies.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on.
    tb_ham : TBHamType
        The tight-binding Hamiltonian.
    energy_unit : str, optional
        The unit of energy to plot, by default '100meV'.
    """
    eigv, _ = tb_ham.get_eigensystem()
    eigv *= get_conversion(tb_ham.unit, energy_unit)
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
    tb_ham : TBHamType
        The tight-binding Hamiltonian.
    eigenstate_idx : int
        The index of the eigenstate to plot.
    """
    _, eigs = tb_ham.get_eigensystem()
    if tb_ham.description == "2P":
        for particle in tb_ham.particles:
            reduced_dm_eigs = get_reduced_dm_eigs(tb_ham, particle, eigenstate_idx)
            eigs_distribution = np.diag(reduced_dm_eigs)
            eigs_prob_distribution = np.multiply(
                eigs_distribution, eigs_distribution.conj()
            ).real
            ax.plot(
                tb_ham.tb_basis,
                eigs_prob_distribution,
                label=particle,
                color=COLORS_PARTICLES[particle],
            )
    elif tb_ham.description == "1P":
        for particle in tb_ham.particles:
            eigs_distribution = np.diag(eigs)
            eigs_prob_distribution = np.multiply(
                eigs_distribution, eigs_distribution.conj()
            ).real
            ax.plot(
                tb_ham.tb_basis,
                eigs_prob_distribution,
                label=particle,
                color=COLORS_PARTICLES[particle],
            )
    ax.set_ylim(0, 1.02)
    ax.set_title(f"Distribution of Eigenstate {eigenstate_idx}")
    ax.legend()
