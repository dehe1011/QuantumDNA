"""
This module provides utility functions for analyzing density matrices.
The functions include calculations for trace distance, purity, coherence, and inverse participation ratio (IPR).
"""

import numpy as np

__all__ = ["calc_trace_distance", "calc_purity", "calc_coherence", "calc_ipr_dm"]

# ------------------------------------------------


def calc_trace_distance(dm_1, dm_2):
    """
    Calculate the trace distance between two density matrices.
    The trace distance is a measure of the distinguishability between two quantum states
    represented by their density matrices. It is defined as half the trace of the absolute
    difference between the two density matrices.

    Parameters
    ----------
    dm_1 : numpy.ndarray
        The first density matrix.
    dm_2 : numpy.ndarray
        The second density matrix.

    Returns
    -------
    float
        The trace distance between the two density matrices.

    Raises
    ------
    ValueError
        If the density matrices do not have the same shape.
    """

    # Check if the density matrices have the same shape
    if dm_1.shape != dm_2.shape:
        raise ValueError("Density matrices must have the same shape.")

    # Calculate the trace distance
    eigenvals = np.linalg.eig(dm_1 - dm_2)[0]
    diag = np.diag(np.absolute(eigenvals))
    return 0.5 * np.trace(diag)


def calc_purity(dm):
    """
    Calculates the purity of a density matrix.

    Parameters
    ----------
    dm : np.ndarray
        Density matrix.

    Returns
    -------
    float
        The purity of the density matrix.
    """

    return np.real(np.trace(np.matmul(dm, dm)))


def calc_coherence(dm):
    """
    Calculates the coherence (absolute sum of the off-diagonals) of a density matrix.

    Parameters
    ----------
    dm : np.ndarray
        Density matrix.

    Returns
    -------
    float
        The coherence of the density matrix.
    """

    coherence = sum(abs(element) for element in dm.flatten()) - np.sum(np.diag(dm))
    return coherence.real


def calc_ipr_dm(dm):
    r"""
    Calculate the inverse participation ratio (IPR) of a density matrix.

    Parameters
    ----------
    dm : np.ndarray
        The density matrix for which the IPR is calculated. Should be a square, complex-valued matrix.

    Returns
    -------
    float
        The inverse participation ratio of the density matrix, providing an indication of state localization.

    Notes
    -----
    .. note::

        The IPR provides a measure of the localization of a quantum state

        .. math::

            \mathrm{IPR} = \frac{\left( \sum_{i,j} |\rho_{i,j}| \right)^2}{N \sum_{i,j} |\rho_{i,j}|^2}

        - Localized state: :math:`IPR = 1/N`, where `N` is the dimension of the density matrix.
        - Delocalized (coherent) state: :math:`IPR = N`, indicating full coherence.
        - Maximally mixed state: :math:`IPR = 1`, representing maximal mixing.

    """

    dims = dm.shape[0]
    numerator = sum(abs(element) for element in dm.flatten())
    denominator = sum(abs(element) ** 2 for element in dm.flatten())
    return numerator**2 / (dims * denominator)
