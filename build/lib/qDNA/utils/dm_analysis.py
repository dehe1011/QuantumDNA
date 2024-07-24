"""
Module for calculating properties of density matrices (quantum states).
"""

import numpy as np

__all__ = ["calc_trace_distance", "calc_purity", "calc_coherence", "calc_ipr_dm"]

# ----------------------- properties of density matrices (quantum states) ----------------------------


def calc_trace_distance(dm_1, dm_2):
    """
    Calculates the trace distance of two density matrices.

    Parameters
    ----------
    dm_1 : np.ndarray
        First density matrix.
    dm_2 : np.ndarray
        Second density matrix.

    Returns
    -------
    float
        The trace distance between the two density matrices.
    """
    if dm_1.shape != dm_2.shape:
        raise ValueError("Density matrices must have the same shape.")
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
    """
    Calculates the inverse participation ratio (IPR) of a density matrix.

    IPR = :math:`( \sum_{i,j} \rho_{i,j} )^2 / ( N \sum_{i,j} \rho_{i,j}^2 )` .
    Localized state: IPR = 1/N
    Delocalized coherent state: IPR = N
    Maximally mixed state (1/N on diagonal, no off-diagonal elements): IPR = 1

    Parameters
    ----------
    dm : np.ndarray
        Density matrix.

    Returns
    -------
    float
        The inverse participation ratio of the density matrix.
    """
    dims = dm.shape[0]
    numerator = sum(abs(element) for element in dm.flatten())
    denominator = sum(abs(element) ** 2 for element in dm.flatten())
    return numerator**2 / (dims * denominator)
