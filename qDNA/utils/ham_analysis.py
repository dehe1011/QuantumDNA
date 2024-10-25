"""
This module provides functions for analyzing Hamiltonians in quantum systems, including calculations of average population,
amplitudes, frequencies, population using Fourier series, and inverse participation ratio (IPR).

Notes
-----
- IPR: inverse participation ratio
- calc: calculate
- pop: population
- eigs: eigenstates
- val: value
- dims: dimensions
- init: initial
"""

from itertools import product

import numpy as np

__all__ = [
    "calc_average_pop",
    "calc_amplitudes",
    "calc_frequencies",
    "get_pop_fourier",
    "calc_ipr_hamiltonian",
]

# ------------------------------------------------


def calc_average_pop(eigs, init_state, end_state):
    """
    Calculates the time-averaged/ static population of the end state when initially in the initial state.

    Parameters
    ----------
    eigs : np.ndarray
        Eigenvector matrix.
    init_state : int
        Initial state in the basis.
    end_state : int
        End state in the basis.

    Returns
    -------
    float
        The static population of the end state.
    """
    average_pop = np.sum(eigs[end_state, :] ** 2 * eigs[init_state, :] ** 2)
    return average_pop


def calc_amplitudes(eigs, init_state, end_state):
    """
    Calculates the amplitudes for transitions between states.

    Parameters
    ----------
    eigs : np.ndarray
        Eigenvector matrix.
    init_state : int
        Initial state in the basis.
    end_state : int
        End state in the basis.

    Returns
    -------
    np.ndarray
        A list of amplitudes for transitions.
    """
    matrix_dimension = eigs.shape[0]
    amplitudes = [
        2
        * eigs[end_state, i]
        * eigs[init_state, i]
        * eigs[end_state, j]
        * eigs[init_state, j]
        for i, j in product(range(matrix_dimension), repeat=2)
        if i < j
    ]
    return np.real(np.array(amplitudes))


def calc_frequencies(eigv):
    """
    Calculates the frequencies corresponding to energy level differences.

    Parameters
    ----------
    eigv : np.ndarray
        Eigenvalue vector.

    Returns
    -------
    np.ndarray
        A list of frequencies.
    """
    matrix_dimension = len(eigv)
    frequencies = [
        abs(eigv[i] - eigv[j]).real
        for i, j in product(range(matrix_dimension), repeat=2)
        if i < j
    ]
    return np.array(frequencies)


def get_pop_fourier(t, average_pop, amplitudes, frequencies):
    """
    Calculates the population using Fourier series.

    Parameters
    ----------
    t : float
        Time variable.
    average_pop : float
        Average population.
    amplitudes : np.ndarray
        Amplitudes for transitions.
    frequencies : np.ndarray
        Frequencies corresponding to energy level differences.

    Returns
    -------
    float
        The population at time t.
    """
    population = average_pop
    for amplitude, frequency in zip(amplitudes, frequencies):
        population += amplitude * np.cos(frequency * t)
    return population


def calc_ipr_hamiltonian(eigs):
    r"""
    Calculates the inverse participation ratio (IPR) for each eigenstate of the Hamiltonian.

    Parameters
    ----------
    eigs : np.ndarray
        Eigenvector matrix.

    Returns
    -------
    list of float
        The IPR values for each eigenstate.

    Notes
    -----
    .. note::

        The IPR is a measure of the localization of an eigenstate. It is defined as the sum of the squared coefficients of the
        eigenvector. The IPR ranges from 1 to N, where N is the dimension of the Hilbert space. The IPR values can be used to
        distinguish between localized and delocalized eigenstates.

        .. math::
            \mathrm{IPR} = 1/\sum_{i=1}^{N} |c_i|^4.

        - Localized state: :math:`\mathrm{IPR} = 1/N`
        - Delocalized coherent state: :math:`\mathrm{IPR} = N`
    """

    dims = eigs.shape[0]
    ratios = []
    for idx in range(dims):
        ipr_val = 1 / np.sum(coeff**4 for coeff in eigs[:, idx])
        ratios.append(ipr_val)
    return ratios
