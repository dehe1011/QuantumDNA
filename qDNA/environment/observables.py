"""
This module provides functions to calculate various observables.
It includes utilities to construct observables for different particle types (electrons, holes, and excitons) as transitions between states in a given basis,
especially their populations and coherences. The observables are represented in a specified basis of tight-binding site states.
"""

import numpy as np

__all__ = [
    "get_tb_observable",
    "get_eh_observable",
    "get_pop_particle",
    "get_coh_particle",
]

# ------------------------------------------------------------------------------


def get_observable(basis, start_state, end_state):
    """
    Creates a matrix element for the transition between a given start and end state in the provided basis.

    Parameters
    ----------
    basis : List[str]
        The list of basis states.
    start_state : str
        The starting state.
    end_state : str
        The ending state.

    Returns
    -------
    np.ndarray
        The density matrix.

    Raises
    ------
    ValueError
        If the start or end state is not found in the basis list.

    Examples
    --------
    >>> get_observable(['(0, 0)', '(0, 1)', '(0, 2)'], '(0, 1)', '(0, 2)')
    array([[0., 0., 0.],
           [0., 0., 1.],
           [0., 0., 0.]])
    """

    start_state_idx = basis.index(start_state)
    end_state_idx = basis.index(end_state)

    num_basis = len(basis)
    matrix = np.zeros((num_basis, num_basis))
    matrix[start_state_idx, end_state_idx] = 1
    return matrix


def get_tb_observable(tb_basis, start_state, end_state):
    """
    Wrapper function of get_observable. Creates a matrix element for the transition between a given start and end state in the provided basis.

    Parameters
    ----------
    tb_basis : List[str]
        The list of tight-binding site basis states.
    start_state : str
        The starting state.
    end_state : str
        The ending state.

    Returns
    -------
    np.ndarray
        The density matrix.

    Examples
    --------
    >>> get_tb_observable(['(0, 0)', '(0, 1)', '(0, 2)'], '(0, 1)', '(0, 2)')
    array([[0., 0., 0.],
           [0., 0., 1.],
           [0., 0., 0.]])
    """

    return get_observable(tb_basis, start_state, end_state)


def get_eh_observable(tb_basis, particle, start_state, end_state):
    """
    Creates a electron-hole matrix element for the given particle type.

    Parameters
    ----------
    tb_basis : List[str]
        The list of tight-binding site basis states.
    particle : str
        The type of particle ('electron', 'hole', or 'exciton').
    start_state : str
        The starting state.
    end_state : str
        The ending state.

    Returns
    -------
    np.ndarray
        The electron-hole density matrix.

    Raises
    ------
    ValueError
        If the particle type is not recognized.

    Examples
    --------
    >>> get_eh_observable(['(0, 0)', '(1, 0)'], 'electron', '(0, 0)', '(1, 0)')
    array([[0., 0., 1., 0.],
           [0., 0., 0., 1.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]])
    """
    num_sites = len(tb_basis)

    # Create the electron observable
    if particle == "electron":
        eh_observable = np.kron(
            get_observable(tb_basis, start_state, end_state), np.eye(num_sites)
        )

    # Create the hole observable
    if particle == "hole":
        eh_observable = np.kron(
            np.eye(num_sites), get_observable(tb_basis, start_state, end_state)
        )

    # Create the exciton observable
    if particle == "exciton":
        eh_observable = np.kron(
            get_observable(tb_basis, start_state, end_state),
            get_observable(tb_basis, start_state, end_state),
        )
    return eh_observable


def get_pop_particle(tb_basis, particle, state):
    """
    Creates the population observable for a specific particle and state.

    Parameters
    ----------
    tb_basis : List[str]
        The list of tight-binding site basis states.
    particle : str
        The type of particle ('electron', 'hole', or 'exciton').
    state : str
        The state of the particle.

    Returns
    -------
    np.ndarray
        The population density matrix.

    Examples
    --------
    >>> get_pop_particle(['(0, 0)', '(1, 0)'], 'electron', '(0, 0)')
    array([[1., 0., 0., 0.],
           [0., 1., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]])
    """

    return get_eh_observable(tb_basis, particle, state, state)


def get_coh_particle(tb_basis, particle, state1, state2):
    """
    Creates the coherence observable for a specific particle and pair of states.

    Parameters
    ----------
    tb_basis : List[str]
        The list of tight-binding site basis states.
    particle : str
        The type of particle ('electron', 'hole', or 'exciton').
    state1 : str
        The first state.
    state2 : str
        The second state.

    Returns
    -------
    np.ndarray
        The coherence density matrix.

    Examples
    --------
    >>> get_coh_particle(['(0, 0)', '(1, 0)'], 'electron', '(0, 0)', '(1, 0)')
    array([[0., 0., 1., 0.],
           [0., 0., 0., 1.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]])
    """

    return get_eh_observable(tb_basis, particle, state1, state2)
