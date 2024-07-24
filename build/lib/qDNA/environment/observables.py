"""
Module calculates observables. 
"""

import numpy as np

__all__ = [
    "get_tb_observable",
    "get_eh_observable",
    "get_pop_particle",
    "get_coh_particle",
]


def get_observable(basis, start_state, end_state):
    """
    Creates a density matrix for the given start and end states in the provided basis.

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

    Example
    -------
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
    Wrapper function to get the density matrix for tight-binding sites.

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

    Example
    -------
    >>> get_tb_observable(['(0, 0)', '(0, 1)', '(0, 2)'], '(0, 1)', '(0, 2)')
    array([[0., 0., 0.],
           [0., 0., 1.],
           [0., 0., 0.]])
    """
    return get_observable(tb_basis, start_state, end_state)


def get_eh_observable(tb_basis, particle, start_state, end_state):
    """
    Constructs the electron-hole density matrix for the given particle type.

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

    Example
    -------
    >>> get_eh_observable(['(0, 0)', '(1, 0)'], 'electron', '(0, 0)', '(1, 0)')
    array([[0., 0., 1., 0.],
           [0., 0., 0., 1.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]])
    """
    num_sites = len(tb_basis)
    if particle == "electron":
        eh_observable = np.kron(
            get_observable(tb_basis, start_state, end_state), np.eye(num_sites)
        )
    if particle == "hole":
        eh_observable = np.kron(
            np.eye(num_sites), get_observable(tb_basis, start_state, end_state)
        )
    if particle == "exciton":
        eh_observable = np.kron(
            get_observable(tb_basis, start_state, end_state),
            get_observable(tb_basis, start_state, end_state),
        )
    return eh_observable


def get_pop_particle(tb_basis, particle, state):
    """
    Gets the population density matrix for a specific particle and state.

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

    Example
    -------
    >>> get_pop_particle(['(0, 0)', '(1, 0)'], 'electron', '(0, 0)')
    array([[1., 0., 0., 0.],
           [0., 1., 0., 0.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]])
    """
    return get_eh_observable(tb_basis, particle, state, state)


def get_coh_particle(tb_basis, particle, state1, state2):
    """
    Gets the coherence density matrix for a specific particle and pair of states.

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

    Example
    -------
    >>> get_coh_particle(['(0, 0)', '(1, 0)'], 'electron', '(0, 0)', '(1, 0)')
    array([[0., 0., 1., 0.],
           [0., 0., 0., 1.],
           [0., 0., 0., 0.],
           [0., 0., 0., 0.]])
    """
    return get_eh_observable(tb_basis, particle, state1, state2)
