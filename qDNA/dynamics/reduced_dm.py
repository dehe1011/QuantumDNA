"""Module for reducing density matrices to the electron, hole or exciton subspace."""

from itertools import product
import numpy as np

from ..hamiltonian import delete_groundstate
from ..environment import get_eh_observable, get_tb_observable

__all__ = ["get_reduced_dm", "get_reduced_dm_eigs"]

# ----------------------------------------------------------------------


def get_reduced_dm(dm, particle, tb_basis):
    """Reduces the density matrix for a specific particle type.

    Parameters
    ----------
    dm : np.ndarray
        The initial density matrix.
    particle : str
        The type of particle ('electron' or 'hole').
    tb_basis : List[str]
        The list of tight-binding site basis states.

    Returns
    -------
    np.ndarray
        The reduced density matrix for the specified particle.

    Raises
    ------
    ValueError
        If the particle type is not recognized.

    Examples
    --------
    >>> dm = np.eye(4)
    >>> tb_basis = ['(0, 0)', '(1, 0)']
    >>> get_reduced_dm(dm, 'electron', tb_basis)
    array([[1.+0.j, 0.+0.j],
           [0.+0.j, 1.+0.j]])
    """

    num_sites = len(tb_basis)

    # Before taking the trace the groundstate needs to be removed
    if dm.shape[0] != num_sites**2:
        dm = delete_groundstate(dm)

    reduced_dm = np.zeros((num_sites, num_sites), dtype=complex)
    for start_state, end_state in product(tb_basis, repeat=2):
        # Calculate expectation using the full density matrix
        observable = get_eh_observable(tb_basis, particle, start_state, end_state)
        value = np.trace(observable @ dm)

        # Multiply expectation value with the corresponding element of the reduced density matrix
        reduced_dm += value * get_tb_observable(tb_basis, start_state, end_state).T
    return reduced_dm


def get_reduced_dm_eigs(eigs, eigenstate_idx, tb_basis, particle):
    """Reduces the density matrix of a selected eigenstate of the Hamiltonian for a
    specific particle type.

    Parameters
    ----------
    eigs : np.ndarray
        The eigenvectors of the Hamiltonian.
    eigenstate_idx : int
        The index of the eigenstate to reduce.
    tb_basis : List[str]
        The list of tight-binding site basis states.
    particle : str
        The type of particle ('electron' or 'hole').

    Returns
    -------
    np.ndarray
        The reduced density matrix for the specified eigenstate and particle type.
    """

    dm = np.outer(eigs[:, eigenstate_idx], eigs[:, eigenstate_idx].conj())
    return get_reduced_dm(dm, particle, tb_basis)


# ----------------------------------------------------------------------
