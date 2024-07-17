"""
Module for reducing density matrices to the electron, hole or exciton subspace.
"""

from itertools import product
import numpy as np

from qDNA.model import delete_groundstate
from qDNA.environment import get_eh_observable, get_tb_observable

__all__ = ["get_reduced_dm", "get_reduced_dm_eigs"]


def get_reduced_dm(dm, particle, tb_basis):
    """
    Reduces the density matrix for a specific particle type.

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

    Example
    -------
    >>> dm = np.eye(4)
    >>> get_reduced_dm(dm, 'electron', ['(0, 0)', '(1, 0)'])
    array([[2., 0.],
           [0., 2.]])
    """
    num_sites = len(tb_basis)
    if dm.shape[0] != num_sites**2:
        dm = delete_groundstate(dm)
    reduced_dm = np.zeros((num_sites, num_sites), dtype=np.complex128)
    for start_state, end_state in product(tb_basis, repeat=2):
        value = np.trace(
            get_eh_observable(tb_basis, particle, start_state, end_state) @ dm
        )
        reduced_dm += value * get_tb_observable(tb_basis, start_state, end_state).T
    return reduced_dm


def get_reduced_dm_eigs(tb_ham, particle, eigenstate_idx):
    """
    Reduces the density matrix using the eigenstates of the Hamiltonian for a specific particle type.

    Parameters
    ----------
    tb_ham : TBHamType
        The Hamiltonian object containing the matrix and site basis.
    particle : str
        The type of particle ('electron' or 'hole').
    eigenstate_idx : int
        The index of the eigenstate.

    Returns
    -------
    np.ndarray
        The reduced density matrix.
    """
    _, eigs = tb_ham.get_eigensystem()
    dm = np.outer(eigs[:, eigenstate_idx], eigs[:, eigenstate_idx].conj())
    return get_reduced_dm(dm, particle, tb_ham.tb_basis)
