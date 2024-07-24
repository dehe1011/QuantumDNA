"""
This module provides functions for working with tight-binding (tb) basis and 
electron-hole (eh) basis in quantum systems. It includes functions to generate 
symbolic representations of basis indices, convert between different 
representations, and perform basis changes.
"""

from itertools import product
from ast import literal_eval
import numpy as np

# Shortcuts
# tb: tight-binding
# dim: dimension
# eh: electron-hole

__all__ = [
    "get_tb_basis",
    "get_eh_basis",
    "get_eh_distance",
    "get_particle_eh_states",
    "basis_change",
    "local_to_global",
    "global_to_local",
]

# -------------------------------------------- Tight-binding basis --------------------------------


def get_tb_basis(tb_dims):
    """
    Generates a symbolic string representation of the indices '(i, j)' where i indicates the strand
    number and j the site number inside the strand.

    Parameters
    ----------
    tb_dims: tuple
        A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.

    Returns
    -------
    list
        List of string representations of site indices.

    Examples
    --------
    >>> get_tb_basis((2, 2))
    ['(0, 0)', '(0, 1)', '(1, 0)', '(1, 1)']
    """
    num_sites = tb_dims[0] * tb_dims[1]
    tb_basis = [int_to_str(tb_dims, site_num) for site_num in range(num_sites)]
    return tb_basis


def str_to_int(tb_dims, tb_basis_str):
    """
    Converts a string representation of a site index to an integer representation.

    Parameters
    ----------
    tb_dims: tuple
        A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.
    tb_basis_str: str
        String representation of the site index.

    Returns
    -------
    int
        Integer representation of the site index.

    Example
    -------
    >>> str_to_int((1,3), '(0, 2)')
    2
    """
    return tuple_to_int(tb_dims, str_to_tuple(tb_basis_str))


def int_to_str(tb_dims, tb_basis_int):
    """
    Converts the integer representation of a basis state to the string representation.

    Parameters
    ----------
    tb_dims: tuple
        A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.
    tb_basis_int: int
        Integer representation of the basis state.

    Returns
    -------
    str
        String representation of the basis state.
    """
    return tuple_to_str(int_to_tuple(tb_dims, tb_basis_int))


def str_to_tuple(tb_basis_str):
    """
    Converts the string representation of a basis state to the tuple representation.

    Parameters
    ----------
    tb_basis_str: str
        String representation of the basis state.

    Returns
    -------
    tuple
        Tuple representation of the basis state.
    """
    return literal_eval(tb_basis_str)


def tuple_to_str(tb_basis_tuple):
    """
    Converts the tuple representation of a basis state to the string representation.

    Parameters
    ----------
    tb_basis_tuple: tuple
        Tuple representation of the basis state.

    Returns
    -------
    str
        String representation of the basis state.
    """
    return str(tb_basis_tuple)


def int_to_tuple(tb_dims, tb_basis_int):
    """
    Converts the integer representation of a basis state to the tuple representation.

    Parameters
    ----------
    tb_dims: tuple
        A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.
    tb_basis_int: int
        Integer representation of the basis state.

    Returns
    -------
    tuple
        Tuple representation of the basis state.
    """
    _, num_sites_per_strand = tb_dims
    return (tb_basis_int // num_sites_per_strand, tb_basis_int % num_sites_per_strand)


def tuple_to_int(tb_dims, tb_basis_tuple):
    """
    Converts the tuple representation of a basis state to the integer representation.

    Parameters
    ----------
    tb_dims: tuple
        A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.
    tb_basis_tuple: tuple
        Tuple representation of the basis state.

    Returns
    -------
    int
        Integer representation of the basis state.
    """
    _, num_sites_per_strand = tb_dims
    strand_num, site_in_strand_num = tb_basis_tuple
    return strand_num * num_sites_per_strand + site_in_strand_num


# --------------------------------------- Electron-hole basis ----------------------------------


def get_eh_basis(tb_dims):
    """
    Generates a symbolic string representation of electron-hole basis indices.

    Parameters
    ----------
    tb_dims: tuple
        A tuple (num_strands, num_sites_per_strand) representing the tight-binding dimensions.

    Returns
    -------
    list
        List of tuples representing electron-hole basis indices.

    Examples
    --------
    >>> get_eh_basis((2, 2))
    [('(0, 0)', '(0, 0)'), ('(0, 0)', '(0, 1)'), ('(0, 0)', '(1, 0)'), ('(0, 0)', '(1, 1)'),
     ('(0, 1)', '(0, 0)'), ('(0, 1)', '(0, 1)'), ('(0, 1)', '(1, 0)'), ('(0, 1)', '(1, 1)'),
     ('(1, 0)', '(0, 0)'), ('(1, 0)', '(0, 1)'), ('(1, 0)', '(1, 0)'), ('(1, 0)', '(1, 1)'),
     ('(1, 1)', '(0, 0)'), ('(1, 1)', '(0, 1)'), ('(1, 1)', '(1, 0)'), ('(1, 1)', '(1, 1)')]
    """
    tb_basis = get_tb_basis(tb_dims)
    eh_basis = list(product(tb_basis, tb_basis))
    return eh_basis


def get_eh_distance(eh_basis):
    """
    Calculates the distance between electron and hole for each state in the basis.

    Parameters
    ----------
    eh_basis: list
        List of electron and hole positions as tuples of strings.

    Returns
    -------
    ndarray
        Array of distances between electron and hole for each state.

    Example
    -------
    >>> get_distance([('(0, 0)', '(1, 1)'), ('(1, 0)', '(0, 0)')])
    array([1.41421356, 1.        ])
    """
    distance_list = []
    for position_electron, position_hole in eh_basis:
        position_electron = literal_eval(position_electron)
        position_hole = literal_eval(position_hole)
        distance = np.sqrt(
            sum(
                (idx_electron - idx_hole) ** 2
                for idx_electron, idx_hole in zip(position_electron, position_hole)
            )
        )
        distance_list.append(distance)
    return np.array(distance_list)


def get_particle_eh_states(particle, tb_basis_element, tb_basis):
    """
    Generates a list of electron-hole states for a given particle and site basis element.

    Parameters
    ----------
    particle: str
        The type of particle ('electron', 'hole', or 'exciton').
    tb_basis_element: str
        The site basis element for which to generate the states.
    tb_basis: list
        The list of all site basis elements.

    Returns
    -------
    list
        List of tuples representing the electron-hole states.

    Examples
    --------
    >>> get_particle_eh_states('electron', '(0, 2)', get_tb_basis((1,3)))
    [('(0, 2)', '(0, 0)'), ('(0, 2)', '(0, 1)'), ('(0, 2)', '(0, 2)')]
    """
    if tb_basis_element not in tb_basis:
        raise ValueError(
            f"Given basis element {tb_basis_element} not in tight-binding basis {tb_basis}"
        )

    if particle == "electron":
        return list(product([tb_basis_element], tb_basis))
    if particle == "hole":
        return list(product(tb_basis, [tb_basis_element]))
    if particle == "exciton":
        return list(product([tb_basis_element], [tb_basis_element]))

    raise ValueError(f"Invalid particle type: {particle}")


# -------------------------- Basis change from local to global basis (eigenbasis) ------------------


def basis_change(matrix, states, liouville=False):
    """
    Performs a basis change of the given matrix.

    Parameters
    ----------
    matrix: ndarray
        The matrix to change basis.
    states: ndarray
        The old basis expressed as a vector in the new basis.
    liouville: bool
        Set to True for an open quantum system. The matrix must have dimension N**2 instead of N.

    Returns
    -------
    ndarray
        Matrix in the new basis.
    """
    if liouville:
        states = np.kron(states, np.conj(states))
    return np.matmul(states, np.matmul(matrix, np.conj(states).T))


def global_to_local(matrix, eigs, liouville=False):
    """
    Performs a basis change from the eigenbasis to the site basis.

    Parameters
    ----------
    matrix: ndarray
        The matrix to change basis.
    eigs: ndarray
        The eigenbasis (old) expressed in the site basis (new).
    liouville: bool
        Set to True for an open quantum system. The matrix must have dimension N**2 instead of N.

    Returns
    -------
    ndarray
        Matrix in the site basis.
    """
    return basis_change(matrix, eigs, liouville=liouville)


def local_to_global(matrix, eigs, liouville=False):
    """
    Performs a basis change from the site basis to the eigenbasis.

    Parameters
    ----------
    matrix: ndarray
        The matrix to change basis.
    eigs: ndarray
        The site basis (old) expressed in the eigenbasis (new).
    liouville: bool
        Set to True for an open quantum system. The matrix must have dimension N**2 instead of N.

    Returns
    -------
    ndarray
        Matrix in the eigenbasis.
    """
    return basis_change(matrix, np.conj(eigs).T, liouville=liouville)
