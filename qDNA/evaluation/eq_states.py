"""This module provides functions to compute equilibrium states for quantum DNA
models."""

import numpy as np
import scipy.constants as c

from ..utils import get_conversion
from ..model import global_to_local, local_to_global
from ..hamiltonian import delete_groundstate

__all__ = ["get_therm_eq_state", "get_deph_eq_state"]

# ---------------------------------------------------------------------


def get_therm_eq_state(me_solver):
    """Calculate the thermal equilibrium state.

    Parameters
    ----------
    me_solver : object
        An instance of a master equation solver which contains the tight-binding Hamiltonian (`tb_ham`)
        and the Lindblad dissipation parameters (`lindblad_diss`).

    Returns
    -------
    numpy.ndarray
        The thermal equilibrium state of the system. If the temperature is zero, returns the ground state.
        Otherwise, returns the thermal equilibrium state as a density matrix in the local basis.

    Notes
    -----
    - The function first checks if the temperature is zero. If so, it returns the ground state.
    - For non-zero temperatures, it calculates the equilibrium values for each eigenvalue of the Hamiltonian,
      normalizes them, and transforms the resulting diagonal matrix to the local basis.
    """

    tb_ham = me_solver.tb_ham
    tb_ham.relaxation = False
    eigv, eigs = tb_ham.get_eigensystem()
    temperature = me_solver.lindblad_diss.temperature

    # Ground state
    if temperature == 0:
        return tb_ham.eigs[:, 0]

    # Thermal equilibrium state
    # Calculate the equilibrium values for each eigenvalue
    eq_values = np.zeros(tb_ham.matrix_dim)
    for i in range(tb_ham.matrix_dim):
        energy = eigv[i] * get_conversion(tb_ham.unit, "J")
        eq_values[i] = np.exp(-energy / (c.k * temperature))

    # Normalize the equilibrium values, convert to a diagonal matrix and transform to the local basis
    eq_values /= np.sum(eq_values)
    eq_values = np.diag(eq_values)
    eq_state = global_to_local(eq_values, eigs)

    return eq_state


def get_deph_eq_state(me_solver):
    """Calculate the dephasing equilibrium state.

    Parameters
    ----------
    me_solver : object
        An instance of :class:`ME_Solver` solver which contains the following attributes:

    Returns
    -------
    numpy.ndarray
        The dephasing equilibrium state as a density matrix.
    """
    deph_eq_state = None

    # Local dephasing
    if me_solver.lindblad_diss.loc_deph_rate:
        # maximally mixed state
        dim = me_solver.tb_ham.matrix_dim
        if me_solver.tb_ham.relaxation:
            dim -= 1
        deph_eq_state = np.eye(dim) / dim

    # Global dephasing
    if me_solver.lindblad_diss.glob_deph_rate:
        loc_init_matrix = me_solver.init_matrix.full().real
        if me_solver.tb_ham.relaxation:
            loc_init_matrix = delete_groundstate(loc_init_matrix)
        _, eigs = me_solver.tb_ham.get_eigensystem()
        glob_init_matrix = local_to_global(loc_init_matrix, eigs)

        # cancel all off-diagonal elements
        glob_init_matrix = np.diag(np.diag(glob_init_matrix))
        loc_init_matrix = global_to_local(glob_init_matrix, eigs)
        deph_eq_state = loc_init_matrix

    assert deph_eq_state is not None, "Dephasing equilibrium state not calculated."
    return deph_eq_state
