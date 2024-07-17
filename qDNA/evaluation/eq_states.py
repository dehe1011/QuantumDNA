"""
This module provides functions to compute equilibrium states for quantum DNA models.
"""

import numpy as np
import scipy.constants as c

from qDNA.utils import get_conversion
from qDNA.model import global_to_local, local_to_global

# Shortcuts:
# eq: equilibrium

__all__ = ["get_therm_eq_state", "get_deph_eq_state"]

# ----------------------------------- Equilibrium States --------------------------------------


def get_therm_eq_state(me_solver):
    """
    Computes the thermal equilibrium state.

    Parameters
    ----------
    me_solver : MESolverType
        The master equation solver instance.

    Returns
    -------
    np.ndarray
        The thermal equilibrium state.
    """
    tb_ham = me_solver.tb_ham
    eigv, eigs = tb_ham.get_eigensystem()
    temperature = me_solver.lindblad_diss.temperature
    if temperature == 0:
        return tb_ham.eigs[:, 0]  # ground state

    eq_values = np.zeros(tb_ham.matrix_dim)
    for i in range(tb_ham.matrix_dim):
        eq_values[i] = np.exp(
            -eigv[i] * get_conversion(tb_ham.unit, "J") / (c.k * temperature)
        )
    eq_values = np.diag(eq_values / np.sum(eq_values))
    eq_state = global_to_local(eq_values, eigs)
    return eq_state


def get_deph_eq_state(me_solver):
    """
    Computes the dephasing equilibrium state.

    Parameters
    ----------
    me_solver : MESolverType
        The master equation solver instance.

    Returns
    -------
    np.ndarray
        The dephasing equilibrium state.
    """
    if me_solver.lindblad_diss.loc_deph_rate:
        return (
            np.eye(me_solver.tb_ham.matrix_dim) / me_solver.tb_ham.matrix_dim
        )  # maximally mixed state
    if me_solver.lindblad_diss.glob_deph_rate:
        loc_init_matrix = me_solver.init_matrix.full()
        _, eigs = me_solver.tb_ham.get_eigensystem()
        glob_init_matrix = local_to_global(loc_init_matrix, eigs)
        glob_init_matrix = np.diag(
            np.diag(glob_init_matrix)
        )  # cancel all off-diagonal elements
        loc_init_matrix = global_to_local(glob_init_matrix, eigs)
        return loc_init_matrix
