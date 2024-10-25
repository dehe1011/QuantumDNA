import numpy as np
import qutip as q

from ..model import global_to_local
from ..hamiltonian import add_groundstate
from .observables import get_pop_particle

# ----------------------------------------------


def get_loc_deph_ops(tb_basis, dephasing_rate, relaxation):
    """
    Local dephasing operators. In total :math:`2N` operators (where :math:`N` is the number of tight-binding sites).

    Parameters
    ----------
    tb_basis : list
        List of tight-binding site basis states.
    dephasing_rate : float
        Dephasing rate.
    relaxation : bool
        Flag for relaxation.

    Returns
    -------
    list
        List of local dephasing operators.
    """

    c_ops = []
    for tb_site in tb_basis:
        op_electron = get_pop_particle(tb_basis, "electron", tb_site)
        op_hole = get_pop_particle(tb_basis, "hole", tb_site)
        if relaxation:
            op_electron = add_groundstate(op_electron)
            op_hole = add_groundstate(op_hole)
        c_ops.append(np.sqrt(dephasing_rate) * q.Qobj(op_electron))
        c_ops.append(np.sqrt(dephasing_rate) * q.Qobj(op_hole))
    return c_ops


def get_glob_deph_ops(eigs, dephasing_rate, relaxation):
    """
    Global dephasing operators. In total :math:`N^2` operators (where :math:`N` is the number of eigenstates).

    Parameters
    ----------
    eigs : np.ndarray
        Eigensystem.
    dephasing_rate : float
        Dephasing rate.
    relaxation : bool
        Flag for relaxation.

    Returns
    -------
    list
        List of global dephasing operators.
    """
    num_eigenstates = eigs.shape[0]
    c_ops = []
    for i in range(num_eigenstates):
        matrix = np.zeros((num_eigenstates, num_eigenstates))
        matrix[i, i] = 1
        op = global_to_local(matrix, eigs)
        if relaxation:
            op = add_groundstate(op)
        c_ops.append(np.sqrt(dephasing_rate) * q.Qobj(op))
    return c_ops


def get_loc_deph_p_ops(tb_basis, dephasing_rate):
    """
    Local dephasing operators for particle description. In total :math:`N` operators (where :math:`N` is the number of tight-binding sites).

    Parameters
    ----------
    tb_basis : list
        List of tight-binding site basis states.
    dephasing_rate : float
        Dephasing rate.

    Returns
    -------
    list
        List of local dephasing operators for particle description.
    """

    num_sites = len(tb_basis)
    c_ops = []
    for i in range(num_sites):
        op = np.sqrt(dephasing_rate) * q.fock_dm(num_sites, i)
        c_ops.append(op)
    return c_ops


def get_glob_deph_p_ops(eigs, dephasing_rate):
    """
    Global dephasing operators for particle description. In total :math:`N` operators (where :math:`N` is the number of eigenstates).

    Parameters
    ----------
    eigs : np.ndarray
        Eigensystem.
    dephasing_rate : float
        Dephasing rate.

    Returns
    -------
    list
        List of global dephasing operators for particle description.
    """

    num_eigenstates = eigs.shape[0]
    c_ops = []
    for i in range(num_eigenstates):
        local_op = q.fock_dm(num_eigenstates, i).full()
        global_op = global_to_local(local_op, eigs)
        c_ops.append(np.sqrt(dephasing_rate) * q.Qobj(global_op))
    return c_ops
