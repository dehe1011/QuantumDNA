import numpy as np
import qutip as q

from ..model import global_to_local
from ..hamiltonian import add_groundstate
from .observables import get_pop_observable

# ----------------------------------------------------------------------


def get_loc_deph_ops(tb_basis, description, dephasing_dict, relaxation):
    """Generate a list of local dephasing collapse operators.

    Parameters
    ----------
    tb_basis : list
        Tight-binding basis.
    description : object
        A descriptor object that provides information about the system's structure.
    dephasing_dict : dict
        A dictionary where keys are particle identifiers and values are their corresponding dephasing rates.
    relaxation : bool
        If True, the ground state is added to the operators.

    Returns
    -------
    list
       List of local dephasing operators.
    """

    c_ops = []
    for tb_site in tb_basis:
        for particle, dephasing_rate in dephasing_dict.items():
            op = get_pop_observable(tb_basis, description, particle, tb_site)
            if relaxation:
                op = add_groundstate(op)
            c_ops.append(np.sqrt(dephasing_rate) * q.Qobj(op))
    return c_ops


# ----------------------------------------------------------------------


def get_glob_deph_ops(eigs, dephasing_rate, relaxation):
    """Generate a list of global dephasing collapse operators. In total :math:`N^2` operators (where :math:`N` is
    the number of eigenstates).

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
        local_op = q.fock_dm(num_eigenstates, i).full()
        global_op = global_to_local(local_op, eigs)
        if relaxation:
            global_op = add_groundstate(global_op)
        c_ops.append(np.sqrt(dephasing_rate) * q.Qobj(global_op))
    return c_ops


# ----------------------------------------------------------------------
