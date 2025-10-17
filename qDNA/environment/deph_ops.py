import numpy as np
import qutip as q

from ..model import global_to_local
from ..hamiltonian import add_groundstate
from .observables import get_pop_observable

# ----------------------------------------------------------------------


def get_loc_deph_ops(tb_basis, description, dephasing_dict, relaxation):
    """Generate a list of local dephasing collapse operators.

    This function creates local dephasing operators based on the provided tight-binding basis,
    system description, dephasing rates, and relaxation flag.

    Parameters
    ----------
    tb_basis : list
        List representing the tight-binding basis states.
    description : object
        Descriptor object containing structural information about the system.
    dephasing_dict : dict
        Dictionary mapping particle identifiers to their respective dephasing rates.
    relaxation : bool
        If True, the ground state is included in the operators.

    Returns
    -------
    list
        A list of local dephasing operators as Qobj instances.
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
    """Generate a list of global dephasing collapse operators.

    This function creates :math:`N^2` operators, where :math:`N` is the number of eigenstates.

    Parameters
    ----------
    eigs : np.ndarray
        Array representing the eigensystem of the Hamiltonian.
    dephasing_rate : float
        The rate of dephasing applied to the system.
    relaxation : bool
        If True, the ground state is included in the operators.

    Returns
    -------
    list
        A list containing the global dephasing operators as Qobj instances.
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
