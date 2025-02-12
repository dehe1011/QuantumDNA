import numpy as np
import qutip as q

# ----------------------------------------------------


def get_relax_op(tb_basis, tb_site):
    """Annihilation operator of an exciton on a given tight-binding site. Relaxation of
    the DNA to its ground state.

    Parameters
    ----------
    tb_basis : list
        List of tight-binding site basis states.
    tb_site : str
        Tight-binding site.

    Returns
    -------
    qutip.Qobj
        Relaxation operator.
    """

    tb_site_idx = tb_basis.index(tb_site)
    num_sites = len(tb_basis)
    relax_op = np.zeros((num_sites**2 + 1, num_sites**2 + 1))
    relax_op[0, 1 + tb_site_idx * (num_sites + 1)] = 1
    return q.Qobj(relax_op)


def get_relax_ops(tb_basis, tb_basis_sites_dict, relax_rates):
    """
    Generate relaxation operators for a given tight-binding basis.
    Parameters
    ----------
    tb_basis : list
        List of tight-binding basis states.
    tb_basis_sites_dict : dict
        Dictionary mapping tight-binding basis states to site indices.
    relax_rates : dict
        Dictionary mapping site indices to relaxation rates.
    Returns
    -------
    relax_ops : list
        List of relaxation operators, each scaled by the square root of the corresponding relaxation rate.
    """

    relax_ops = []
    for tb_site in tb_basis:
        relax_rate = relax_rates[tb_basis_sites_dict[tb_site]]
        if relax_rate != 0:
            relax_op = get_relax_op(tb_basis, tb_site)
            relax_ops.append(np.sqrt(relax_rate) * relax_op)
    return relax_ops
