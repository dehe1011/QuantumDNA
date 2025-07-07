from itertools import permutations, product
import numpy as np
import qutip as q

from ..model import global_to_local
from ..hamiltonian import add_groundstate
from .therm_rates import rate_constant_redfield

# ----------------------------------------------------------------------


def get_glob_therm_op(eigs, eigenstate_i, eigenstate_j, relaxation, matrix_dim):
    """Global thermalizing operator.

    Parameters
    ----------
    eigs : np.ndarray
        Eigensystem.
    eigenstate_i : int
        Index of the initial eigenstate.
    eigenstate_j : int
        Index of the final eigenstate.
    relaxation : bool
        Flag for relaxation.
    matrix_dim : int
        Dimension of the matrix.

    Returns
    -------
    qutip.Qobj
        Thermalizing operator.
    """

    op = np.zeros((matrix_dim, matrix_dim), dtype=complex)
    op[eigenstate_i, eigenstate_j] = 1
    op = global_to_local(op, eigs)
    if relaxation:
        op = add_groundstate(op)
    return q.Qobj(op)


def get_glob_therm_ops(eigv, eigs, relaxation, **kwargs):
    """Generate global thermalizing operators.

    Parameters
    ----------
    eigv : array_like
        Eigenvalues of the system Hamiltonian.
    eigs : array_like
        Eigenvectors of the system Hamiltonian.
    relaxation : float
        Relaxation rate for the system.
    deph_rate : float, optional
        Dephasing rate (default is 7).
    cutoff_freq : float, optional
        Cutoff frequency for the spectral density (default is 20).
    reorg_energy : float, optional
        Reorganization energy (default is 1).
    temperature : float, optional
        Temperature of the thermal bath (default is 300).
    spectral_density : str, optional
        Type of spectral density function (default is "debye").
    exponent : float, optional
        Exponent for the spectral density function (default is 1).

    Returns
    -------
    list
        List of collapse operators for the Lindblad master equation.
    """

    matrix_dim = eigs.shape[0]
    c_ops = []
    for eigs_i, eigs_j in permutations(range(matrix_dim), 2):
        # Calculate Lindblad rate
        omega_i, omega_j = eigv[eigs_i], eigv[eigs_j]
        omega = omega_i - omega_j
        lind_rate = rate_constant_redfield(omega, **kwargs)
        # Calculate thermalizing operator
        lind_op = get_glob_therm_op(eigs, eigs_i, eigs_j, relaxation, matrix_dim)
        # Append to the list
        c_ops.append(np.sqrt(lind_rate) * lind_op)

    return c_ops


# ----------------------------------------------------------------------


def get_loc_therm_op(eigv, eigs, unique, site_m, relaxation, matrix_dim):
    """Local thermalizing operator.

    Parameters
    ----------
    eigv : np.ndarray
        Eigenvalues.
    eigs : np.ndarray
        Eigensystem.
    unique : float
        Unique frequency gap.
    site_m : int
        Local site index.
    relaxation : bool
        Flag for relaxation.
    matrix_dim : int
        Dimension of the matrix.

    Returns
    -------
    qutip.Qobj
        Thermalizing operator.
    """

    op = np.zeros((matrix_dim, matrix_dim), dtype=complex)
    for i, j in product(range(matrix_dim), repeat=2):
        omega_i, omega_j = eigv[i], eigv[j]
        state_i, state_j = eigs[:, i], eigs[:, j]
        if omega_i - omega_j == unique:
            op += (
                state_j[site_m].conjugate()
                * state_i[site_m]
                * np.outer(state_j, state_i)
            )
    if relaxation:
        op = add_groundstate(op)

    return q.Qobj(op)


def get_loc_therm_ops(eigv, eigs, relaxation, **kwargs):
    """Generate local thermalizing operators.

    Parameters
    ----------
    eigv : np.ndarray
        Eigenvalues.
    eigs : np.ndarray
        Eigensystem.
    relaxation : bool
        Flag for relaxation.
    deph_rate : float
        Dephasing rate.
    cutoff_freq : float
        Cutoff frequency.
    reorg_energy : float
        Reorganization energy.
    temperature : float
        Temperature in Kelvin.
    spectral_density : str
        Type of spectral density.
    exponent : float
        Exponent for Ohmic spectral density.

    Returns
    -------
    list
        List of thermalizing operators.
    """

    matrix_dim = len(eigv)
    c_ops = []

    # Calculate unique frequency gaps
    gaps = eigv.reshape(matrix_dim, 1) - eigv
    unique = np.unique(gaps.flatten())

    for unique, site_m in product(unique, range(matrix_dim)):
        # Calculate Lindblad rate
        omega = unique
        lind_rate = rate_constant_redfield(omega, **kwargs)
        # Calculate local thermalizing operator
        lind_op = get_loc_therm_op(eigv, eigs, unique, site_m, relaxation, matrix_dim)
        # Append to the list
        c_ops.append(np.sqrt(lind_rate) * lind_op)

    return c_ops


# ----------------------------------------------------------------------
