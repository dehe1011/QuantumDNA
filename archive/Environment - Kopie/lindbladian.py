from itertools import permutations, product
import numpy as np
from bath import rate_constant_redfield
from utils import basis_change

LINDBLAD_MODELS = ['local dephasing lindblad', 'global thermalising lindblad', 'local thermalising lindblad']

def loc_deph_lindblad_op(dims, site_j):
    l_op = np.zeros((dims, dims), dtype=complex)
    l_op[site_j][site_j] = 1
    return l_op

def glob_therm_lindblad_op(dims, state_a, state_b):
    # returns Lindblad operator in the eigenbasis -> perform basis change
    # eigenstates a and b have to be different!
    l_op = np.zeros((dims, dims), dtype=complex)
    l_op[state_a][state_b] = 1
    return l_op

def loc_therm_lindblad_op(eigv, eigs, unique, site_m):
    dims = len(eigv)
    l_op = np.zeros((dims, dims), dtype=complex)
    for idx_i, idx_j in product(range(dims), repeat=2):
        omega_i, omega_j = eigv[idx_i], eigv[idx_j]
        state_i, state_j = eigs[:, idx_i], eigs[:, idx_j]
        if omega_i - omega_j == unique:
            l_op += (state_j[site_m].conjugate() * state_i[site_m] * np.outer(state_j, state_i))
    return l_op

def lindblad_superop_element(l_op):
    # constructs a single element of the Libndblad dissipator 
    l_op_dag = l_op.T.conjugate()
    id = np.eye(l_op.shape[0])
    return np.kron(l_op.conjugate(), l_op) - 0.5 * (np.kron(np.matmul(l_op_dag, l_op).conjugate(), id) + np.kron(id, np.matmul(l_op_dag, l_op)))

def lindbladian_superop(dims, dynamics_model, hamiltonian = None, deph_rate = None, cutoff_freq = None, 
                        reorg_energy = None, temperature = None, spectral_density = None, exponent = 1):
    # constructs the whole Lindblad dissipator 
    lindbladian = np.zeros((dims ** 2, dims ** 2), dtype=complex)

    if dynamics_model == 'local dephasing lindblad':
        for site_j in range(dims):
            l_op = loc_deph_lindblad_op(dims, site_j)
            superop_element = lindblad_superop_element(l_op)
            lindbladian += superop_element
        return deph_rate * lindbladian  # unit: rad ps^-1

    eigv, eigs = np.linalg.eig(hamiltonian)

    if dynamics_model == 'global thermalising lindblad':
        for state_a, state_b in permutations(range(dims), 2):
            omega_a, omega_b = eigv[state_a], eigv[state_b]
            k_ab = rate_constant_redfield( (omega_a - omega_b), deph_rate, cutoff_freq, reorg_energy, temperature, spectral_density, exponent )
            l_op = glob_therm_lindblad_op(dims, state_a, state_b)
            superop_element = lindblad_superop_element(l_op)
            # we perform a basis change from the eigenstate (global) basis to the site (local) basis  
            superop_element = basis_change(superop_element, eigs, liouville=True)
            lindbladian += k_ab * superop_element
        return lindbladian  # unit: rad ps^-1

    if dynamics_model == 'local thermalising lindblad':
        gaps = eigv.reshape(dims, 1) - eigv # matrix that contains all possible eigenenergy differences
        unique = np.unique(gaps.flatten()) 
        for unique, site_m in product(unique, range(dims)):
            k_omega = rate_constant_redfield( unique, deph_rate, cutoff_freq, reorg_energy, temperature, spectral_density, exponent )
            l_op = loc_therm_lindblad_op(eigv, eigs, unique, site_m)
            superop_element = lindblad_superop_element(l_op)
            lindbladian += k_omega * superop_element
        return lindbladian  # unit: rad ps^-1

def hamiltonian_superop(hamiltonian):
    # constructs the unitary evolution of the system
    dims = hamiltonian.shape[0]
    id = np.identity(dims)
    return -1.0j * ( np.kron(hamiltonian, id) - np.kron(id, hamiltonian.T.conjugate()) ) 

