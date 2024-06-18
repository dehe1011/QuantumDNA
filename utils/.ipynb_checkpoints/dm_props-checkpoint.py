import numpy as np

__all__ = ['calc_trace_distance', 'calc_purity', 'calc_coherence', 'calc_ipr_dm', 'calc_ipr_hamiltonian']

# ----------------------- properties of density matrices (quantum states) ----------------------------

def calc_trace_distance(density_matrix_1, density_matrix_2): 
    # trace distance of two density matrices 
    eigenvals = np.linalg.eig(density_matrix_1 - density_matrix_2)[0]
    diag = np.diag(np.absolute(eigenvals))
    return 0.5 * np.trace(diag)

def calc_purity(density_matrix):
    # calculates the purity of a density matrix
    return np.real( np.trace( np.matmul(density_matrix,density_matrix) ) )

def calc_coherence(density_matrix):
    # calculates the coherence (absolute sum of the off-diagonals) of a density matrix 
    coherence = 0
    for row in density_matrix:
        for element in row:
            coherence += abs(element)
    return coherence - 1 

def calc_ipr_dm(density_matrix):
    # calculates the inverse participation ratio (IPR) of a density matrix
    # 1/N (localized, coherences in some basis) < IPR < N (delocalized, coherent), IPR = 1 for maximally mixed state (incoherent in any basis)
    dims = density_matrix.shape[0]
    numer, denom = 0, 0
    for row in density_matrix:
        for element in row:
            numer += abs(element)
            denom += abs(element)**2
    return numer**2 / (dims * denom)

def calc_ipr_hamiltonian(hamiltonian):
    # 1 (localized eigenstate/ exciton state) < IPR < N (delocalized eigenstate/ exciton state)
    # calculates the inverse participation ratio (IPR) for each eigenstate of the Hamiltonian
    dims = hamiltonian.shape[0]
    states = np.linalg.eig(hamiltonian)[1]
    ratios = np.zeros(dims)
    for idx, state in enumerate(states):
        tmp = 0
        for coeff in state:
            tmp += coeff ** 4
        ratios[idx] = tmp ** -1
    return ratios
    