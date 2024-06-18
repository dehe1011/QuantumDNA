import numpy as np

# ----------------------- properties of density matrices (quantum states) ----------------------------

def trace_distance(density_matrix_1, density_matrix_2): 
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
    dims = density_matrix.shape[0]
    numer, denom = 0, 0
    for row in density_matrix:
        for element in row:
            numer += abs(element)
            denom += abs(element)**2
    return numer**2 / (dims * denom)

def calc_ipr_hamiltonian(hamiltonian):
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

def basis_change(matrix, states, liouville = False):
    # performs a basis change of the given matrix  
    # states contain the new basis expressed as vector in the old basis
    if liouville: states = np.kron(states, states.conjugate())
    # for open quantum systems the dimenison of the matrices is N**2 instead of N
    return np.matmul(states, np.matmul(matrix, states.conjugate().T))


