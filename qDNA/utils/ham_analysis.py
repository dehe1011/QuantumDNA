from itertools import chain, product
import numpy as np

# Shortcuts:
# IPR: inverse participation ratio
# calc: calculate
# pop: population
# eigs: eigenstates
# val: value
# dims: dimensions
# init: initial

__all__ = ['calc_average_pop', 'calc_amplitudes', 'calc_frequencies', 'get_pop_fourier', 'calc_ipr_hamiltonian']

def calc_average_pop(eigs: np.ndarray, init_state: int, end_state: int) -> float:
    """
    Calculates the time-averaged/ static population of the end state when initially in the initial state.

    Args:
        eigs: Eigenvector matrix.
        init_state: Initial state in the basis.
        end_state: End state in the basis.
        basis: The list of basis elements.

    Returns:
        The static population of the end state.
    """
    average_pop = np.sum(eigs[end_state, :] ** 2 * eigs[init_state, :] ** 2)
    return average_pop

def calc_amplitudes(eigs: np.ndarray, init_state: int, end_state: int) -> np.ndarray:
    """
    Calculates the amplitudes for transitions between states.

    Args:
        eigs: Eigenvector matrix.
        init_state: Initial state in the basis.
        end_state: End state in the basis.
        basis: The list of basis elements.

    Returns:
        A list of amplitudes for transitions.
    """
    matrix_dimension = eigs.shape[0]
    amplitudes = [2 * eigs[end_state, i] * eigs[init_state, i] * eigs[end_state, j] * eigs[init_state, j]
                  for i, j in product(range(matrix_dimension), repeat=2) if i<j]
    return np.real(np.array(amplitudes))

def calc_frequencies(eigv: np.ndarray) -> np.ndarray:
    """
    Calculates the frequencies corresponding to energy level differences.

    Args:
        eigv: Eigenvalue vector.

    Returns:
        A list of frequencies.
    """
    matrix_dimension = len(eigv)
    frequencies = [abs(eigv[i] - eigv[j]).real for i, j in product(range(matrix_dimension), repeat=2) if i<j]
    return np.array(frequencies)

def get_pop_fourier(t: float, average_pop: float, amplitudes: np.ndarray, frequencies: np.ndarray) -> float:
    return average_pop + np.sum( [amplitude * np.cos(frequency * t) for amplitude, frequency in zip(amplitudes, frequencies)] )

def calc_ipr_hamiltonian(eigs):
    """
    1 (localized eigenstate/ exciton state) < IPR < N (delocalized eigenstate/ exciton state)
    calculates the inverse participation ratio (IPR) for each eigenstate of the Hamiltonian
    """
    dims = eigs.shape[0]
    ratios = []
    for idx in range(dims):
        ipr_val = 1/np.sum([coeff**4 for coeff in eigs[:,idx] ])
        ratios.append( ipr_val )
    return ratios
    
