from itertools import product, chain
import numpy as np
from typing import Any, List, Tuple, Union
from .tb_basis import basis_converter, get_particle_eh_states
from .tb_ham import TBHamType, delete_basis_dimension

def check_input(state: Any, basis: List[Any]) -> Tuple[int, int]:
    """
    Validates that the state is present in the provided basis and converts them to the appropriate representation.

    """
    if state not in basis:
        raise ValueError(f"State {state} not in basis {basis}")
    state = basis_converter(state, basis)
    return state

def calc_average_pop(eigs: np.ndarray, init_state: int, end_state: int) -> float:
    """
    Calculates the static population of the end state when initially in the initial state.

    Args:
        eigs: Eigenvector matrix.
        init_state: Initial state in the basis.
        end_state: End state in the basis.
        basis: The list of basis elements.

    Returns:
        The static population of the end state.
    """
    average_pop = np.sum(eigs[end_state, :] ** 2 * eigs[init_state, :] ** 2)
    return float(average_pop)

def calc_amplitudes(eigs: np.ndarray, init_state: int, end_state: int) -> List[float]:
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
                  for i, j in product(range(matrix_dimension), repeat=2) if i < j]
    return [float(amp.real) for amp in amplitudes]

def calc_frequencies(eigv: np.ndarray) -> List[float]:
    """
    Calculates the frequencies corresponding to energy level differences.

    Args:
        eigv: Eigenvalue vector.

    Returns:
        A list of frequencies.
    """
    matrix_dimension = len(eigv)
    frequencies = [abs(eigv[i] - eigv[j]).real for i, j in product(range(matrix_dimension), repeat=2) if i < j]
    return [float(freq) for freq in frequencies]

# -----------------------------------------------------------------------------------------------------------

def fourier_analysis_1P(tb_ham: TBHamType, init_state: Any, end_state: Any) -> Tuple[float, List[float], List[float]]:
    """
    Performs Fourier analysis for a single particle system.

    Args:
        hamiltonian: The Hamiltonian of the system.
        init_state: The initial state in the tight-binding basis.
        end_state: The end state in the tight-binding basis.

    Returns:
        A tuple containing the average population, amplitudes, and frequencies.

    Examples:
        >>> tb_model = TB_Model(4, 'WM')
        >>> hamiltonian = Hamiltonian('GGGG', tb_model, Ham_kwargs={'relaxation': True, 'particle': 'electron'})
        >>> fourier_analysis_one_particle(hamiltonian, '(0, 0)', '(0, 3)')
    """
    tb_basis = tb_ham.tb_basis
    init_state = check_input(init_state, tb_basis)
    end_state = check_input(end_state, tb_basis)
    average_pop = calc_average_pop(tb_ham.eigs, init_state, end_state)
    amplitudes = calc_amplitudes(tb_ham.eigs, init_state, end_state)
    frequencies = calc_frequencies(tb_ham.eigv)
    return average_pop, amplitudes, frequencies

def fourier_analysis_2P(tb_ham: TBHamType, init_state: Any, end_state: Any, particle: str) -> Tuple[float, List[float], List[float]]:
    """
    Performs Fourier analysis for a two-particle system.

    Args:
        hamiltonian: The Hamiltonian of the system.
        init_state: The initial state in the electron-hole basis.
        end_state: The end state in the tight-binding basis.
        particle: The type of particle ('electron', 'hole', or 'exciton').

    Returns:
        A tuple containing the average population, amplitudes, and frequencies.

    Examples:
        >>> tb_model = TB_Model(4, 'WM')
        >>> hamiltonian = Hamiltonian('GGGG', tb_model, Ham_kwargs={'relaxation': True, 'particle': 'exciton'})
        >>> fourier_analysis_two_particles(hamiltonian, ('(0, 0)', '(0, 0)'), '(0, 3)', 'electron')
    """
    assert tb_ham.particle == "exciton", "You must choose the electron-hole Hamiltonian. Please use fourier_analysis_1P() if you only want to describe the dynamics of one particle. "
    eh_basis = tb_ham.eh_basis
    tb_basis = tb_ham.tb_basis
    init_state = check_input(init_state, eh_basis)
    end_states = [check_input(end_state, eh_basis) for end_state in get_particle_eh_states(particle, end_state, tb_basis)]
    average_pop = np.sum([calc_average_pop(tb_ham.eigs, init_state, es) for es in end_states])
    amplitudes = list(chain.from_iterable([calc_amplitudes(tb_ham.eigs, init_state, es) for es in end_states]))
    frequencies = calc_frequencies(tb_ham.eigv) * len(end_states)
    return average_pop, amplitudes, frequencies

def get_pop_fourier(t: float, average_pop: float, amplitudes: List[float], frequencies: List[float]) -> float:
    return average_pop + np.sum( [amplitude * np.cos(frequency * t) for amplitude, frequency in zip(amplitudes, frequencies)] )

    
