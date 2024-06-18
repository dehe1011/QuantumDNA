from scipy import linalg, constants
import numpy as np
from qutip.nonmarkov.dlheom_solver import HSolverDL
from qutip import Qobj
from utils import basis_change, calc_purity, trace_distance

TEMP_INDEP_MODELS = ['local dephasing lindblad']
TEMP_DEP_MODELS = ['global thermalising lindblad', 'local thermalising lindblad', 'HEOM']
DYNAMICS_MODELS = TEMP_INDEP_MODELS + TEMP_DEP_MODELS

def get_initial_dm(dims, init_site_pop):
    # init_site_pop is a list, e.g., [1,6] 
    num_states = len(init_site_pop)
    state = np.zeros(dims)
    for site in init_site_pop:
        state[site - 1] = 1
    state /= np.sqrt(num_states) # normalization
    return np.outer(state, state)

def get_eq_state(dynamics_model, dims, hamiltonian, temperature):
    # returns the equilibrium state of the system
    # maximally-mixed state (infinite temperature limit)
    if dynamics_model in TEMP_INDEP_MODELS:
        return np.eye(dims, dtype=complex) / dims
        
    # thermal equilibrium: frac{e^{- H / k_B T}}{tr(e^{- H / k_B T})}
    eigv, eigs = np.linalg.eig(hamiltonian)
    eq_values = np.zeros_like(eigv)
    for idx, _ in enumerate(eq_values):
        eq_values[idx] = np.exp(- eigv[idx] * 1e12 * constants.hbar / (constants.k * temperature))
    eq_values = np.diag( eq_values / np.sum(eq_values) )
    eq_state = basis_change(eq_values, eigs, liouville = False)
    return eq_state

def evolve_matrix_one_step(dm, superop, time_interval):
    dims = dm.shape[0]
    propagator = linalg.expm(superop * time_interval)
    evolved = np.matmul(propagator, dm.flatten('C'))
    evolved = evolved.reshape((dims, dims), order='C')
    return evolved

# time in fs
def time_evo_lindblad(dm, superop, timesteps, time_interval, dynamics_model, hamiltonian, temperature):

    dims = dm.shape[0]
    # convert from fs to ps (to match superoperator units)
    time_interval = time_interval * 1e-3
    
    evolution = []
    time, dm = 0., dm
    squared = calc_purity(dm)
    eq_state = get_eq_state(dynamics_model, dims, hamiltonian, temperature)
    distance = trace_distance(dm, eq_state)
    evolution.append( [time, dm, squared, distance] )
    
    for step in range(1, timesteps + 1):
        time += time_interval
        dm = evolve_matrix_one_step(dm, superop, time_interval)
        dm /= np.trace(dm) # renormalize dm
        # convert time back from ps to fs
        evolution.append( [time * 1e3, dm, calc_purity(dm), trace_distance(dm, eq_state)] ) 
    return evolution

# time in ps
def time_evo_heom(dm, timesteps, time_interval, hamiltonian, coupling_op, reorg_energy, temperature, bath_cutoff,
                  matsubara_terms, cutoff_freq, matsubara_coeffs, matsubara_freqs):

    dims = dm.shape[0]
    hsolver = HSolverDL( Qobj(hamiltonian), Qobj(coupling_op), reorg_energy, temperature, bath_cutoff, matsubara_terms,
                        cutoff_freq, planck=1.0, boltzmann=1.0, renorm=True, stats=False )

    hsolver.exp_coeff = matsubara_coeffs
    hsolver.exp_freq = matsubara_freqs
    times = np.array(range(timesteps + 1)) * time_interval  # ps
    result = hsolver.run(Qobj(dm), times)
    
    temperature = temperature *(constants.hbar * 1e12)/constants.k # convert rad/ps to K
    eq_state = get_eq_state('HEOM', dims, hamiltonian, temperature)
    evolution = []
    for i in range(0, len(result.states)):
        dm = np.array(result.states[i]).T
        dm /= np.trace(dm) # renormalize
        # convert from ps to fs
        evolution.append( [float(result.times[i]) * 1e3, dm, calc_purity(dm), trace_distance(dm, eq_state)] )
    return evolution, np.array(hsolver.exp_coeff), np.array(hsolver.exp_freq)

def process_evo_data(time_evolution, elements, trace_measure):
    # returns elements like ['11', '21', ...]
    times = np.empty(len(time_evolution), dtype=float)
    matrix_data = ({element: np.empty(len(time_evolution), dtype=complex)
                    for element in elements} if elements else None)
    squared = (np.empty(len(time_evolution), dtype=float)
               if 'squared' in trace_measure else None)
    distance = (np.empty(len(time_evolution), dtype=float)
                if 'distance' in trace_measure else None)
    for idx, (time, rho_t, squ, dist) in enumerate(time_evolution, start=0):
        times[idx] = time  # already in fs
        if matrix_data is not None:
            for element in elements:
                n, m = int(element[0]) - 1, int(element[1]) - 1
                matrix_data[element][idx] = rho_t[n][m]
        # Process trace measure data
        if squared is not None:
            squared[idx] = squ
        if distance is not None:
            distance[idx] = dist
    return times, matrix_data, squared, distance