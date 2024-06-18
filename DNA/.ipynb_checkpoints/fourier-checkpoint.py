from utils import pop_electrons, pop_holes, pop_excitons, real_units
from .hamiltonian import Hamiltonian 
import numpy as np

__all__ = ['static_pop_tot', 'calc_static_pop_dir', 'calc_fourier', 'calc_fourier_dir', 'get_pop_lists']

def get_basis_states_list(N, base, particle):
    # returns all of the basis states that describe the presence of a particle at a given base 
    i = base
    if particle == 'electron': basis = pop_electrons(N,i,False).diag()
    elif particle == 'hole': basis = pop_holes(N,i,False).diag()
    elif particle == 'exciton': basis = pop_excitons(N,i,False).diag()
    basis_states_list=[basis_state for basis_state in range(N**2) if basis[basis_state]==1]
    return basis_states_list

# ----------------------------------- Average populations -------------------------------------------

def static_pop_tot(DNAstring, init_base = 0, Ham_kwargs={}):
    # returns the time-averaged population for all particles and bases 
    Ham_kwargs['relaxation'] = False
    H = Hamiltonian(DNAstring, Ham_kwargs=Ham_kwargs)
    w,v = np.linalg.eigh(H.matrix)
    pop_static_e=[calc_static_pop(w,v, 'electron', base, init_base = init_base) for base in range(H.N)]
    pop_static_h=[calc_static_pop(w,v, 'hole', base, init_base = init_base) for base in range(H.N)]
    pop_static_ex=[calc_static_pop(w,v, 'exciton', base, init_base = init_base) for base in range(H.N)]
    return pop_static_e, pop_static_h, pop_static_ex

def calc_static_pop(w,v, particle, base, init_base = 0):
    # returns the time-averaged population for a selected base and particle
    N = int( np.sqrt(len(w)) )
    init_state = init_base * (N+1)
    basis_states_list = get_basis_states_list(N, base, particle)
    pop_static = 0
    for basis_state in basis_states_list:
        for i in range(N**2):
            pop_static += v[basis_state,i]**2*v[init_state,i]**2
    return pop_static.real

def calc_static_pop_dir(DNAstring, particle, base, init_base=0, Ham_kwargs={}):
    # wrapper of calc_static_pop()
    Ham_kwargs['relaxation'] = False
    H = Hamiltonian(DNAstring, Ham_kwargs=Ham_kwargs)
    w,v = np.linalg.eigh(H.matrix)
    pop_static = calc_static_pop(w,v, particle, base, init_base=init_base)
    return pop_static

# ------------------------------------- Fourier spectrum ---------------------------------------------

def calc_fourier(w,v, particle, base, init_base=0):
    # be careful: sometimes the absolute value of the amplitude is feasible
    # returns the amplitude and frequencies (Fourier spectrum) for a selected base and particle 
    N = int( np.sqrt(len(w)) )
    init_state = init_base * (N+1)
    basis_states_list = get_basis_states_list(N, base, particle)
    amplitude,frequency=[],[]
    for basis_state in basis_states_list:
        amplitude_state=[(2*v[basis_state,i]*v[init_state,i]*v[basis_state,j]*v[init_state,j]).real for i in range(N**2) for j in range(N**2) if i<j]
        frequency_state=[ (abs(w[i]-w[j])/(2*np.pi*real_units/1000)).real for i in range(N**2) for j in range(N**2) if i<j]
        amplitude.extend(amplitude_state)
        frequency.extend(frequency_state)
    return amplitude, frequency

def calc_fourier_dir(DNAstring, particle, base, init_base=0, Ham_kwargs={}):
    # wrapper of calc_fourier()
    Ham_kwargs['relaxation'] = False
    H = Hamiltonian(DNAstring, Ham_kwargs=Ham_kwargs)
    w,v = np.linalg.eigh(H.matrix)
    amplitude, frequency = calc_fourier(w,v, particle, base, init_base=init_base)
    return amplitude, frequency

# --------------------------------- data preparation for plotting --------------------------------------

def cumulative_sum(input_list):
    # helper function to arange the dataset for plotting
    input_list = list(zip(*input_list))
    result = []
    running_total = [0] * len(input_list[0])
    for sublist in input_list:
        running_total = [x + y for x, y in zip(running_total, sublist)] 
        result.append(running_total[:])  
    return result

def get_pop_lists(DNAstring, J_list, init_base=0, Ham_kwargs={}):
    # returns the average populations for all interactions specified in J_list
    pop_e_list, pop_h_list, pop_ex_list = [],[],[]
    for J in J_list:
        Ham_kwargs['J'] = J
        pop_e, pop_h, pop_ex = static_pop_tot(DNAstring, init_base=init_base, Ham_kwargs=Ham_kwargs)
        pop_e_list.append(pop_e)
        pop_h_list.append(pop_h)
        pop_ex_list.append(pop_ex)
    pop_e_list, pop_h_list, pop_ex_list = cumulative_sum(pop_e_list), cumulative_sum(pop_h_list), cumulative_sum(pop_ex_list)
    zeros = [0]*len(J_list)
    pop_e_list.append(zeros), pop_h_list.append(zeros), pop_ex_list.append(zeros)
    return pop_e_list, pop_h_list, pop_ex_list
