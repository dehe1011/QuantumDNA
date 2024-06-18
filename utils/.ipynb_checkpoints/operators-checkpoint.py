import qutip as q
import numpy as np 
from itertools import product
from .basis import add_basis_dimension, global_to_local

__all__ = ['H_particle', 'H_eh_config', 'H_exciton', 'reduced_dm_op', 'reduced_dm', 'reduced_dm_electron', 'reduced_dm_hole',
          'reduced_dm_exciton', 'pop_electrons', 'pop_holes', 'pop_excitons', 'coh_electrons', 'coh_holes', 'coh_excitons', 'pop_eig']

# ------------------------------------------ Hamiltonian matrix elements -------------------------------------

def H_particle(N,i,j): 
    # matrix element of the electron/ hole space, i.e., |e_i><e_j|, |h_i><h_j|
    return q.basis(N,i) * q.basis(N,j).dag()

def H_eh_config(N,i,j,k,l):
    # matrix element for the e-h space, i.e., |e_i h_k><e_j h_l|
    return q.tensor( H_particle(N,i,j), H_particle(N,k,l) ) 
    
def H_exciton(N,i,j):  
    # matrix element for the exciton space, i.e., |e_i h_i><e_j h_j|
    return H_eh_config(N,i,j,i,j)

# ------------------------------------------ reduced density matrix  ----------------------------------------

def reduced_dm_op(N,i,j, particle, relaxation):
    # returns all elements of the dm representing the particle transition from position j to i (e.g., |e_i><e_j| )
    # for the electron (hole) this corresponds to tracing out the hole (electron) degrees of freedom
    if particle == 'electron': element = q.Qobj( q.tensor( H_particle(N,i,j), q.qeye(N) ).full() )
    elif particle == 'hole': element = q.Qobj( q.tensor( q.qeye(N), H_particle(N,i,j) ).full() )
    elif particle == 'exciton': element = q.Qobj( H_exciton(N,i,j).full() )
    return add_basis_dimension(N, element, relaxation=relaxation )

def reduced_dm(dm, particle, relaxation):
    # returns the reduced density matrix for electron, hole or exciton of dimension N 
    if relaxation: N = int( np.sqrt( dm.shape[0]-1 ) )
    else: N = int( np.sqrt( dm.shape[0] ) )
    reduced_dm = 0
    for i,j in product(range(N), repeat = 2):
        op = reduced_dm_op(N,i,j, particle, relaxation)
        dm = q.Qobj( dm.full() )
        reduced_dm += q.expect( op, dm ) * H_particle(N,i,j)
    return reduced_dm

def reduced_dm_electron(dm, relaxation):
    return reduced_dm(dm, 'electron', relaxation)

def reduced_dm_hole(dm, relaxation):
    return reduced_dm(dm, 'hole', relaxation)

def reduced_dm_exciton(dm, relaxation):
    return reduced_dm(dm, 'exciton', relaxation)

# -------------------------------- observables for population and coherence -------------------------

def pop_electrons(N, i, relaxation): 
    # observable for the electron population of site i
    return reduced_dm_op(N,i,i, 'electron', relaxation)
    
def pop_holes(N, i, relaxation):  
    # observable for the hole population of site i
    return reduced_dm_op(N,i,i, 'hole', relaxation)
    
def pop_excitons(N, i, relaxation): 
    # observable for the exciton population of site i
    return reduced_dm_op(N,i,i, 'exciton', relaxation)

def coh_electrons(N, i, j, relaxation): 
    # observable for the electron coherence between sites i and j
    return reduced_dm_op(N,i,j, 'electron', relaxation)
    
def coh_holes(N, i, j, relaxation):  
    # observable for the hole coherence between sites i and j
    return reduced_dm_op(N,i,j, 'hole', relaxation)
    
def coh_excitons(N, i, j, relaxation): 
    # observable for the exciton coherence between sites i and j
    return reduced_dm_op(N,i,j, 'exciton', relaxation)

def pop_eig(N, i, eigs, relaxation):
    # observable for the population of eigenstate i 
    pop_eig = add_basis_dimension(N, q.fock_dm(N**2,i), relaxation=relaxation )
    pop_eig = global_to_local(pop_eig.full(), eigs)
    return q.Qobj(pop_eig)
