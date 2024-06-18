import numpy as np
import qutip as q

__all__ = ['get_site_basis', 'site_basis_converter', 'get_eh_basis', 'eh_basis_converter', 'add_basis_dimension',
          'basis_change', 'global_to_local', 'local_to_global']

# -------------------------------------------- site basis (local) -----------------------------------------------
def get_site_basis(N):
    # returns particle sites |p_ij> forming a basis for the particle Hilbert space
    return [( i//(N//2), i%(N//2) ) for i in range(N)]

def site_basis_converter(N, state):
    if type(state) == str: # '02' for |e_02> -> 2
        i,j = int(state[0]), int(state[1])
        return int( (i*N//2)+j ) 
    elif type(state) == int: # 2 -> [0,2]
        DNA_basis = get_site_basis(N)
        return DNA_basis[state]

# -------------------------------------- eh configuration basis (local) -----------------------------------------

def get_eh_basis(N):
    # returns all e-h configurations |e_ij h_kl> forming a basis for the e-h Hilbert space  
    basis_e, basis_h = get_site_basis(N), get_site_basis(N)
    basis = [(e,h) for e in basis_e for h in basis_h] # ((i,j),(k,l))
    basis = [[x for tpl in tpl_pair for x in tpl] for tpl_pair in basis] # [i,j,k,l]
    return basis

def eh_basis_converter(N, state):
    if type(state) == str: # '0213' for |e_02 h_13> -> 23 
        i,j,k,l = int(state[0]), int(state[1]), int(state[2]), int(state[3])
        return int( ( (i*N//2)+j )*N + ( (k*N//2)+l ) )
    elif type(state) == int: # 23 -> [0,2,1,3]
        DNA_basis = get_eh_basis(N)
        return DNA_basis[state]

# ------------------------------------------- add basis elements ------------------------------------------------

def add_basis_dimension(N, matrix, relaxation=True):
    # includes the system's ground state to account for relaxation by adding one line and column 
    # representing the ground state basis element |0>
    if relaxation==True:
        matrix = matrix.full()
        matrix = np.r_[np.zeros((1,N**2)),matrix]
        matrix = np.c_[np.zeros((N**2+1,1)),matrix]
        matrix = q.Qobj(matrix)
    return matrix

# -------------------------- basis change from local to global basis (eigenbasis) -------------------------------

def basis_change(matrix, states, liouville = False):
    # performs a basis change of the given matrix  
    # states contain the old basis expressed as vector in the new basis
    if liouville: states = np.kron(states, states.conjugate())
    # for open quantum systems the dimenison of the matrices is N**2 instead of N
    return np.matmul(states, np.matmul(matrix, states.conjugate().T))

def global_to_local(matrix, eigs, liouville = False):
    # performs a basis change from the eigenbasis to the site basis
    # eigs: eigenbasis (old) expressed in site basis (new) > transforms from global to local
    states = eigs
    return basis_change(matrix, states, liouville = liouville)

def local_to_global(matrix, eigs, liouville = False):
    # performs a basis change from the site basis to the eigenbasis
    # eigs.conj().T: site basis (old) expressed in eigenbasis (new) > transforms from local to global
    states = eigs.conjugate().T
    return basis_change(matrix, states, liouville = liouville)
