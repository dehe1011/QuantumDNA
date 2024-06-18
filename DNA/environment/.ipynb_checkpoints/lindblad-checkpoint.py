import numpy as np
import qutip as q
from DNA.model import TBHamType, add_basis_dimension, add_basis_dimension, global_to_local, basis_converter
from DNA.observables import get_pop_particle, get_eh_observable
from typing import List, Type

from utils import get_conversion, get_config
from .lindblad_rates import rate_constant_redfield
from itertools import product, permutations, chain    

# -------------------------------------- Lindblad operators -------------------------------------

# collect available Lindblad dissipators 
class Lindblad_Diss:
    def __init__(self, tb_ham: TBHamType, dissipator_kwargs={}):
        self.tb_ham = tb_ham

        self.dissipator_kwargs = get_config()["dissipator_kwargs_default"]
        self.dissipator_kwargs.update(dissipator_kwargs)
        
        self.loc_deph_rate = self.dissipator_kwargs.get('loc_deph_rate') # rad/ps
        self.glob_deph_rate = self.dissipator_kwargs.get('glob_deph_rate') # rad/ps
        self.relax_rate = self.dissipator_kwargs.get('relax_rate') # rad/ps
        
        self.loc_therm = self.dissipator_kwargs.get('loc_therm')
        self.glob_therm = self.dissipator_kwargs.get('glob_therm')
        self.deph_rate = self.dissipator_kwargs.get('deph_rate') # in rad/ps
        self.cutoff_freq = self.dissipator_kwargs.get('cutoff_freq') # in rad/ps
        self.reorg_energy = self.dissipator_kwargs.get('reorg_energy') # in rad/ps
        self.temperature = self.dissipator_kwargs.get('temperature') # in K
        self.spectral_density = self.dissipator_kwargs.get('spectral_density')
        self.exponent = self.dissipator_kwargs.get('exponent')

        self.relax_ops = self.get_relax_ops()
        self.deph_ops = self.get_deph_ops()
        self.therm_ops = self.get_therm_ops()
        self.c_ops = self.relax_ops + self.deph_ops + self.therm_ops

    def get_deph_ops(self):
        if self.loc_deph_rate and self.glob_deph_rate:
            print("Error: dephasing must be either global or local.")
        if self.loc_deph_rate:
            if self.tb_ham.particle == 'exciton':
                return loc_deph_ops(self.tb_ham.tb_basis, self.loc_deph_rate, self.tb_ham.relaxation)
            if self.tb_ham.particle in ('electron', 'hole'):
                return loc_deph_p_ops(self.tb_ham.tb_basis, self.loc_deph_rate)
        if self.glob_deph_rate:
            if self.tb_ham.particle == 'exciton':
                return glob_deph_ops(self.tb_ham.eigs, self.glob_deph_rate, self.tb_ham.relaxation)
            if self.tb_ham.particle in ('electron', 'hole'):
                return glob_deph_p_ops(self.tb_ham.eigs, self.glob_deph_rate)
        else: return []

    def get_relax_ops(self):
        if self.tb_ham.relaxation and self.relax_rate: 
            return [ np.sqrt(self.relax_rate) * get_relax_op(self.tb_ham.tb_basis, tb_site) for tb_site in self.tb_ham.tb_basis]
        else: return []

    def get_therm_ops(self):
        if self.loc_therm and self.glob_therm:
            print("Error: thermalizing must be either global or local.")
        eigv = self.tb_ham.eigv * get_conversion(self.tb_ham.unit, 'rad/ps')
        if self.loc_therm:
            return get_loc_therm_ops(eigv, self.tb_ham.eigs, self.tb_ham.relaxation, deph_rate=self.deph_rate, cutoff_freq=self.cutoff_freq, reorg_energy=self.reorg_energy, temperature=self.temperature, spectral_density=self.spectral_density, exponent=self.exponent)
        if self.glob_therm:
            return get_glob_therm_ops(eigv, self.tb_ham.eigs, self.tb_ham.relaxation, deph_rate=self.deph_rate, cutoff_freq=self.cutoff_freq, reorg_energy=self.reorg_energy, temperature=self.temperature, spectral_density=self.spectral_density, exponent=self.exponent)
        else: return []

    def __vars__(self) -> dict:
        return vars(self)

LindbladDissType = Type[Lindblad_Diss]

# -------------------------------------- relaxation operator ------------------------------------

def get_relax_op(tb_basis, tb_site):
    """
    Annihilation operator of an exciton on a given tight-binding site. Relaxation of the DNA to its ground state.
    """
    tb_site = basis_converter(tb_site, tb_basis)
    num_sites = len(tb_basis)
    relax_op = np.zeros((num_sites**2+1,num_sites**2+1))
    relax_op[0,1+tb_site*(num_sites+1)] = 1
    return q.Qobj(relax_op)

# ------------------------------------- thermalizing operators -----------------------------------

def get_glob_therm_op(eigs, eigenstate_i, eigenstate_j, relaxation, matrix_dim):
    # op(e_i,e_j) = |e_i><e_j| 
    op = np.zeros((matrix_dim, matrix_dim), dtype=complex)
    op[eigenstate_i][eigenstate_j] = 1
    op = global_to_local(op, eigs)
    if relaxation:
        op = add_basis_dimension(op)
    return q.Qobj(op)

def get_glob_therm_ops(eigv, eigs, relaxation, deph_rate=7, cutoff_freq=20, reorg_energy=1, temperature=300, spectral_density='debye', exponent=1):
    # c_ops = [rate(omega_i-omega_j)*op(e_i,e_j) for i,j]
    matrix_dim = eigs.shape[0]
    c_ops = []
    # the sum runs over all different eigenstates
    for eigenstate_i, eigenstate_j in permutations(range(matrix_dim), 2):
        omega_i, omega_j = eigv[eigenstate_i], eigv[eigenstate_j]
        lind_rate = rate_constant_redfield( (omega_i - omega_j), deph_rate, cutoff_freq, reorg_energy, temperature, spectral_density, exponent )
        lind_op = get_glob_therm_op(eigs, eigenstate_i, eigenstate_j, relaxation, matrix_dim)
        c_ops.append( np.sqrt(lind_rate) * lind_op )
    return c_ops

# local thermalising 
def get_loc_therm_op(eigv, eigs, unique, site_m, relaxation, matrix_dim):
    # op(unique,m) = sum_{i,j} <e_j|m> <m|e_i> |e_j><e_i| \delta_{\omega_i - omega_j, unique}
    # note that if you sum over all unique frequency gaps you obtain the local state |m><m|, but all of them have different frequencies!
    op = np.zeros((matrix_dim, matrix_dim), dtype=complex)
    for i, j in product(range(matrix_dim), repeat=2):
        omega_i, omega_j = eigv[i], eigv[j]
        state_i, state_j = eigs[:, i], eigs[:, j]
        if omega_i - omega_j == unique:
            op += (state_j[site_m].conjugate() * state_i[site_m] * np.outer(state_j, state_i))
    if relaxation: 
        op = add_basis_dimension(op)
    return q.Qobj(op)

def get_loc_therm_ops(eigv, eigs, relaxation, deph_rate=7, cutoff_freq=20, reorg_energy=1, temperature=300, spectral_density='debye', exponent=1):
    # c_ops = [ rate(unique)*op(unique,m) for unique,m in freq_gaps, local_sites]
    matrix_dim = len(eigv)
    c_ops = []
    gaps = eigv.reshape(matrix_dim, 1) - eigv # matrix that contains all possible eigenenergy differences
    unique = np.unique(gaps.flatten()) 
    # the sum runs over all available frequency gaps of the system and all local sites 
    for unique, site_m in product(unique, range(matrix_dim)):
        # for each unique frequency difference the rate is different 
        lind_rate = rate_constant_redfield( unique, deph_rate, cutoff_freq, reorg_energy, temperature, spectral_density, exponent)
        lind_op = get_loc_therm_op(eigv, eigs, unique, site_m, relaxation, matrix_dim)
        c_ops.append( np.sqrt(lind_rate) * lind_op )
    return c_ops

# ---------------------- local and global dephasing operators (for e-h and particle description) -------------------

def loc_deph_ops(tb_basis, dephasing_rate, relaxation):
    c_ops = []
    for tb_site in tb_basis:
        op_electron = get_pop_particle(tb_basis, 'electron', tb_site)
        op_hole = get_pop_particle(tb_basis, 'hole', tb_site)
        if relaxation: 
            op_electron = add_basis_dimension(op_electron)
            op_hole = add_basis_dimension(op_hole)
        c_ops.append(np.sqrt(dephasing_rate) * q.Qobj(op_electron) )
        c_ops.append(np.sqrt(dephasing_rate) * q.Qobj(op_hole) )
    return c_ops

def glob_deph_ops(eigs, dephasing_rate, relaxation):
    num_eigenstates = eigs.shape[0]
    c_ops = []
    for i in range(num_eigenstates):
        matrix = np.zeros((num_eigenstates, num_eigenstates))
        matrix[i,i] = 1
        op = global_to_local(matrix, eigs)
        if relaxation:
            op = add_basis_dimension(op)
        c_ops.append(np.sqrt(dephasing_rate) * q.Qobj(op) )
    return c_ops

def loc_deph_p_ops(tb_basis, dephasing_rate): 
    num_sites = len(tb_basis)
    return [ np.sqrt(dephasing_rate) * q.fock_dm(num_sites,i) for i in range(num_sites) ]
    
def glob_deph_p_ops(eigs, dephasing_rate):
    num_eigenstates = eigs.shape[0]
    return [np.sqrt(dephasing_rate) * q.Qobj(global_to_local( q.fock_dm(num_eigenstates,i).full(), eigs)) for i in range(num_eigenstates) ]

# ----------------------------------- Lindblad observables -----------------------------------

def get_e_ops(tb_basis: List[str], relaxation: bool):
    pop_dict, coh_dict, bath_dict = {},{},{}
    for particle in PARTICLES:
        for tb_site1, tb_site2 in product(tb_basis, repeat=2):
            observable = get_eh_observable(tb_basis, particle, tb_site1, tb_site2)
            if relaxation: 
                observable = add_basis_dimension(observable)
            if tb_site1 == tb_site2:
                pop_dict[particle+'_'+tb_site1+'_'+tb_site2] = q.Qobj(observable)
            else:
                coh_dict[particle+'_'+tb_site1+'_'+tb_site2] = q.Qobj(observable)
    if relaxation:
        groundstate_pop_dict["groundstate"] = q.fock_dm( len(tb_basis)**2+1, 0 )
    return pop_dict, coh_dict, groundstate_pop_dict
