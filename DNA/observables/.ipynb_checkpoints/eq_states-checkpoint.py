import numpy as np
from utils import get_conversion
from DNA.model import global_to_local, local_to_global
from DNA.model import TBHamType
from DNA.dynamics import MESolverType
import scipy.constants as c

__all__ = ['get_therm_eq_state', 'get_deph_eq_state']

# ----------------------------------- equilibrium states --------------------------------------

def get_therm_eq_state(tb_ham: TBHamType, temperature: float) -> np.ndarray:
    """
    Computes the thermal equilibrium state. 
    """
    if temperature == 0: return tb_ham.eigs[:,0] # ground state 
    # thermal equilibrium: frac{e^{- H / k_B T}}{tr(e^{- H / k_B T})}
    eq_values = np.zeros(tb_ham.matrix_dim)
    for i in range(tb_ham.matrix_dim):
        eq_values[i] = np.exp(- tb_ham.eigv[i] * get_conversion(tb_ham.unit, 'J') / (c.k * temperature))
    eq_values = np.diag( eq_values / np.sum(eq_values) )
    eq_state = global_to_local(eq_values, tb_ham.eigs)
    return eq_state

def get_deph_eq_state(me_solver: MESolverType):
    """ 
    Computes the dephasing equilibrium state.
    """
    if me_solver.dephasing_type == 'local': 
        # maximally mixed state (equals inifinite temperature)
        return np.eye(me_solver.tb_ham.matrix_dim) / me_solver.tb_ham.matrix_dim 
    elif me_solver.dephasing_type == 'global':
        loc_init_matrix = me_solver.init_matrix
        glob_init_matrix = local_to_global( loc_init_matrix, me_solver.tb_ham.eigs )
        glob_init_matrix = np.diag(np.diag( glob_init_matrix )) # this cancels all off-diagonal elements 
        loc_init_matrix = global_to_local( glob_init_matrix, me_solver.tb_ham.eigs)
        return loc_init_matrix