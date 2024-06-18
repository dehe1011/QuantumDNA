import qutip as q
import numpy as np
from typing import Dict, List, Any, Type
from DNA.model import TBHamType, TBModelType, basis_converter
from .reduced_dm import get_reduced_dm
from DNA.environment import Lindblad_Diss, LindbladDissType
from utils import get_config, get_conversion

class ME_Solver:
    def __init__(self, tb_ham: TBHamType, tb_model: TBModelType, lindblad_diss: LindbladDissType, me_kwargs: Dict = {}):
        self.tb_ham = tb_ham
        self.tb_model = tb_model 
        self.lindblad_diss = lindblad_diss
        self.me_kwargs = get_config()["me_kwargs_default"]
        self.me_kwargs.update(me_kwargs)
        
        self._t_steps = self.me_kwargs.get('t_steps')
        self._t_end = self.me_kwargs.get('t_end') 
        self.times = np.linspace(0, self._t_end, self._t_steps) # in rad/ps
        self._init_matrix = self.get_init_matrix()

        self.result = None
        self.bath_pop = None
        
    @property
    def init_matrix(self):
        return self._init_matrix

    @init_matrix.setter
    def init_matrix(self, new_init_matrix):
        old_init_matrix = self._init_matrix
        self._init_matrix = new_init_matrix
            
    @property 
    def t_end(self):
        return self._t_end

    @t_end.setter
    def t_end(self, new_t_end):
        old_t_end = self._t_end
        self._t_end = new_t_end
        if new_t_end != old_t_end:
            self.times = np.linspace(0, self._t_end, self._t_steps)
            
    @property 
    def t_steps(self):
        return self._t_steps

    @t_steps.setter
    def t_steps(self, new_t_steps):
        old_t_steps = self._t_steps
        self._t_steps = new_t_steps
        if new_t_steps != old_t_steps:
            self.times = np.linspace(0, self._t_end, self._t_steps)

    def get_init_matrix(self) -> q.Qobj:
        if self.tb_ham.particle == 'exciton':
            self.init_state = (self.me_kwargs.get('init_e_state'), self.me_kwargs.get('init_h_state') )
            init_state_num = basis_converter(self.init_state, self.tb_ham.eh_basis)
            if self.tb_ham.relaxation:
                return q.fock_dm(self.tb_ham.matrix_dim, init_state_num+1)
            else:
                return q.fock_dm(self.tb_ham.matrix_dim, init_state_num)
        else:
            self.init_state = self.me_kwargs.get('init_tb_state')
            init_state_num = basis_converter(self.init_state, self.tb_ham.tb_basis)
            return q.fock_dm(self.tb_ham.matrix_dim, init_state_num)  

    def get_result(self) -> List[np.ndarray]:
        if self.result: return self.result
        ham_matrix = q.Qobj( self.tb_ham.matrix )
        result = q.mesolve(ham_matrix, self.init_matrix, self.times, self.lindblad_diss.c_ops, [] , progress_bar = None).states 
        self.times *= get_conversion('rad/ps', self.tb_ham.unit)
        self.result = result
        return result

    def get_result_particle(self, particle):
        if tb_ham.particle == "exciton":
            if not self.result: get_result()
            result_particle = [get_reduced_dm(dm, particle, self.tb_model.tb_basis) for dm in self.result]
            return result_particle
        else: 
            raise("This function can only be called in electron-hole description")
                
    def get_bath_pop(self) -> List[float]:
        if self.bath_pop: return self.bath_pop
        if self.tb_ham.relaxation:
            ham_matrix = q.Qobj( self.tb_ham.matrix )
            result = q.mesolve(ham_matrix, self.init_matrix, self.times, self.lindblad_diss.c_ops, q.fock_dm(self.tb_ham.matrix_dim, 0) )
            self.times *= get_conversion('rad/ps', self.tb_ham.unit)
            self.bath_pop = result.expect[0]
            return self.bath_pop

MESolverType = Type[ME_Solver]
