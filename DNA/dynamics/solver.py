import numpy as np
import qutip as q
from typing import Dict, List, Any, Type
from itertools import permutations

from utils import get_config, get_conversion
from .reduced_dm import get_reduced_dm
from DNA import check_me_kwargs, DNA_Seq
from DNA.model import TBHamType, TBModelType, TB_Ham, TB_Model
from DNA.environment import Lindblad_Diss, LindbladDissType

# Shortcuts
# me: master equation
# diss: dissipator
# t: time
# init: initial
# pop: population
# coh: coherence

__all__ = ['ME_Solver', 'MESolverType', 'get_me_solver']

# ----------------------------------------------------------------

class ME_Solver:
    def __init__(self, tb_ham: TBHamType, lindblad_diss: LindbladDissType, **me_kwargs):
        # check the inputs
        assert isinstance(tb_ham, TB_Ham), "tb_ham must be an instance of the class TB_Ham"
        assert isinstance(lindblad_diss, Lindblad_Diss), "lindblad_diss must be an instance of the class Lindblad_Diss"
        self.me_kwargs = get_config()['me_kwargs_default']
        self.me_kwargs.update(me_kwargs)
        check_me_kwargs(**self.me_kwargs)
        self.verbose = get_config()['verbose']
        if self.verbose: 
            print("Successfully checked all inputs for the ME_Solver instance.")
        
        self.tb_ham = tb_ham
        self.tb_model = self.tb_ham.tb_model 
        self.lindblad_diss = lindblad_diss
        
        self._t_steps = int( self.me_kwargs.get('t_steps') )
        self._t_end = int( self.me_kwargs.get('t_end') )
        self.times = np.linspace(0, self.t_end, self.t_steps ) 
        self.t_unit = self.me_kwargs.get('t_unit')
        assert self.t_steps/self.t_end > 1/2, f"t_end {self.t_end} cannot be sufficiently resolved by t_steps {self.t_steps}. Pleare increase the number of steps or reduce the timespan. Alternative: change the unit of time from fs to ps (the mesolver does not know about the unit, but you do ;) )"
        # if self.t_unit == 'fs':
        #     self.times *= 1e-3 # converts times back from fs to ps
        self.tb_ham.unit = 'rad/'+self.t_unit

        if self.tb_ham.description == '2P':
            self.init_state = (self.me_kwargs.get('init_e_state'), self.me_kwargs.get('init_h_state') )
        if self.tb_ham.description == '1P':
            if self.tb_ham.particles == ['electron']:
                self.init_state = self.me_kwargs.get('init_e_state')
            if self.tb_ham.particles == ['hole']:
                self.init_state = self.me_kwargs.get('init_h_state')
    
        self.init_matrix = self.get_init_matrix()
        self.reset()
            
        if self.verbose:
            print("Successfully initialized the ME_Solver instance.")

    def __vars__(self) -> dict:
        return vars(self)

    def __repr__(self) -> str:
        return f"ME_Solver({self.tb_ham}, {self.lindblad_diss}, {self.me_kwargs})"

    def __eq__(self, other):
        return self.__repr__() == other.__repr__()

    # ------------------------------------------------------------------
            
    @property 
    def t_end(self):
        return self._t_end

    @t_end.setter
    def t_end(self, new_t_end):
        old_t_end = self._t_end
        self._t_end = new_t_end
        if new_t_end != old_t_end:
            self.times = np.linspace(0, self._t_end, self._t_steps)
            self.reset()
            
    @property 
    def t_steps(self):
        return self._t_steps

    @t_steps.setter
    def t_steps(self, new_t_steps):
        old_t_steps = self._t_steps
        self._t_steps = new_t_steps
        if new_t_steps != old_t_steps:
            self.times = np.linspace(0, self._t_end, self._t_steps)
            self.reset()

    # --------------------------------------------------------------------

    def reset(self):
        self.result = []
        if self.tb_ham.description == '2P':
            self.groundstate_pop = {}
            self.pop = {}
            self.coh = {}
            for particle in self.tb_ham.particles:
                vars(self)['result_'+particle] = []

    def get_init_matrix(self) -> q.Qobj:
        if self.tb_ham.description == '2P':
            init_state_idx = self.tb_ham.eh_basis.index(self.init_state)
            if self.tb_ham.relaxation:
                return q.fock_dm(self.tb_ham.matrix_dim, init_state_idx+1)
            else:
                return q.fock_dm(self.tb_ham.matrix_dim, init_state_idx)
        if self.tb_ham.description == '1P':
            init_state_idx = self.tb_ham.tb_basis.index(self.init_state)
            return q.fock_dm(self.tb_ham.matrix_dim, init_state_idx)  

    def get_result(self) -> List[np.ndarray]:
        if not self.result:          
            ham_matrix = q.Qobj( self.tb_ham.matrix )
            result = q.mesolve(ham_matrix, self.init_matrix, self.times, self.lindblad_diss.c_ops, [] , progress_bar = None).states 
            self.result = result
        return self.result

    def get_result_particle(self, particle):
        """
        Function only defined in '2P' description.
        """
        if not self.result: 
            self.get_result()
        if not vars(self)['result_'+particle]:
            vars(self)['result_'+particle] = [get_reduced_dm(dm, particle, self.tb_model.tb_basis) for dm in self.result]
        return vars(self)['result_'+particle]

    def get_pop(self) -> Dict[str, float]:
        assert self.tb_ham.description == '2P', "only available for 2P description"
        if not self.pop:
            ham_matrix = q.Qobj( self.tb_ham.matrix )
            result = q.mesolve(ham_matrix, self.init_matrix, self.times, self.lindblad_diss.c_ops, self.lindblad_diss.pop_ops )
            for particle in self.tb_ham.particles:
                for tb_site in self.tb_ham.tb_basis:
                    self.pop[particle+'_'+tb_site] = result.expect[particle+'_'+tb_site]
        return self.pop 

    def get_coh(self) -> Dict[str, float]:
        assert self.tb_ham.description == '2P', "only available for 2P description"
        if not self.coh:
            ham_matrix = q.Qobj( self.tb_ham.matrix )
            result = q.mesolve(ham_matrix, self.init_matrix, self.times, self.lindblad_diss.c_ops, self.lindblad_diss.coh_ops )
            for particle in self.tb_ham.particles:
                self.coh[particle] = 0
                for tb_site1, tb_site2 in permutations(self.tb_ham.tb_basis, 2):
                    self.coh[particle] += abs( result.expect[particle+'_'+tb_site1+'_'+tb_site2] )
        return self.coh 
                    
    def get_groundstate_pop(self) -> Dict[str, float]:
        assert self.tb_ham.description == '2P', "only available for 2P description"
        if not self.groundstate_pop:
            assert  self.tb_ham.relaxation, "function only defined if relaxation is True, otherwise the groundstate is not populated "
            ham_matrix = q.Qobj( self.tb_ham.matrix )
            result = q.mesolve(ham_matrix, self.init_matrix, self.times, self.lindblad_diss.c_ops, self.lindblad_diss.groundstate_pop_ops )
            self.groundstate_pop['groundstate'] = result.expect['groundstate']
        return self.groundstate_pop 

MESolverType = Type[ME_Solver]

# --------------------------------------------------------------------------------------

def get_me_solver(upper_strand, tb_model_name, **kwargs):
    dna_seq = DNA_Seq(upper_strand, tb_model_name)
    tb_ham = TB_Ham(dna_seq, **kwargs)
    lindblad_diss = Lindblad_Diss(tb_ham, **kwargs)
    me_solver = ME_Solver(tb_ham, lindblad_diss, **kwargs)
    return me_solver
