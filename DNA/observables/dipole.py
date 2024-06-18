import numpy as np
from DNA.model import TB_Model, TB_Ham
from DNA.dynamics import ME_Solver
from DNA.environment import Lindblad_Diss
from DNA.model import get_distance
from utils import my_save, my_load
import multiprocessing
from functools import partial
from itertools import chain 

# -------------------------------- average charge separation/ dipole moment ------------------------------------------
    
def calc_dipole(DNAstring, kwargs={}):
    tb_model = TB_Model('ELM', (2, len(DNAstring[0]) ) )
    tb_ham = TB_Ham(DNAstring, tb_model, Ham_kwargs=kwargs)
    lindblad_diss = Lindblad_Diss(tb_ham,dissipator_kwargs=kwargs)
    me_solver = ME_Solver(tb_ham, tb_model, lindblad_diss, me_kwargs=kwargs)
    distance_list = 3.4*get_distance(tb_ham.eh_basis)
    distance_list = np.array([ 0. ,  3.4,  6.8,  3.4,  6.8, 10.2,  3.4,  0. ,  3.4,  6.8,  3.4,
        6.8,  6.8,  3.4,  0. , 10.2,  6.8,  3.4,  3.4,  6.8, 10.2,  0. ,
        3.4,  6.8,  6.8,  3.4,  6.8,  3.4,  0. ,  3.4, 10.2,  6.8,  3.4,
        6.8,  3.4,  0. ])
    distances = [distance_list @ dm.diag()[1:] for dm in me_solver.get_result()]
    return np.mean(distances).real

def calc_dipole_wrapper(DNAstring, lifetime_dict={}, kwargs={}):
    kwargs['t_end'] = lifetime_dict[''.join(DNAstring)]
    kwargs['relax_rate'] = 0
    return calc_dipole(DNAstring, kwargs = kwargs)

def calc_dipole_dict(DNA_sequence_list, filename, **kwargs):
    lifetime_dict, kwargs1 = my_load('lifetime_'+filename, directory = 'stored_data\stored_results', load_metadata=True)
    num_cpu = multiprocessing.cpu_count()
    partial_calc_dipole = partial(calc_dipole_wrapper, lifetime_dict=lifetime_dict, kwargs=kwargs)
    with multiprocessing.Pool(processes=num_cpu) as pool:
        dipole_list = pool.map(partial_calc_dipole, DNA_sequence_list)
    DNA_sequence_list = [''.join(DNA_sequence) for DNA_sequence in DNA_sequence_list]
    dipole_dict = dict(zip(DNA_sequence_list, dipole_list))
    my_save(dipole_dict, kwargs, 'dipole_'+filename, directory = 'stored_data\stored_results')
    return dipole_dict