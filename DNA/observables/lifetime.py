import numpy as np
from DNA.model import TB_Model, TB_Ham
from DNA.dynamics import ME_Solver
from DNA.environment import Lindblad_Diss
from utils import my_save, my_load
import multiprocessing
from functools import partial
from itertools import chain 

# --------------------------------------- estimated exciton lifetime ----------------------------------------------

def calc_lifetime(DNAstring, kwargs={}):
    tb_model = TB_Model('ELM', (2, len(DNAstring[0]) ) )
    tb_ham = TB_Ham(DNAstring, tb_model, Ham_kwargs=kwargs)
    lindblad_diss = Lindblad_Diss(tb_ham,dissipator_kwargs=kwargs)
    me_solver = ME_Solver(tb_ham, tb_model, lindblad_diss, me_kwargs=kwargs)
    gs_pop = me_solver.get_bath_pop()
    _, index = next( (val, i) for i, val in enumerate(gs_pop) if val >= 1-1/np.e)
    return me_solver.times[index]

def calc_lifetime_dict(DNA_sequence_list, filename, **kwargs):    
    num_cpu = multiprocessing.cpu_count()
    partial_calc_lifetime = partial(calc_lifetime, kwargs = kwargs)
    with multiprocessing.Pool(processes=num_cpu) as pool:
        lifetime_list = pool.map(partial_calc_lifetime, DNA_sequence_list)
    DNA_sequence_list = [''.join(DNA_sequence) for DNA_sequence in DNA_sequence_list]
    lifetime_dict = dict(zip(DNA_sequence_list, lifetime_list))
    my_save(lifetime_dict, kwargs, 'lifetime_'+filename, directory = 'stored_data\stored_results', save_excel = False)
    return lifetime_dict

# --------------------------------- relative distribution of DNA nucleobases ------------------------------------

def base_counter(my_dict):
    # counts the relative amount of nucleobases in the given dict
    A, T, G, C = 0, 0, 0, 0
    A_list, T_list, G_list, C_list = [],[],[],[]
    val_list = []
    for sequence, val in my_dict.items():
        val_list.append(val)
        N = len(sequence)
        A += sequence.count('A')
        summe=A+T+G+C
        number_sequences = summe//N
        A_list.append(A/summe)
        T_list.append(T/summe)
        G_list.append(G/summe)
        C_list.append(C/summe)
    return val_list, np.array(A_list),np.array(T_list),np.array(G_list),np.array(C_list)
