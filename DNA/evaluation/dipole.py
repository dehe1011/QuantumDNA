import numpy as np
import multiprocessing
from functools import partial
from itertools import chain 
from tqdm import tqdm
from typing import Dict

from utils import my_save, my_load, get_config
from DNA.model import get_eh_distance
from DNA.dynamics import get_me_solver

# -------------------------------- average charge separation/ dipole moment ------------------------------------------
    
def calc_dipole(upper_strand: str, tb_model_name: str, **kwargs) -> float:
    """
    Calculates the average charge separation.
    """
    kwargs['relax_rate'] = 0
    me_solver = get_me_solver(upper_strand, tb_model_name, **kwargs)
    distance_list = 3.4 * get_eh_distance(me_solver.tb_ham.eh_basis)
    distances = [distance_list @ dm.diag()[1:] for dm in me_solver.get_result()]
    return np.mean(distances).real

def calc_dipole_wrapper(upper_strand: str, tb_model_name: str, lifetime_dict: Dict[str, float], **kwargs) -> float:
    """
    Calculates the average charge separation in the exciton lifetime.
    """
    kwargs['t_end'] = lifetime_dict[upper_strand]
    kwargs['t_steps'] = kwargs['t_end']//2 + 2
    kwargs['t_unit'] = 'fs'
    return calc_dipole(upper_strand, tb_model_name, **kwargs)

def calc_dipole_dict(tb_model_name: str, filename: str, num_cpu: int = None) -> Dict[str, float]:
    try:
        lifetime_dict, kwargs = my_load('lifetime_'+filename+'_version_0', directory = 'stored_data\stored_results', load_metadata=True)
    except:
        print('could not load lifetime_dict from stored_data\stored_results\lifetime_'+filename+'_version_0.json')

    if not num_cpu:
        num_cpu = multiprocessing.cpu_count() - 1
    upper_strands = list(lifetime_dict.keys())
    partial_calc_dipole = partial(calc_dipole_wrapper, tb_model_name=tb_model_name, lifetime_dict = lifetime_dict, **kwargs)
    with multiprocessing.Pool(processes=num_cpu) as pool:
        dipole_list = list(tqdm(pool.imap(partial_calc_dipole, upper_strands), total=len(upper_strands)))
        
    dipole_dict = dict(zip(upper_strands, dipole_list))
    my_save(dipole_dict, kwargs, 'dipole_'+filename, directory = 'stored_data\stored_results')
    if get_config()['verbose']:
        print(f"Saved result at stored_data\stored_results\dipole_"+filename+'.json')
    return dipole_dict
