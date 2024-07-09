import numpy as np
import multiprocessing
from functools import partial 
from itertools import chain 
from tqdm import tqdm
from typing import List, Dict

from utils import my_save, get_config
from DNA.dynamics import get_me_solver

# --------------------------------------- estimated exciton lifetime ----------------------------------------------

def calc_lifetime(upper_strand: str, tb_model_name: str, **kwargs) -> float:
    """
    Calculates lifetime in fs. 
    """
    me_solver = get_me_solver(upper_strand, tb_model_name, **kwargs)
    gs_pop = me_solver.get_groundstate_pop()['groundstate']
    try: 
        _, index = next( (val, i) for i, val in enumerate(gs_pop) if val >= 1-1/np.e)
        lifetime = me_solver.times[index]
        if me_solver.t_unit == 'ps': 
            lifetime *= 1000
        return lifetime
    except: 
        return "no relaxation in the given time"

def calc_lifetime_dict(upper_strands: List[str], tb_model_name: str, filename: str, num_cpu: int = None, **kwargs) -> Dict[str, float]:    
    if not num_cpu:
        num_cpu = multiprocessing.cpu_count() - 1
    partial_calc_lifetime = partial(calc_lifetime, tb_model_name=tb_model_name, **kwargs)
    with multiprocessing.Pool(processes=num_cpu) as pool:
        lifetime_list = list(tqdm(pool.imap(partial_calc_lifetime, upper_strands), total=len(upper_strands)))
        
    lifetime_dict = dict(zip(upper_strands, lifetime_list))
    my_save(lifetime_dict, kwargs, 'lifetime_'+filename, directory = 'stored_data\stored_results', save_excel = False)
    if kwargs.get('verbose') or get_config()['verbose']:
        print(f"Saved result at stored_data\stored_results\lifetime_"+filename+'.json')
    return lifetime_dict
