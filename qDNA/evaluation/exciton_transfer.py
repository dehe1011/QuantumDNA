"""
This module provides functions to calculate the average exciton population for quantum DNA models.
"""

import os
import multiprocessing
from functools import partial
import numpy as np
from tqdm import tqdm

from qDNA.tools import my_save, my_load, ROOT_DIR
from qDNA.dynamics import get_me_solver

__all__ = ["calc_exciton_transfer", "calc_exciton_transfer_wrapper", "calc_exciton_transfer_dict"]

def calc_exciton_transfer(upper_strand, tb_model_name, average=True, **kwargs):
    """
    Calculates the average exciton population on the upper and lower strand.

    Parameters
    ----------
    upper_strand : str
        The upper strand of DNA sequence.
    tb_model_name : str
        The name of the tight-binding model.
    average : bool
        Indicates if the exciton_population should be time-averaged.
    kwargs : dict
        Additional keyword arguments for the master equation solver.

    Returns
    -------
    tuple or List[tuple]
        The average exciton population on the upper and lower strand.
    """
    me_solver = get_me_solver(upper_strand, tb_model_name, **kwargs)
    upper_strand_pop = [me_solver.get_pop()[f'exciton_(0, {j})'] for j in range(me_solver.tb_model.num_sites_per_strand)]
    lower_strand_pop = [me_solver.get_pop()[f'exciton_(1, {j})'] for j in range(me_solver.tb_model.num_sites_per_strand)]
    upper_strand_pop = np.sum( upper_strand_pop, axis=0)
    lower_strand_pop = np.sum( lower_strand_pop, axis=0)
    if average:
        return np.mean(upper_strand_pop), np.mean(lower_strand_pop)
    else:
        return upper_strand_pop, lower_strand_pop

def calc_exciton_transfer_wrapper(upper_strand, tb_model_name, lifetime_dict, **kwargs):
    """
    Calculates the average exciton population on the upper and lower strand.

    Parameters
    ----------
    upper_strand : str
        The upper strand of DNA sequence.
    tb_model_name : str
        The name of the tight-binding model.
    kwargs : dict
        Additional keyword arguments for the master equation solver.

    Returns
    -------
    tuple
        The average exciton population on the upper and lower strand.
    """
    kwargs["t_end"] = lifetime_dict[upper_strand]
    kwargs["t_steps"] = kwargs["t_end"] // 2 + 2
    kwargs["t_unit"] = "fs"
    return calc_exciton_transfer(upper_strand, tb_model_name, **kwargs)


def calc_exciton_transfer_dict(tb_model_name, filename, num_cpu=None):
    """
    Calculates the average exciton population for multiple upper strands using multiprocessing.

    Parameters
    ----------
    tb_model_name : str
        The name of the tight-binding model.
    filename : str
        The filename to load the lifetime dictionary from.
    num_cpu : int, optional
        The number of CPU cores to use. Defaults to the total number of CPUs minus one.

    Returns
    -------
    Dict[str, float]
        Dictionary containing the average exciton population for each upper strand.
    """
    try:
        lifetime_dict, kwargs = my_load(
            "lifetime_" + filename,
            directory=os.path.join(ROOT_DIR, "data", "processed"),
            load_metadata=True,
        )
    except:
        print("Could not load lifetime_dict")

    if not num_cpu:
        num_cpu = multiprocessing.cpu_count() - 1
    upper_strands = list(lifetime_dict.keys())
    partial_calc_dexciton_transfer = partial(
        calc_exciton_transfer_wrapper,
        tb_model_name=tb_model_name,
        lifetime_dict=lifetime_dict,
        **kwargs,
    )
    with multiprocessing.Pool(processes=num_cpu) as pool:
        exciton_transfer_list = list(
            tqdm(
                pool.imap(partial_calc_dexciton_transfer, upper_strands), total=len(upper_strands)
            )
        )

    exciton_transfer_dict = dict(zip(upper_strands, exciton_transfer_list))
    my_save(
        exciton_transfer_dict,
        kwargs,
        "exciton_transfer_" + filename,
        directory=os.path.join(ROOT_DIR, "data", "processed"),
        version_index=False,
    )
    return exciton_transfer_dict