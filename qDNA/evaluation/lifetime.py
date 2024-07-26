"""
This module provides functions to calculate the estimated exciton lifetime for quantum DNA models.
"""

import multiprocessing
from functools import partial
import numpy as np
from tqdm import tqdm

from qDNA.tools import my_save, timeit, ROOT_DIR
from qDNA.dynamics import get_me_solver

__all__ = ["calc_lifetime", "calc_lifetime_dict"]

# --------------------------------------- Estimated Exciton Lifetime ------------------------


@timeit
def calc_lifetime(upper_strand, tb_model_name, **kwargs):
    """
    Calculates the exciton lifetime in femtoseconds (fs).

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
    float or str
        The exciton lifetime in femtoseconds, or a message indicating no relaxation in the given time.
    """
    me_solver = get_me_solver(upper_strand, tb_model_name, **kwargs)
    gs_pop = me_solver.get_groundstate_pop()["groundstate"]
    try:
        _, index = next((val, i) for i, val in enumerate(gs_pop) if val >= 1 - 1 / np.e)
        lifetime = me_solver.times[index]
        if me_solver.t_unit == "ps":
            lifetime *= 1000
        return lifetime
    except StopIteration:
        return "no relaxation in the given time"


def calc_lifetime_dict(upper_strands, tb_model_name, filename, num_cpu=None, **kwargs):
    """
    Calculates the exciton lifetime for multiple upper strands using multiprocessing.

    Parameters
    ----------
    upper_strands : List[str]
        List of upper strands of DNA sequences.
    tb_model_name : str
        The name of the tight-binding model.
    filename : str
        The filename to save the lifetime dictionary.
    num_cpu : int, optional
        The number of CPU cores to use. Defaults to the total number of CPUs minus one.
    kwargs : dict
        Additional keyword arguments for the master equation solver.

    Returns
    -------
    Dict[str, float]
        Dictionary containing the exciton lifetime for each upper strand.
    """
    if not num_cpu:
        num_cpu = multiprocessing.cpu_count() - 1
    partial_calc_lifetime = partial(
        calc_lifetime, tb_model_name=tb_model_name, **kwargs
    )
    with multiprocessing.Pool(processes=num_cpu) as pool:
        lifetime_list = list(
            tqdm(
                pool.imap(partial_calc_lifetime, upper_strands),
                total=len(upper_strands),
            )
        )

    lifetime_dict = dict(zip(upper_strands, lifetime_list))
    my_save(
        lifetime_dict,
        kwargs,
        "lifetime_" + filename,
        directory=ROOT_DIR + "/data/processed",
        save_excel=False,
        version_index=False,
    )
    return lifetime_dict
