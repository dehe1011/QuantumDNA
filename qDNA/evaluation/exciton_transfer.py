"""
This module provides functions to calculate the average exciton population for quantum DNA models.
"""

# import sys
import multiprocessing
import os
from functools import partial

import numpy as np
from tqdm import tqdm

from qDNA.dynamics import get_me_solver
from qDNA.tools import ROOT_DIR, my_load, my_save

__all__ = [
    "calc_backbone_transfer",
    "calc_average_transfer",
    "calc_exciton_transfer",
    "calc_exciton_transfer_wrapper",
    "calc_exciton_transfer_dict",
]


def calc_average_transfer(tb_sites, me_solver, average=True):
    """
    Calculates the average populations on the given TB sites in a given time period.

    Parameters
    ----------
    tb_sites : List[str]
        The tight-binding sites for which the average population should be calculated.
    me_solver: MESolverType
        Instance of the ME_Solver class.
    average : bool
        Indicates if the exciton_population should be time-averaged.

    Returns
    -------
    dict
        The average populations for each particle.
    """
    average_pop = dict(
        zip(me_solver.tb_ham.particles, [0] * len(me_solver.tb_ham.particles))
    )
    for particle in me_solver.tb_ham.particles:
        avg_pop = np.sum(
            [me_solver.get_pop()[particle + "_" + tb_site] for tb_site in tb_sites],
            axis=0,
        )
        if average:
            avg_pop = np.mean(avg_pop)
        average_pop[particle] = avg_pop
    return average_pop


def calc_backbone_transfer(upper_strand, tb_model_name, **kwargs):
    """
    Calculates the average population of the backbone sites in a given time period.

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
    dict
        The average backbone population for each particle.
    """
    me_solver = get_me_solver(upper_strand, tb_model_name, **kwargs)
    tb_model = me_solver.tb_ham.tb_model
    assert (
        tb_model.tb_model_name[0] == "F"
    ), "Backbone population can only be calculated for Fishbone models"
    upper_backbone_sites = [
        f"(0, {site})" for site in range(tb_model.num_sites_per_strand)
    ]
    lower_backbone_sites = [
        f"({tb_model.num_strands-1}, {site})"
        for site in range(tb_model.num_sites_per_strand)
    ]
    backbone_sites = upper_backbone_sites + lower_backbone_sites
    return calc_average_transfer(backbone_sites, me_solver)


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
    Tuple[dict] or List[Tuple[dict]]
        The average exciton population on the upper and lower strand.

    """
    me_solver = get_me_solver(upper_strand, tb_model_name, **kwargs)
    tb_model = me_solver.tb_ham.tb_model
    assert tb_model.tb_model_name[0] != "F", "Not definded for Fishbone models"
    assert tb_model.tb_model_name[0] != "W", "Not definded for Wire models"

    upper_strand_sites = [
        f"(0, {j})" for j in range(me_solver.tb_model.num_sites_per_strand)
    ]
    lower_strand_sites = [
        f"(1, {j})" for j in range(me_solver.tb_model.num_sites_per_strand)
    ]
    upper_strand_pop = calc_average_transfer(
        upper_strand_sites, me_solver, average=average
    )
    lower_strand_pop = calc_average_transfer(
        lower_strand_sites, me_solver, average=average
    )
    return upper_strand_pop, lower_strand_pop


def calc_backbone_transfer_wrapper(
    upper_strand, tb_model_name, lifetime_dict, **kwargs
):
    """
    Calculates the average backbone population on the upper and lower strand.

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
    dict
        The average backbone population on the upper and lower strand.
    """
    kwargs["t_end"] = lifetime_dict[upper_strand]
    kwargs["t_steps"] = kwargs["t_end"] // 2 + 2
    kwargs["t_unit"] = "fs"
    return calc_backbone_transfer(upper_strand, tb_model_name, **kwargs)


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
    Tuple[dict] or List[Tuple[dict]]
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
            directory=os.path.join(ROOT_DIR, "qDNA", "data", "processed"),
            load_metadata=True,
        )
    except:
        print("Could not load lifetime_dict")

    if not num_cpu:
        num_cpu = multiprocessing.cpu_count() - 1
    upper_strands = list(lifetime_dict.keys())
    partial_calc_exciton_transfer = partial(
        calc_exciton_transfer_wrapper,
        tb_model_name=tb_model_name,
        lifetime_dict=lifetime_dict,
        **kwargs,
    )
    with multiprocessing.Pool(processes=num_cpu) as pool:
        exciton_transfer_list = list(
            tqdm(
                pool.imap(partial_calc_exciton_transfer, upper_strands),
                total=len(upper_strands),
                # file=sys.stdout,
            )
        )

    exciton_transfer_dict = dict(zip(upper_strands, exciton_transfer_list))
    my_save(
        exciton_transfer_dict,
        kwargs,
        "exciton_transfer_" + filename,
        directory=os.path.join(ROOT_DIR, "qDNA", "data", "processed"),
        version_index=False,
    )
    return exciton_transfer_dict
