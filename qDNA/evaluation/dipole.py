"""
This module provides functions to calculate the average charge separation (dipole moment) for quantum DNA models.
"""

# import sys
import multiprocessing
import os
from functools import partial

import numpy as np
from tqdm import tqdm

from qDNA.dynamics import get_me_solver
from qDNA.model import get_eh_distance
from qDNA.tools import ROOT_DIR, my_load, my_save
from qDNA.utils import convert_to_debye

__all__ = [
    "calc_dipole",
    "calc_dipole_moment",
    "calc_dipole_wrapper",
    "calc_dipole_dict",
    "calc_dipole_moment_dict",
]

# -------------------------------- Average Charge Separation/ Dipole Moment ------------------------------------------


def calc_dipole(upper_strand, tb_model_name, average=True, **kwargs):
    """
    Calculates the average charge separation.

    Parameters
    ----------
    upper_strand : str
        The upper strand of DNA sequence.
    tb_model_name : str
        The name of the tight-binding model.
    average : bool
        Indicates if the charge separation should be time-averaged.
    kwargs : dict
        Additional keyword arguments for the master equation solver.

    Returns
    -------
    float or List[float]
        The average charge separation.
    """
    kwargs["relax_rate"] = 0
    me_solver = get_me_solver(upper_strand, tb_model_name, **kwargs)
    distance_list = 3.4 * get_eh_distance(me_solver.tb_ham.eh_basis)
    distances = [distance_list @ dm.diag()[1:] for dm in me_solver.get_result()]
    if average:
        return np.mean(distances).real
    else:
        return distances


def calc_dipole_moment(upper_strand, tb_model_name, average=True, **kwargs):
    """
    Calculates the dipole moment.

    Parameters
    ----------
    upper_strand : str
        The upper strand of DNA sequence.
    tb_model_name : str
        The name of the tight-binding model.
    average : bool
        Indicates if the dipole moment should be time-averaged.
    kwargs : dict
        Additional keyword arguments for the master equation solver.

    Returns
    -------
    float or List[float]
        The average dipole moment.
    """
    return convert_to_debye(
        calc_dipole(upper_strand, tb_model_name, average=average, **kwargs)
    )


def calc_dipole_wrapper(upper_strand, tb_model_name, lifetime_dict, **kwargs):
    """
    Calculates the average charge separation in the exciton lifetime.

    Parameters
    ----------
    upper_strand : str
        The upper strand of DNA sequence.
    tb_model_name : str
        The name of the tight-binding model.
    lifetime_dict : Dict[str, float]
        Dictionary containing lifetimes for each upper strand.
    kwargs : dict
        Additional keyword arguments for the master equation solver.

    Returns
    -------
    float
        The average charge separation during the exciton lifetime.
    """
    kwargs["t_end"] = lifetime_dict[upper_strand]
    kwargs["t_steps"] = kwargs["t_end"] // 2 + 2
    kwargs["t_unit"] = "fs"
    return calc_dipole(upper_strand, tb_model_name, **kwargs)


def calc_dipole_dict(tb_model_name, filename, num_cpu=None):
    """
    Calculates the average charge separation for multiple upper strands using multiprocessing.

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
        Dictionary containing the average charge separation for each upper strand.
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
    partial_calc_dipole = partial(
        calc_dipole_wrapper,
        tb_model_name=tb_model_name,
        lifetime_dict=lifetime_dict,
        **kwargs,
    )
    with multiprocessing.Pool(processes=num_cpu) as pool:
        dipole_list = list(
            tqdm(
                pool.imap(partial_calc_dipole, upper_strands),
                total=len(upper_strands),
                # file=sys.stdout,
            )
        )

    dipole_dict = dict(zip(upper_strands, dipole_list))
    my_save(
        dipole_dict,
        kwargs,
        "dipole_" + filename,
        directory=os.path.join(ROOT_DIR, "qDNA", "data", "processed"),
        version_index=False,
    )
    return dipole_dict


def calc_dipole_moment_dict(tb_model_name, filename, num_cpu=None):
    """
    Calculates the dipole moment for multiple upper strands using multiprocessing.

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
        Dictionary containing the dipole moment for each upper strand.
    """
    try:
        dipole_dict, kwargs = my_load(
            "dipole_" + filename,
            directory=os.path.join(ROOT_DIR, "qDNA", "data", "processed"),
            load_metadata=True,
        )
        upper_strands = list(dipole_dict.keys())
        dipole_list = list(dipole_dict.values())
    except:
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
        partial_calc_dipole = partial(
            calc_dipole_wrapper,
            tb_model_name=tb_model_name,
            lifetime_dict=lifetime_dict,
            **kwargs,
        )
        with multiprocessing.Pool(processes=num_cpu) as pool:
            dipole_list = list(
                tqdm(
                    pool.imap(partial_calc_dipole, upper_strands),
                    total=len(upper_strands),
                    # file=sys.stdout,
                )
            )

    dipole_moment_list = [convert_to_debye(dipole) for dipole in dipole_list]
    dipole_moment_dict = dict(zip(upper_strands, dipole_moment_list))
    my_save(
        dipole_moment_dict,
        kwargs,
        "dipole_moment_" + filename,
        directory=os.path.join(ROOT_DIR, "qDNA", "data", "processed"),
        version_index=False,
    )
    return dipole_moment_dict
