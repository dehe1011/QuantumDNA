"""This module provides functions to calculate the average charge separation (dipole
moment) for quantum DNA models."""

import multiprocessing
from functools import partial

import numpy as np
from tqdm import tqdm

from ..dynamics import get_me_solver
from ..model import get_eh_distance
from ..tools import load_json, save_json
from ..utils import convert_to_debye

__all__ = [
    "calc_dipole",
    "calc_dipole_moment",
    "calc_dipole_wrapper",
    "calc_dipole_dict",
    "calc_dipole_moment_dict",
]

# ------------------------------------------------


def calc_dipole(upper_strand, tb_model_name, average=True, **kwargs):
    """Calculates the average charge separation.

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

    Examples
    --------
    >>> calc_dipole("GCG", "ELM")
    2.951734389657976
    """

    kwargs["relax_rate"] = 0
    me_solver = get_me_solver(upper_strand, tb_model_name, **kwargs)
    distance_list = 3.4 * get_eh_distance(me_solver.tb_ham.eh_basis)
    distances = [distance_list @ dm.diag()[1:] for dm in me_solver.get_result()]
    if average:
        return np.mean(distances).real

    return distances


def calc_dipole_moment(upper_strand, tb_model_name, **kwargs):
    """Calculates the dipole moment.

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
    float or List[float]
        The average dipole moment.

    Examples
    --------
    >>> calc_dipole_moment("GCG", "ELM")
    14.177784530660903
    """
    return convert_to_debye(
        calc_dipole(upper_strand, tb_model_name, average=True, **kwargs)
    )


def calc_dipole_wrapper(upper_strand, tb_model_name, lifetime_dict, **kwargs):
    """Calculates the average charge separation in the exciton lifetime.

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


def calc_dipole_dict(tb_model_name, filename, directory, num_cpu=None):
    """Calculates the average charge separation for multiple upper strands using
    multiprocessing.

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
        lifetime_dict, kwargs = load_json(
            "lifetime_" + filename,
            directory,
            load_metadata=True,
        )
    except FileNotFoundError as e:
        print(f"Could not load lifetime_dict: {e}")

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
    save_json(
        dipole_dict,
        kwargs,
        "dipole_" + filename,
        directory,
    )
    return dipole_dict


def calc_dipole_moment_dict(tb_model_name, filename, directory, num_cpu=None):
    """Calculates the dipole moment for multiple upper strands using multiprocessing.

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
        dipole_dict, kwargs = load_json(
            "dipole_" + filename,
            directory,
            load_metadata=True,
        )
        upper_strands = list(dipole_dict.keys())
        dipole_list = list(dipole_dict.values())
    except FileNotFoundError as e:
        print(f"Could not load dipole_dict: {e}")
        try:
            lifetime_dict, kwargs = load_json(
                "lifetime_" + filename,
                directory,
                load_metadata=True,
            )
        except FileNotFoundError as e1:
            print(f"Could not load lifetime_dict: {e1}")

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
                )
            )

    dipole_moment_list = [convert_to_debye(dipole) for dipole in dipole_list]
    dipole_moment_dict = dict(zip(upper_strands, dipole_moment_list))
    save_json(
        dipole_dict,
        kwargs,
        "dipole_moment_" + filename,
        directory,
    )
    return dipole_moment_dict
