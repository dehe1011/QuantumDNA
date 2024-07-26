"""
This module provides functions to save and load tight-binding parameters for quantum DNA models.
"""

import os
from copy import deepcopy
from qDNA.tools import ROOT_DIR, my_save, my_load

__all__ = [
    "save_tb_params",
    "load_tb_params",
    "wrap_save_tb_params",
    "wrap_load_tb_params",
]


def save_tb_params(
    tb_param_dict,
    info_dict,
    directory=os.path.join(ROOT_DIR, "data", "raw", "tb_params"),
    notes=None,
):
    """
    Save tight-binding parameters to a file.

    Parameters
    ----------
    tb_param_dict : dict
        Dictionary containing the tight-binding parameters.
    info_dict : dict
        Dictionary with metadata, e.g., 'source', 'particle', and 'tb_model_name'.
    directory : str, optional
        Directory to save the file, by default 'data/raw/tb_params'.
    notes : str, optional
        Additional notes to include in the info_dict, by default None.

    Example
    -------
    >>> save_tb_params(
            {'couple': 5, 'bff': 3, 'enemies': -2, 'very_cool': 3.5},
            {'source': 'author2024', 'particle': 'particle', 'tb_model_name': 'WM'},
            directory='stored_data/my_params',
            notes='here you can add some notes.'
        )
    """
    filename = "_".join(info_dict.values())
    info_dict_notes = deepcopy(info_dict)
    if notes:
        info_dict_notes["notes"] = notes
    my_save(
        tb_param_dict,
        info_dict_notes,
        filename,
        directory=directory,
        save_excel=False,
        version_index=False,
    )


def load_tb_params(
    info_dict,
    directory=os.path.join(ROOT_DIR, "data", "raw", "tb_params"),
    load_metadata=False,
):
    """
    Load tight-binding parameters from a file.

    Parameters
    ----------
    info_dict : dict
        Dictionary with metadata, e.g., 'source', 'particle', and 'tb_model_name'.
    directory : str, optional
        Directory to load the file from, by default 'data/raw/tb_params'.
    load_metadata : bool, optional
        Whether to load metadata along with the parameters, by default False.

    Returns
    -------
    dict
        Loaded tight-binding parameters.

    Example
    -------
    >>> load_tb_params(
            {'source': 'author2024', 'particle': 'particle', 'tb_model_name': 'WM'},
            directory='stored_data/my_params',
            load_metadata=True
        )
    """
    filename = "_".join(info_dict.values())
    return my_load(filename, load_metadata=load_metadata, directory=directory)


def wrap_save_tb_params(
    tb_param_dict,
    source,
    particle,
    tb_model_name,
    unit,
    directory=os.path.join(ROOT_DIR, "data", "raw", "tb_params"),
    notes=None,
):
    """
    Wrapper function for save_tb_params().

    Parameters
    ----------
    tb_param_dict : dict
        Dictionary containing the tight-binding parameters.
    source : str
        Source of the parameters, e.g., 'Hawke2010'.
    particle : str
        Type of particle, e.g., 'electron', 'hole', 'exciton'.
    tb_model_name : str
        Name of the tight-binding model.
    unit : str
        The unit of the parameters.
    directory : str, optional
        Directory to save the file, by default 'data/raw/tb_params'.
    notes : str, optional
        Additional notes to include in the info_dict, by default None.
    """
    info_dict = {
        "source": source,
        "particle": particle,
        "tb_model_name": tb_model_name,
        "unit": unit,
    }
    save_tb_params(tb_param_dict, info_dict, directory=directory, notes=notes)


def wrap_load_tb_params(
    source,
    particle,
    tb_model_name,
    directory=os.path.join(ROOT_DIR, "data", "raw", "tb_params"),
    load_metadata=False,
):
    """
    Wrapper function for load_tb_params().

    Parameters
    ----------
    source : str
        Source of the parameters, e.g., 'Hawke2010'.
    particle : str
        Type of particle, e.g., 'electron', 'hole', 'exciton'.
    tb_model_name : str
        Name of the tight-binding model.
    directory : str, optional
        Directory to load the file from, by default 'data/raw/tb_params'.
    load_metadata : bool, optional
        Whether to load metadata along with the parameters, by default False.

    Returns
    -------
    dict
        Loaded tight-binding parameters.
    """
    info_dict = {"source": source, "particle": particle, "tb_model_name": tb_model_name}
    return load_tb_params(info_dict, directory=directory, load_metadata=load_metadata)
