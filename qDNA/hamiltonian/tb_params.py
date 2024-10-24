"""
This module provides utility functions to save and load tight-binding parameters as JSON files.
By default, the parameters are stored in the "qDNA/data/raw/tb_params" directory.
The module includes functions to handle metadata associated with the parameters, ensuring that the data is well-organized and easily retrievable.
"""

import os
from .. import DATA_DIR
from ..tools import load_json, save_json

__all__ = [
    "save_tb_params",
    "load_tb_params",
    "wrap_save_tb_params",
    "wrap_load_tb_params",
]

# ------------------------------------------------


def save_tb_params(
    tb_params,
    metadata,
    directory,
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
    >>> save_tb_params({'t_AB': 5, 't_AC': 3, 't_BC': -2}, {'source': 'author2024', 'particle': 'particle', 'tb_model_name': 'model'}, "delete_this_folder")
    """

    filename = "_".join(
        [metadata[key] for key in ["source", "particle", "tb_model_name"]]
    )
    save_json(tb_params, metadata, filename, directory)


def load_tb_params(
    metadata,
    directory,
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
    >>> load_tb_params({'source': 'author2024', 'particle': 'particle', 'tb_model_name': 'model'}, "delete_this_folder", load_metadata=False)
    """

    filename = "_".join(
        [metadata[key] for key in ["source", "particle", "tb_model_name"]]
    )
    return load_json(filename, directory, load_metadata=load_metadata)


def wrap_save_tb_params(
    tb_params,
    source,
    particle,
    tb_model_name,
    unit=None,
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
    unit : str, optional
        The unit of the parameters.
    notes : str, optional
        Additional notes to include in the info_dict, by default None.

    Example
    -------
    >>> wrap_save_tb_params({'t_AB': 5, 't_AC': 3, 't_BC': -2}, 'author2024', 'particle', 'model', "delete_this_folder")
    """

    directory = os.path.join(DATA_DIR, "raw", "tb_params")
    metadata = {
        "source": source,
        "particle": particle,
        "tb_model_name": tb_model_name,
        "unit": unit,
        "notes": notes,
    }
    save_tb_params(tb_params, metadata, directory)


def wrap_load_tb_params(
    source,
    particle,
    tb_model_name,
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
    load_metadata : bool, optional
        Whether to load metadata along with the parameters, by default False.

    Returns
    -------
    dict
        Loaded tight-binding parameters.

    Example
    -------
    >>> wrap_load_tb_params('author2024', 'particle', 'model', load_metadata=False)
    """

    directory = os.path.join(DATA_DIR, "raw", "tb_params")
    metadata = {"source": source, "particle": particle, "tb_model_name": tb_model_name}
    return load_tb_params(metadata, directory, load_metadata=load_metadata)
