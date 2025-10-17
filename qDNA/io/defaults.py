import os

from .. import ROOT_DIR, DATA_DIR
from .io_yaml import load_yaml
from .io_json import load_json

# ----------------------------------------------------------------------


def load_defaults(filepath=None):
    """
    Load default settings from a YAML file.

    Parameters
    ----------
    filepath : str, optional
        Path to the YAML file containing default settings. If None,
        the function loads the defaults from 'defaults.yaml' in the
        qDNA directory.

    Returns
    -------
    dict
        Parsed contents of the YAML file as a dictionary.
    """

    if filepath is None:
        filepath = os.path.join(ROOT_DIR, "qDNA", "defaults.yaml")
    return load_yaml(filepath)


def load_lcao_param(param_id, filepath=None):
    """
    Load Linear Combination of Atomic Orbitals (LCAO) parameters from a JSON file.

    Parameters
    ----------
    param_id : str
        Identifier for the parameter file to load.
    filepath : str, optional
        Path to the JSON file. If not provided, the default path is constructed
        using `DATA_DIR` and the `param_id`.

    Returns
    -------
    dict
        Parsed JSON data containing the LCAO parameters.
    """

    if filepath is None:
        filepath = os.path.join(DATA_DIR, "lcao_params", param_id + ".json")
    return load_json(filepath)


def load_tb_model_props(filepath=None):
    """
    Load tight-binding model properties from a JSON file.

    Parameters
    ----------
    filepath : str, optional
        Path to the JSON file containing the model properties. If None,
        defaults to the file "tb_models_props.json" in the DATA_DIR directory.

    Returns
    -------
    dict
        Dictionary containing the loaded model properties.
    """

    if filepath is None:
        filepath = os.path.join(DATA_DIR, "tb_models_props.json")
    return load_json(filepath)


def load_options(filepath=None):
    """
    Load options from a JSON file.

    Parameters
    ----------
    filepath : str, optional
        Path to the JSON file containing options. If None, defaults to
        'options.json' in the DATA_DIR directory.

    Returns
    -------
    dict
        Dictionary containing the loaded options.
    """

    if filepath is None:
        filepath = os.path.join(DATA_DIR, "options.json")
    return load_json(filepath)


DEFAULTS = load_defaults()
OPTIONS = load_options()

# ----------------------------------------------------------------------
