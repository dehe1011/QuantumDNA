import os

from .. import ROOT_DIR, DATA_DIR
from .io_yaml import load_yaml
from .io_json import load_json

# ----------------------------------------------------------------------


def load_defaults(filepath=None):
    if filepath is None:
        filepath = os.path.join(ROOT_DIR, "qDNA", "defaults.yaml")
    return load_yaml(filepath)


def load_lcao_param(param_id, filepath=None):
    if filepath is None:
        filepath = os.path.join(DATA_DIR, "lcao_params", param_id + ".json")
    return load_json(filepath)


def load_tb_model_props(filepath=None):
    if filepath is None:
        filepath = os.path.join(DATA_DIR, "tb_models_props.json")
    return load_json(filepath)


def load_options(filepath=None):
    if filepath is None:
        filepath = os.path.join(DATA_DIR, "options.json")
    return load_json(filepath)


DEFAULTS = load_defaults()
OPTIONS = load_options()

# ----------------------------------------------------------------------
