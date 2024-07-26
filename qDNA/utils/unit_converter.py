"""
Module for unit conversion utilities.
"""

from itertools import permutations
import numpy as np
import scipy.constants as c

from qDNA.tools import get_config

__all__ = ["UNITS", "get_conversion", "get_all_conversions", "get_conversion_dict"]

UNITS = get_config()["UNITS"]


def get_conversion(start_unit, end_unit):
    """
    Converts a value from the start_unit to the end_unit.

    Parameters
    ----------
    start_unit : str
        The unit to convert from.
    end_unit : str
        The unit to convert to.

    Returns
    -------
    float
        The conversion factor between the two units.

    Notes
    -----
    DNA parameters are saved in 100 meV. We mainly use the conversion '100meV_to_rad/fs'.
    Example: real_units = 1 / get_conversion('100meV', 'rad/fs').
    """
    convert_to_Joule = {
        "J": 1,
        "eV": c.e,
        "meV": 1e-3 * c.e,
        "100meV": 100 * 1e-3 * c.e,
        "rad/fs": 1 / 1e-15 * c.hbar,
        "rad/ps": 1 / 1e-12 * c.hbar,
        "1/cm": c.c / 1e-2 * 2 * np.pi * c.hbar,
        "K": c.k,
        "1/fs": 2 * np.pi / 1e-15 * c.hbar,
        "1/ps": 2 * np.pi / 1e-12 * c.hbar,
    }
    if start_unit not in UNITS:
        raise ValueError(f"Unknown unit: {start_unit}. Defined units: {UNITS}")
    if end_unit not in UNITS:
        raise ValueError(f"Unknown unit: {end_unit}. Defined units: {UNITS}")
    return convert_to_Joule[start_unit] / convert_to_Joule[end_unit]


def get_all_conversions():
    """
    Generates a dictionary with all possible conversions between units.

    Returns
    -------
    Dict[str, float]
        A dictionary where the keys are conversion descriptions (e.g., '100meV_to_rad/fs')
        and the values are the conversion factors.
    """
    conversion_dict = {}
    for start_unit, end_unit in permutations(UNITS, 2):
        conversion = start_unit + "_to_" + end_unit
        conversion_val = get_conversion(start_unit, end_unit)
        conversion_dict[conversion] = conversion_val
    return conversion_dict


def get_conversion_dict(param_dict, start_unit, end_unit):
    """
    Converts the values in param_dict from start_unit to end_unit.

    Parameters
    ----------
    param_dict : Dict[Any, float]
        Dictionary containing parameters to convert.
    start_unit : str
        The unit to convert from.
    end_unit : str
        The unit to convert to.

    Returns
    -------
    Dict[Any, float]
        The parameter dictionary with converted values.
    """
    conversion_val = get_conversion(start_unit, end_unit)
    converted_dict = {key: value * conversion_val for key, value in param_dict.items()}
    return converted_dict
