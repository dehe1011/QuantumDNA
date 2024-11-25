"""This module provides utility functions for unit conversion, particularly for
converting physical quantities to different units of measurement.

It includes functions to convert charge separation to dipole moments, convert between
various units, generate all possible conversion factors, and convert all values in a
dictionary from one unit to another.
"""

from itertools import permutations

import numpy as np
import scipy.constants as c

from ..tools import UNITS

__all__ = [
    "UNITS",
    "get_conversion",
    "get_all_conversions",
    "get_conversion_dict",
    "convert_to_debye",
]

# ------------------------------------------------


def convert_to_debye(charge_separation):
    """Converts the charge separation of two particles with elementary charge in
    Angstrom to a dipole moment in Debye.

    Parameters
    ----------
    charge_separation: float
        The electron and hole separation in Angstrom.

    Returns
    -------
    float
        The electrical dipole moment in Debye.
    """

    dipole_moment = c.e * 1e-10 * charge_separation
    return 100 * c.c * dipole_moment / 1e-19


def get_conversion(start_unit, end_unit):
    """Converts a value from the start_unit to the end_unit.

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

    Examples
    --------
    >>> 1/get_conversion("100meV", "rad/fs").
    """

    # Conversion factors to Joule
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

    # Check if the units are defined
    if start_unit not in UNITS:
        raise ValueError(f"Unknown unit: {start_unit}. Defined units: {UNITS}")
    if end_unit not in UNITS:
        raise ValueError(f"Unknown unit: {end_unit}. Defined units: {UNITS}")

    # Return the conversion factor
    return convert_to_Joule[start_unit] / convert_to_Joule[end_unit]


def get_all_conversions():
    """Generates a dictionary with all possible conversions between units.

    Returns
    -------
    Dict[str, float]
        A dictionary where the keys are conversion descriptions
        and the values are the conversion factors.
    """

    conversion_dict = {}
    for start_unit, end_unit in permutations(UNITS, 2):
        conversion = start_unit + "_to_" + end_unit
        conversion_val = get_conversion(start_unit, end_unit)
        conversion_dict[conversion] = conversion_val
    return conversion_dict


def get_conversion_dict(param_dict, start_unit, end_unit):
    """Convert the values in a dictionary from one unit to another.

    Parameters
    ----------
    param_dict : dict
        Dictionary containing the parameters to be converted. Keys are parameter names and values are the parameter values.
    start_unit : str
        The unit of the input values in `param_dict`.
    end_unit : str
        The unit to which the values in `param_dict` should be converted.

    Returns
    -------
    dict
        A new dictionary with the same keys as `param_dict`, but with values converted to `end_unit`.
    """

    conversion_val = get_conversion(start_unit, end_unit)
    converted_dict = {key: value * conversion_val for key, value in param_dict.items()}
    return converted_dict
