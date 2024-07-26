"""
This module is taken from quantum_HEOM (github.com/jwa7/quantum_HEOM), J.W. Abbott, 2022, DOI: 10.5281/zenodo.7230160.

It provides functions for calculating bath spectral densities and Lindblad rates.
"""

import numpy as np
import scipy.constants as c
from qDNA.tools import get_config

SPECTRAL_DENSITIES = get_config()["SPECTRAL_DENSITIES"]

# --------------------------- Bath Spectral Densities --------------------------------------


def debye_spectral_density(omega, cutoff_freq, reorg_energy):
    """
    Calculates the Debye spectral density.

    Parameters
    ----------
    omega : float
        Frequency.
    cutoff_freq : float
        Cutoff frequency.
    reorg_energy : float
        Reorganization energy.

    Returns
    -------
    float
        Debye spectral density.
    """
    if omega <= 0:
        return 0
    return 2 * reorg_energy * omega * cutoff_freq / (omega**2 + cutoff_freq**2)


def ohmic_spectral_density(omega, cutoff_freq, reorg_energy, exponent):
    """
    Calculates the Ohmic spectral density.

    Parameters
    ----------
    omega : float
        Frequency.
    cutoff_freq : float
        Cutoff frequency.
    reorg_energy : float
        Reorganization energy.
    exponent : float
        Specifies a sub- or superohmic bath.

    Returns
    -------
    float
        Ohmic spectral density.
    """
    if omega <= 0:
        return 0
    return (np.pi * reorg_energy * omega / cutoff_freq) * np.exp(-omega / cutoff_freq)


# ----------------------------- Lindblad Rates -------------------------------------------


def bose_einstein_distrib(omega, temperature):
    """
    Calculates the Bose-Einstein distribution.

    Parameters
    ----------
    omega : float
        Frequency.
    temperature : float
        Temperature.

    Returns
    -------
    float
        Bose-Einstein distribution.
    """
    return 1.0 / (np.exp(c.hbar * omega * 1e12 / (c.k * temperature)) - 1)


def rate_constant_redfield(
    omega,
    deph_rate,
    cutoff_freq,
    reorg_energy,
    temperature,
    spectral_density,
    exponent=None,
):
    """
    Calculates the Redfield rate constant.

    Parameters
    ----------
    omega : float
        Frequency.
    deph_rate : float or None
        Dephasing rate. If None, it will be calculated.
    cutoff_freq : float
        Cutoff frequency.
    reorg_energy : float
        Reorganization energy.
    temperature : float
        Temperature.
    spectral_density : str
        Spectral density type ('debye' or 'ohmic').
    exponent : float, optional
        Exponent for the Ohmic spectral density.

    Returns
    -------
    float
        Redfield rate constant.
    """
    if omega == 0:
        if deph_rate is None:
            deph_rate = dephasing_rate(cutoff_freq, reorg_energy, temperature)
        return deph_rate

    if spectral_density == "debye":
        spec_omega_ij = debye_spectral_density(omega, cutoff_freq, reorg_energy)
        spec_omega_ji = debye_spectral_density(-omega, cutoff_freq, reorg_energy)
    elif spectral_density == "ohmic":
        spec_omega_ij = ohmic_spectral_density(
            omega, cutoff_freq, reorg_energy, exponent
        )
        spec_omega_ji = ohmic_spectral_density(
            -omega, cutoff_freq, reorg_energy, exponent
        )

    n_omega_ij = bose_einstein_distrib(omega, temperature)
    n_omega_ji = bose_einstein_distrib(-omega, temperature)
    return 2 * ((spec_omega_ij * (1 + n_omega_ij)) + (spec_omega_ji * n_omega_ji))


def dephasing_rate(cutoff_freq, reorg_energy, temperature):
    """
    Calculates the dephasing rate in the limit of the Redfield rate equation as the frequency approaches zero.

    Parameters
    ----------
    cutoff_freq : float
        Cutoff frequency.
    reorg_energy : float
        Reorganization energy.
    temperature : float
        Temperature.

    Returns
    -------
    float
        Dephasing rate.
    """
    return (4 * reorg_energy * c.k * temperature) / (c.hbar * cutoff_freq * 1e12)
