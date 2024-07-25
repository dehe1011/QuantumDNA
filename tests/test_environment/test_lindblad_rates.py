import pytest
import numpy as np
import scipy.constants as c

# Importing the functions to be tested
from qDNA.environment import (
    debye_spectral_density,
    ohmic_spectral_density,
    bose_einstein_distrib,
    rate_constant_redfield,
    dephasing_rate,
)


@pytest.mark.parametrize(
    "omega, cutoff_freq, reorg_energy, expected",
    [
        (0, 1.0, 1.0, 0),
        (1.0, 1.0, 1.0, 2.0 * 1.0 * 1.0 / (1.0**2 + 1.0**2)),
        (2.0, 1.0, 1.0, 2.0 * 1.0 * 2.0 * 1.0 / (2.0**2 + 1.0**2)),
    ],
)
def test_debye_spectral_density(omega, cutoff_freq, reorg_energy, expected):
    result = debye_spectral_density(omega, cutoff_freq, reorg_energy)
    assert np.isclose(result, expected)


@pytest.mark.parametrize(
    "omega, cutoff_freq, reorg_energy, exponent, expected",
    [
        (0, 1.0, 1.0, 1.0, 0),
        (1.0, 1.0, 1.0, 1.0, np.pi * 1.0 * 1.0 / 1.0 * np.exp(-1.0 / 1.0)),
        (2.0, 1.0, 1.0, 1.0, np.pi * 1.0 * 2.0 / 1.0 * np.exp(-2.0 / 1.0)),
    ],
)
def test_ohmic_spectral_density(omega, cutoff_freq, reorg_energy, exponent, expected):
    result = ohmic_spectral_density(omega, cutoff_freq, reorg_energy, exponent)
    assert np.isclose(result, expected)


@pytest.mark.parametrize(
    "omega, temperature, expected",
    [
        (0, 300, np.inf),
        (1.0, 300, 1.0 / (np.exp(c.hbar * 1.0 * 1e12 / (c.k * 300)) - 1)),
        (2.0, 300, 1.0 / (np.exp(c.hbar * 2.0 * 1e12 / (c.k * 300)) - 1)),
    ],
)
def test_bose_einstein_distrib(omega, temperature, expected):
    result = bose_einstein_distrib(omega, temperature)
    assert np.isclose(result, expected)


@pytest.mark.parametrize(
    "omega, deph_rate, cutoff_freq, reorg_energy, temperature, spectral_density, exponent, expected",
    [
        (
            0,
            None,
            1.0,
            1.0,
            300,
            "debye",
            None,
            (4 * 1.0 * c.k * 300) / (c.hbar * 1.0 * 1e12),
        ),
        (
            1.0,
            0.5,
            1.0,
            1.0,
            300,
            "debye",
            None,
            (
                2
                * (2.0 * 1.0 * 1.0 / (1.0**2 + 1.0**2))
                * (1 + (1.0 / (np.exp(c.hbar * 1.0 * 1e12 / (c.k * 300)) - 1)))
            ),
        ),
    ],
)
def test_rate_constant_redfield(
    omega,
    deph_rate,
    cutoff_freq,
    reorg_energy,
    temperature,
    spectral_density,
    exponent,
    expected,
):
    result = rate_constant_redfield(
        omega,
        deph_rate,
        cutoff_freq,
        reorg_energy,
        temperature,
        spectral_density,
        exponent,
    )
    assert np.isclose(result, expected)


@pytest.mark.parametrize(
    "cutoff_freq, reorg_energy, temperature, expected",
    [
        (1.0, 1.0, 300, (4 * 1.0 * c.k * 300) / (c.hbar * 1.0 * 1e12)),
    ],
)
def test_dephasing_rate(cutoff_freq, reorg_energy, temperature, expected):
    result = dephasing_rate(cutoff_freq, reorg_energy, temperature)
    assert np.isclose(result, expected)


if __name__ == "__main__":
    pytest.main()
