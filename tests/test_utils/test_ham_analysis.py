import pytest
import numpy as np

from qDNA.utils import (
    calc_average_pop,
    calc_amplitudes,
    calc_frequencies,
    get_pop_fourier,
)

matrix = np.array([[0, 1], [1, 0]])
eigv, eigs = np.linalg.eigh(matrix)


@pytest.mark.parametrize("eigs, state1, state2, expected", [(eigs, 0, 0, 0.5)])
def test_calc_average_pop(eigs, state1, state2, expected):
    assert np.allclose(calc_average_pop(eigs, state1, state2), expected)


@pytest.mark.parametrize(
    "eigs, state1, state2, expected", [(eigs, 0, 0, np.array([0.5]))]
)
def test_calc_amplitudes(eigs, state1, state2, expected):
    assert np.allclose(calc_amplitudes(eigs, state1, state2), expected)


@pytest.mark.parametrize("eigv, expected", [(eigv, np.array([2.0]))])
def test_calc_frequencies(eigv, expected):
    assert np.allclose(calc_frequencies(eigv), expected)


@pytest.mark.parametrize(
    "eigs, eigv, state1, state2, t, expected",
    [(eigs, eigv, 0, 0, 10, 0.7040410309066957)],
)
def test_calc_frequencies(eigs, eigv, state1, state2, t, expected):
    average_pop = calc_average_pop(eigs, state1, state2)
    amplitudes = calc_amplitudes(eigs, state1, state2)
    frequencies = calc_frequencies(eigv)
    pop_fourier = get_pop_fourier(t, average_pop, amplitudes, frequencies)
    assert np.allclose(pop_fourier, expected)
