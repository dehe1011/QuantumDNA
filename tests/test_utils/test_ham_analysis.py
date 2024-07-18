import pytest
import numpy as np
from qDNA import calc_average_pop, calc_amplitudes, calc_frequencies

matrix = np.array([[0, 1], [1, 0]])
eigv, eigs = np.linalg.eigh(matrix)

@pytest.mark.parametrize("eigs, state1, state2, expected", [
    (eigs, 0, 0, 0.5)
])
def test_calc_average_pop(eigs, state1, state2, expected):
    assert np.allclose(calc_average_pop(eigs, state1, state2), expected)

@pytest.mark.parametrize("eigs, state1, state2, expected", [
    (eigs, 0, 0, np.array([0.5]))
])
def test_calc_amplitudes(eigs, state1, state2, expected):
    assert np.allclose(calc_amplitudes(eigs, state1, state2), expected)

@pytest.mark.parametrize("eigv, expected", [
    (eigv, np.array([2.0]))
])
def test_calc_frequencies(eigv, expected):
    assert np.allclose(calc_frequencies(eigv), expected)

if __name__ == "__main__":
    pytest.main()