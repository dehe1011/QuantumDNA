import pytest
import numpy as np
from qDNA import calc_trace_distance, calc_purity, calc_coherence, calc_ipr_dm

dm_1 = np.array([[0.5, 0.5j], [-0.5j, 0.5]])
dm_2 = np.array([[0.5, 0], [0, 0.5]])

@pytest.mark.parametrize("dm1, dm2, expected", [
    (dm_1, dm_2, 0.5)
])
def test_calc_trace_distance(dm1, dm2, expected):
    assert np.allclose(calc_trace_distance(dm1, dm2), expected)

@pytest.mark.parametrize("dm, expected", [
    (dm_1, 1.0),
    (dm_2, 0.5)
])
def test_calc_purity(dm, expected):
    assert np.allclose(calc_purity(dm), expected)

@pytest.mark.parametrize("dm, expected", [
    (dm_1, 1.0),
    (dm_2, 0.0)
])
def test_calc_coherence(dm, expected):
    assert np.allclose(calc_coherence(dm), expected)

@pytest.mark.parametrize("dm, expected", [
    (dm_1, 2.0),
    (dm_2, 1.0)
])
def test_calc_ipr_dm(dm, expected):
    assert np.allclose(calc_ipr_dm(dm), expected)

if __name__ == "__main__":
    pytest.main()