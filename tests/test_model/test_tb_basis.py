import pytest
import numpy as np
from qDNA.model import get_tb_basis, get_eh_basis, get_eh_distance, get_particle_eh_states, basis_change, global_to_local, local_to_global

@pytest.mark.parametrize("input, expected", [
    ((2, 2), ['(0, 0)', '(0, 1)', '(1, 0)', '(1, 1)'])
])
def test_get_tb_basis(input, expected):
    assert get_tb_basis(input) == expected


@pytest.mark.parametrize("input, expected", [
    ((2, 2), [('(0, 0)', '(0, 0)'), ('(0, 0)', '(0, 1)'), ('(0, 0)', '(1, 0)'), ('(0, 0)', '(1, 1)'),
              ('(0, 1)', '(0, 0)'), ('(0, 1)', '(0, 1)'), ('(0, 1)', '(1, 0)'), ('(0, 1)', '(1, 1)'),
              ('(1, 0)', '(0, 0)'), ('(1, 0)', '(0, 1)'), ('(1, 0)', '(1, 0)'), ('(1, 0)', '(1, 1)'),
              ('(1, 1)', '(0, 0)'), ('(1, 1)', '(0, 1)'), ('(1, 1)', '(1, 0)'), ('(1, 1)', '(1, 1)')])
])
def test_get_eh_basis(input, expected):
    assert get_eh_basis(input) == expected


@pytest.mark.parametrize("input, expected", [
    ([('(0, 0)', '(1, 1)'), ('(1, 0)', '(0, 0)')], np.array([1.41421356, 1.0]))
])
def test_get_eh_distance(input, expected):
    assert np.allclose(get_eh_distance(input), expected)


@pytest.mark.parametrize("ptype, state, basis, expected", [
    ('electron', '(0, 2)', get_tb_basis((1, 3)), [('(0, 2)', '(0, 0)'), ('(0, 2)', '(0, 1)'), ('(0, 2)', '(0, 2)')]),
    ('hole', '(0, 2)', get_tb_basis((1, 3)), [('(0, 0)', '(0, 2)'), ('(0, 1)', '(0, 2)'), ('(0, 2)', '(0, 2)')]),
    ('exciton', '(0, 2)', get_tb_basis((1, 3)), [('(0, 2)', '(0, 2)')])
])
def test_get_particle_eh_states(ptype, state, basis, expected):
    assert get_particle_eh_states(ptype, state, basis) == expected


@pytest.mark.parametrize("matrix, states, expected", [
    (np.array([[0, 1], [1, 0]]), 1/np.sqrt(2) * np.array([[1, 1], [1, -1]]), np.array([[1, 0], [0, -1]]))
])
def test_basis_change(matrix, states, expected):
    assert np.allclose(basis_change(matrix, states), expected)


@pytest.mark.parametrize("matrix, eigs, expected", [
    (np.array([[1, 0], [0, -1]]), 1/np.sqrt(2) * np.array([[1, 1], [1, -1]]), np.array([[0, 1], [1, 0]]))
])
def test_global_to_local(matrix, eigs, expected):
    assert np.allclose(global_to_local(matrix, eigs), expected)


@pytest.mark.parametrize("matrix, eigs, expected", [
    (np.array([[0, 1], [1, 0]]), 1/np.sqrt(2) * np.array([[1, 1], [1, -1]]), np.array([[1, 0], [0, -1]]))
])
def test_local_to_global(matrix, eigs, expected):
    assert np.allclose(local_to_global(matrix, eigs), expected)


if __name__ == "__main__":
    pytest.main()
