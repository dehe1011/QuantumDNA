import unittest
import numpy as np
from itertools import product

from DNA.model import (
    basis_change, global_to_local, local_to_global
)

class TestBasisChange(unittest.TestCase):
    def test_basis_change(self):
        matrix = np.array([[0, 1], [1, 0]])
        states = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]])
        expected = np.array([[1, 0], [0, -1]])
        np.allclose(basis_change(matrix, states), expected)

    def test_global_to_local(self):
        matrix = np.array([[1, 0], [0, -1]])
        eigs = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]])
        expected = np.array([[0, 1], [1, 0]])
        np.allclose(global_to_local(matrix, eigs), expected)

    def test_local_to_global(self):
        matrix = np.array([[0, 1], [1, 0]])
        eigs = 1/np.sqrt(2) * np.array([[1, 1], [1, -1]])
        expected = np.array([[1, 0], [0, -1]])
        np.allclose(local_to_global(matrix, eigs), expected)

if __name__ == '__main__':
    unittest.main()
