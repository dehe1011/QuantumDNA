import pytest
import numpy as np

from qDNA.evaluation import calc_dipole, calc_dipole_moment


@pytest.mark.parametrize(
    "upper_strand, tb_model_name, average, expected",
    [("GCG", "ELM", True, 2.951734389657976)],
)
def test_calc_dipole(upper_strand, tb_model_name, average, expected):
    assert np.allclose(calc_dipole(upper_strand, tb_model_name, average), expected)


@pytest.mark.parametrize(
    "upper_strand, tb_model_name, expected", [("GCG", "ELM", 14.177784530660903)]
)
def test_calc_dipole_moment(upper_strand, tb_model_name, expected):
    assert np.allclose(calc_dipole_moment(upper_strand, tb_model_name), expected)
