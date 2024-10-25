import pytest
import numpy as np

from qDNA.dynamics import get_me_solver
from qDNA.evaluation import get_therm_eq_state, get_deph_eq_state


@pytest.mark.parametrize(
    "upper_strand, tb_model_name, expected",
    [
        (
            "GC",
            "WM",
            np.array(
                [
                    [0.25, -0.00966561, 0.0921533, -0.00356287],
                    [-0.00966561, 0.25, -0.00356287, 0.0921533],
                    [0.0921533, -0.00356287, 0.25, -0.00966561],
                    [-0.00356287, 0.0921533, -0.00966561, 0.25],
                ]
            ),
        )
    ],
)
def test_get_therm_eq_state(upper_strand, tb_model_name, expected):
    me_solver = get_me_solver(upper_strand, tb_model_name)
    assert np.allclose(get_therm_eq_state(me_solver), expected)


@pytest.mark.parametrize(
    "upper_strand, tb_model_name, loc_deph_rate, expected",
    [
        (
            "GC",
            "WM",
            True,
            np.array(
                [
                    [0.2, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.2, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.2, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.2, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.2],
                ]
            ),
        )
    ],
)
def test_get_deph_eq_state_loc(upper_strand, tb_model_name, loc_deph_rate, expected):
    me_solver1 = get_me_solver(upper_strand, tb_model_name, loc_deph_rate=loc_deph_rate)
    assert np.allclose(get_deph_eq_state(me_solver1), expected)


@pytest.mark.parametrize(
    "upper_strand, tb_model_name, glob_deph_rate, expected",
    [
        (
            "GC",
            "WM",
            True,
            np.array(
                [
                    [
                        2.50000000e-01 + 0.0j,
                        -1.36383960e-13 + 0.0j,
                        -4.97032970e-14 + 0.0j,
                        -7.98319744e-14 + 0.0j,
                    ],
                    [
                        -1.36383960e-13 + 0.0j,
                        2.50000000e-01 + 0.0j,
                        -7.99152411e-14 + 0.0j,
                        -4.97796249e-14 + 0.0j,
                    ],
                    [
                        -4.97032970e-14 + 0.0j,
                        -7.99152411e-14 + 0.0j,
                        2.50000000e-01 + 0.0j,
                        -1.36335387e-13 + 0.0j,
                    ],
                    [
                        -7.98250355e-14 + 0.0j,
                        -4.97865638e-14 + 0.0j,
                        -1.36349265e-13 + 0.0j,
                        2.50000000e-01 + 0.0j,
                    ],
                ]
            ),
        )
    ],
)
def test_get_deph_eq_state_glob(upper_strand, tb_model_name, glob_deph_rate, expected):
    me_solver2 = get_me_solver(
        upper_strand,
        tb_model_name,
        glob_deph_rate=glob_deph_rate,
        loc_deph_rate=False,
        relaxation=False,
    )
    assert np.allclose(get_deph_eq_state(me_solver2), expected)
