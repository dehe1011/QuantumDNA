import pytest
import numpy as np
import scipy.constants as c

from qDNA.utils.unit_converter import (
    get_conversion,
    convert_to_debye,
    get_conversion_dict,
)


# most used unit conversion
@pytest.mark.parametrize(
    "start_unit, end_unit, expected",
    [
        ("100meV", "rad/fs", (100e-3 * c.e) / c.hbar * 1e-15),
        ("rad/fs", "100meV", c.hbar / (100e-3 * c.e) * 1e15),
        ("rad/fs", "100meV", 6.582119569509065),
    ],
)
def test_get_conversion(start_unit, end_unit, expected):
    assert np.allclose(get_conversion(start_unit, end_unit), expected)


@pytest.mark.parametrize(
    "param_dict, start_unit, end_unit, expected",
    [
        (
            {"A": 2, "B": 3},
            "100meV",
            "rad/fs",
            {"A": 0.3038534895757253, "B": 0.45578023436358794},
        )
    ],
)
def test_get_conversion_dict(param_dict, start_unit, end_unit, expected):
    conversion_dict = get_conversion_dict(param_dict, start_unit, end_unit)
    assert np.allclose(list(conversion_dict.values()), list(expected.values()))


# Dipole moment for a charge separation of one DNA base (3.4 Angstrom)
@pytest.mark.parametrize("charge_separation, expected", [(3.4, 16.330896022738898)])
def test_convert_to_debye(charge_separation, expected):
    assert np.allclose(convert_to_debye(charge_separation), expected)
