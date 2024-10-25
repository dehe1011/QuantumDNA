import pytest

from qDNA.evaluation import calc_lifetime


@pytest.mark.parametrize(
    "upper_strand, tb_model_name, kwargs, expected",
    [
        # Case where ground state population reaches threshold and unit is rad/ps
        ("GCG", "ELM", {"relax_rate": 3, "unit": "rad/ps"}, 775.5511022044088),
        # Case with no relaxation occurring in the given time
        (
            "GCG",
            "ELM",
            {"relax_rate": 0, "unit": "rad/ps"},
            "no relaxation in the given time",
        ),
    ],
)
def test_calc_lifetime(upper_strand, tb_model_name, kwargs, expected):
    assert calc_lifetime(upper_strand, tb_model_name, **kwargs) == expected
