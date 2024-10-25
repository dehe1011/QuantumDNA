import pytest

from qDNA.evaluation import calc_exciton_transfer


@pytest.mark.parametrize(
    "upper_strand, tb_model_name, expected",
    [
        (
            "GC",
            "ELM",
            (
                {
                    "electron": 0.4957243676173506,
                    "hole": 0.45494093255325446,
                    "exciton": 0.2014896064994649,
                },
                {
                    "electron": 0.5042756323826466,
                    "hole": 0.5450590674467428,
                    "exciton": 0.2465940907884935,
                },
            ),
        )
    ],
)
def test_calc_exciton_transfer(upper_strand, tb_model_name, expected):
    assert calc_exciton_transfer(upper_strand, tb_model_name) == expected
