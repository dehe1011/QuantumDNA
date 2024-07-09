import pytest
from DNA.model import get_tb_config

@pytest.mark.parametrize("input, expected", [
    (('WM', (1,2)), [
        ('E', '(0, 0)', '(0, 0)'),
        ('E', '(0, 1)', '(0, 1)'),
        ('t', '(0, 1)', '(0, 0)')
    ])
])
def test_get_tb_config_WM(input, expected):
    assert get_tb_config(*input) == expected

@pytest.mark.parametrize("input, expected", [
    (('LM', (2,2)), [
        ('E', '(0, 0)', '(0, 0)'),
        ('E', '(0, 1)', '(0, 1)'),
        ('t', '(0, 1)', '(0, 0)'),
        ('E', '(1, 0)', '(1, 0)'),
        ('E', '(1, 1)', '(1, 1)'),
        ('t', '(1, 0)', '(1, 1)'),
        ('h', '(0, 0)', '(1, 0)'),
        ('h', '(0, 1)', '(1, 1)')
    ])
])
def test_get_tb_config_LM(input, expected):
    assert get_tb_config(*input) == expected

@pytest.mark.parametrize("input, expected", [
    (('ELM', (2,2)), [
        ('E', '(0, 0)', '(0, 0)'),
        ('E', '(0, 1)', '(0, 1)'),
        ('t', '(0, 1)', '(0, 0)'),
        ('E', '(1, 0)', '(1, 0)'),
        ('E', '(1, 1)', '(1, 1)'),
        ('t', '(1, 0)', '(1, 1)'),
        ('h', '(0, 0)', '(1, 0)'),
        ('h', '(0, 1)', '(1, 1)'),
        ('r+', '(0, 0)', '(1, 1)'),
        ('r-', '(0, 1)', '(1, 0)')
    ])
])
def test_get_tb_config_ELM(input, expected):
    assert get_tb_config(*input) == expected

@pytest.mark.parametrize("input, expected", [
    (('FWM', (3,2)), [
        ('h', '(0, 0)', '(1, 0)'),
        ('h', '(0, 1)', '(1, 1)'),
        ('h', '(1, 0)', '(2, 0)'),
        ('h', '(1, 1)', '(2, 1)'),
        ('E', '(1, 0)', '(1, 0)'),
        ('E', '(1, 1)', '(1, 1)'),
        ('t', '(1, 1)', '(1, 0)')
    ])
])
def test_get_tb_config_FWM(input, expected):
    assert get_tb_config(*input) == expected

@pytest.mark.parametrize("input, expected", [
    (('FLM', (4,2)), [
        ('h', '(0, 0)', '(1, 0)'),
        ('h', '(0, 1)', '(1, 1)'),
        ('h', '(2, 0)', '(3, 0)'),
        ('h', '(2, 1)', '(3, 1)'),
        ('E', '(1, 0)', '(1, 0)'),
        ('E', '(1, 1)', '(1, 1)'),
        ('t', '(1, 1)', '(1, 0)'),
        ('E', '(2, 0)', '(2, 0)'),
        ('E', '(2, 1)', '(2, 1)'),
        ('t', '(2, 0)', '(2, 1)'),
        ('h', '(1, 0)', '(2, 0)'),
        ('h', '(1, 1)', '(2, 1)')
    ])
])
def test_get_tb_config_FLM(input, expected):
    assert get_tb_config(*input) == expected

@pytest.mark.parametrize("input, expected", [
    (('FELM', (4,2)), [
        ('h', '(0, 0)', '(1, 0)'),
        ('h', '(0, 1)', '(1, 1)'),
        ('h', '(2, 0)', '(3, 0)'),
        ('h', '(2, 1)', '(3, 1)'),
        ('E', '(1, 0)', '(1, 0)'),
        ('E', '(1, 1)', '(1, 1)'),
        ('t', '(1, 1)', '(1, 0)'),
        ('E', '(2, 0)', '(2, 0)'),
        ('E', '(2, 1)', '(2, 1)'),
        ('t', '(2, 0)', '(2, 1)'),
        ('h', '(1, 0)', '(2, 0)'),
        ('h', '(1, 1)', '(2, 1)'),
        ('r+', '(1, 0)', '(2, 1)'),
        ('r-', '(1, 1)', '(2, 0)')
    ])
])
def test_get_tb_config_FELM(input, expected):
    assert get_tb_config(*input) == expected

@pytest.mark.parametrize("input, expected", [
    (('FC', (4,2)), [
        ('h', '(0, 0)', '(1, 0)'),
        ('h', '(0, 1)', '(1, 1)'),
        ('h', '(2, 0)', '(3, 0)'),
        ('h', '(2, 1)', '(3, 1)'),
        ('E', '(1, 0)', '(1, 0)'),
        ('E', '(1, 1)', '(1, 1)'),
        ('t', '(1, 1)', '(1, 0)'),
        ('E', '(2, 0)', '(2, 0)'),
        ('E', '(2, 1)', '(2, 1)'),
        ('t', '(2, 0)', '(2, 1)'),
        ('h', '(1, 0)', '(2, 0)'),
        ('h', '(1, 1)', '(2, 1)'),
        ('r+', '(1, 0)', '(2, 1)'),
        ('r-', '(1, 1)', '(2, 0)'),
        ('t', '(0, 0)', '(0, 1)'),
        ('t', '(3, 0)', '(3, 1)')
    ])
])
def test_get_tb_config_FC(input, expected):
    assert get_tb_config(*input) == expected

if __name__ == "__main__":
    pytest.main()
