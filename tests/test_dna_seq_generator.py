import pytest
from DNA import DNA_Seq, create_upper_strands

# Run test in a Jupyter Notebook
# !python -m pytest -vv tests/test_dna_seq_generator.py --disable-pytest-warnings

@pytest.mark.parametrize("seq, mode, methylated, expected", [
    ('GcG', 'ELM', True, ('GcG', 'CGc')),
    ('GcG', 'WM', True, ('GcG',)),
    ('GcG', 'FLM', True, ('BBB', 'GcG', 'CGc', 'BBB'))
])
def test_DNA_Seq(seq, mode, methylated, expected):
    assert DNA_Seq(seq, mode, methylated).dna_seq == expected

@pytest.mark.parametrize("length, bases, expected", [
    (3, ['A', 'T'], ['AAA', 'AAT', 'ATA', 'ATT', 'TAA', 'TAT', 'TTA', 'TTT'])
])
def test_create_upper_strands(length, bases, expected):
    assert create_upper_strands(length, bases) == expected

if __name__ == "__main__":
    pytest.main()