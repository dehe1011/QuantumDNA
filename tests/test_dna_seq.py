import pytest
from qDNA import DNA_Seq, create_upper_strands


@pytest.mark.parametrize(
    "seq, mode, methylated, expected",
    [
        ("GCG", "ELM", True, ("GCG", "CGC")),
        ("GCG", "WM", True, ("GCG",)),
        ("GCG", "FLM", True, ("BBB", "GCG", "CGC", "BBB")),
    ],
)
def test_DNA_Seq(seq, mode, methylated, expected):
    assert DNA_Seq(seq, mode, methylated).dna_seq == expected


@pytest.mark.parametrize(
    "length, bases, expected",
    [(3, ["A", "T"], ["AAA", "AAT", "ATA", "ATT", "TAA", "TAT", "TTA", "TTT"])],
)
def test_create_upper_strands(length, bases, expected):
    assert create_upper_strands(length, bases) == expected
