import pytest

from qDNA import DNA_Seq, create_upper_strands


@pytest.mark.parametrize(
    "seq, mode, methylated, lower_strand, expected",
    [
        ("GCG", "ELM", True, None, ("GCG", "CGC")),
        ("GCG", "WM", True, None, ("GCG",)),
        ("GCG", "FLM", True, None, ("BBB", "GCG", "CGC", "BBB")),
        ("GCGCG", "LM", False, "cGcGC", ("GCGCG", "cGcGC")),
    ],
)
def test_DNA_Seq(seq, mode, methylated, lower_strand, expected):
    assert (
        DNA_Seq(seq, mode, methylated=methylated, lower_strand=lower_strand).dna_seq
        == expected
    )


@pytest.mark.parametrize(
    "length, bases, expected",
    [(3, ["A", "T"], ["AAA", "AAT", "ATA", "ATT", "TAA", "TAT", "TTA", "TTT"])],
)
def test_create_upper_strands(length, bases, expected):
    assert create_upper_strands(length, bases) == expected
