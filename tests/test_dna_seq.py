import pytest

from qDNA import DNA_Seq, create_upper_strands


@pytest.mark.parametrize(
    "seq, mode, methylated, lower_strand, expected",
    [
        (list("GCG"), "ELM", True, None, (list("GCG"), list("CGC"))),
        (list("GCG"), "WM", True, None, (list("GCG"),)),
        (
            list("GCG"),
            "FLM",
            True,
            None,
            (list("BBB"), list("GCG"), list("CGC"), list("BBB")),
        ),
        (list("GCGCG"), "LM", False, list("FGFGC"), (list("GCGCG"), list("FGFGC"))),
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
