"""
Module for generating DNA sequences. 
"""

from itertools import product
import re

from qDNA.tools import get_config

# Shortcuts:
# seq: sequence
# props: properties
# A: adenine
# C: cytosine
# G: guanine
# T: thymine
# c: 5'-methylcytosine

__all__ = ["DNA_Seq", "create_upper_strands", "TB_MODELS_PROPS"]

TB_MODELS_PROPS = {
    "WM": {
        "backbone": False,
        "double_stranded": False,
        "num_strands": 1,
        "diagonal_hopping": False,
    },
    "LM": {
        "backbone": False,
        "double_stranded": True,
        "num_strands": 2,
        "diagonal_hopping": False,
    },
    "ELM": {
        "backbone": False,
        "double_stranded": True,
        "num_strands": 2,
        "diagonal_hopping": True,
    },
    "FWM": {
        "backbone": True,
        "double_stranded": False,
        "num_strands": 3,
        "diagonal_hopping": False,
    },
    "FLM": {
        "backbone": True,
        "double_stranded": True,
        "num_strands": 4,
        "diagonal_hopping": False,
    },
    "FELM": {
        "backbone": True,
        "double_stranded": True,
        "num_strands": 4,
        "diagonal_hopping": True,
    },
    "FC": {
        "backbone": True,
        "double_stranded": True,
        "num_strands": 4,
        "diagonal_hopping": True,
    },
}


class DNA_Seq:
    """
    A class representing a DNA sequence with optional methylation and backbone properties.

    Parameters
    ----------
    upper_strand : str
        The sequence of the upper strand of the DNA.
    tb_model_name : str
        The name of the tight-binding model.
    methylated : bool, optional
        Whether the DNA sequence is methylated (default is True).

    Attributes
    ----------
    upper_strand : str
        The sequence of the upper strand of the DNA.
    lower_strand : str
        The sequence of the lower strand of the DNA.
    methylated : bool
        Indicates if the DNA is methylated.
    tb_model_name : str
        The name of the tight-binding model.
    tb_model_props : dict
        Properties of the tight-binding model.
    backbone : bool
        Indicates if the model includes a backbone.
    double_stranded : bool
        Indicates if the DNA is double-stranded.
    num_strands : int
        Number of strands in the model.
    num_sites_per_strand : int
        Number of sites per strand.
    tb_dims : Tuple[int, int]
        Dimensions of the tight-binding model.
    complementary_base_dict : dict
        Mapping of bases to their complements.
    dna_seq : Tuple[str]
        The generated DNA sequence including backbone if applicable.
    """

    def __init__(self, upper_strand, tb_model_name, methylated=True):
        self.upper_strand = upper_strand
        self.lower_strand = ""
        self.methylated = methylated
        self.tb_model_name = tb_model_name

        self.tb_model_props = TB_MODELS_PROPS[self.tb_model_name]
        self.backbone = self.tb_model_props["backbone"]
        self.double_stranded = self.tb_model_props["double_stranded"]
        self.num_strands = self.tb_model_props["num_strands"]

        self.num_sites_per_strand = len(self.upper_strand)
        self.tb_dims = (self.num_strands, self.num_sites_per_strand)

        self.complementary_base_dict = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "c": "G",
        }
        self.dna_seq = self._create_dna_seq()

    def __vars__(self) -> dict:
        """
        Returns the instance variables as a dictionary.
        """
        return vars(self)

    def __repr__(self) -> str:
        """
        Returns a string representation of the DNA_Seq instance.
        """
        return f"DNA_Seq({self.upper_strand}, {self.tb_model_name}, methylated={self.methylated})"

    def __eq__(self, other) -> bool:
        """
        Compares two DNA_Seq instances for equality.
        """
        return self.__repr__() == other.__repr__()

    def _create_dna_seq(self):
        """
        Creates the DNA sequence based on the properties of the model.
        """
        if self.double_stranded:
            self.lower_strand = "".join(
                self.complementary_base_dict[dna_base] for dna_base in self.upper_strand
            )
            if self.methylated:
                self._add_methylation()
        if self.backbone:
            self.backbone_strand = "B" * len(self.upper_strand)
        if self.double_stranded and self.backbone:
            return (
                self.backbone_strand,
                self.upper_strand,
                self.lower_strand,
                self.backbone_strand,
            )
        if self.double_stranded and not self.backbone:
            return (self.upper_strand, self.lower_strand)
        if not self.double_stranded and self.backbone:
            return (self.backbone_strand, self.upper_strand, self.backbone_strand)
        return (self.upper_strand,)

    def _add_methylation(self):
        """
        Adds methylation to the DNA sequence where 'cG' is found in the upper strand.
        """
        matches = [match.start() for match in re.finditer("cG", self.upper_strand)]
        lower_strand_list = list(self.lower_strand)
        for match in matches:
            lower_strand_list[match + 1] = "c"
        self.lower_strand = "".join(lower_strand_list)


def create_upper_strands(num_dna_bases, dna_bases):
    """
    Creates all possible permutations of the given DNA bases for the specified length.

    Parameters
    ----------
    num_dna_bases : int
        The number of bases in the DNA strand.
    dna_bases : List[str]
        The list of DNA bases to use for permutations.

    Returns
    -------
    List[str]
        A list of all possible permutations of the given DNA bases.
    """
    DNA_BASES = get_config()["DNA_BASES"]
    assert all(
        [dna_base in DNA_BASES for dna_base in dna_bases]
    ), f"Elements of dna_bases must be in {DNA_BASES}"
    return list(
        "".join(upper_strand)
        for upper_strand in product(dna_bases, repeat=num_dna_bases)
    )
