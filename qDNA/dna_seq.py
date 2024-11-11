"""
This module provides functionality for representing and manipulating DNA sequences with optional methylation and backbone properties.

Shortcuts:
----------

- seq: sequence
- props: properties
- A: adenine
- C: cytosine
- G: guanine
- T: thymine
- c: 5'-methylcytosine
"""

import re
from itertools import product

from qDNA.tools import DNA_BASES, TB_MODELS_PROPS

__all__ = ["DNA_Seq", "create_upper_strands"]

# ------------------------------------------------


class DNA_Seq:
    """
    A class to represent a DNA sequence and its properties based on a tight-binding model.

    Parameters
    ----------
    upper_strand : str
        The upper strand of the DNA sequence.
    tb_model_name : str
        The name of the tight-binding model.
    methylated : bool, optional
        Indicates whether the DNA sequence is methylated (default is True).
    lower_strand : str, optional
        The lower strand of the DNA sequence (default is None).

    Attributes
    ----------
    upper_strand : str
        The upper strand of the DNA sequence.
    lower_strand : str
        The lower strand of the DNA sequence.
    methylated : bool
        Indicates whether the DNA sequence is methylated.
    tb_model_name : str
        The name of the tight-binding model.
    tb_model_props : dict
        Properties of the tight-binding model.
    backbone : bool
        Indicates whether the model includes a backbone.
    double_stranded : bool
        Indicates whether the DNA is double-stranded.
    num_strands : int
        Number of strands in the DNA sequence.
    num_sites_per_strand : int
        Number of sites per strand in the DNA sequence.
    tb_dims : tuple
        Dimensions of the tight-binding model.
    complementary_base_dict : dict
        Dictionary mapping each base to its complementary base.
    dna_seq : tuple
        The generated DNA sequence.
    """

    def __init__(self, upper_strand, tb_model_name, methylated=True, lower_strand=None):
        # Initialize the DNA sequence
        self.upper_strand = upper_strand
        self.lower_strand = lower_strand
        self.methylated = methylated
        self.tb_model_name = tb_model_name

        # Get the properties of the tight-binding model
        self.tb_model_props = TB_MODELS_PROPS[self.tb_model_name]
        self.backbone = self.tb_model_props["backbone"]
        self.double_stranded = self.tb_model_props["double_stranded"]
        self.num_strands = self.tb_model_props["num_strands"]

        # Set dimensions of the tight-binding model
        self.num_sites_per_strand = len(self.upper_strand)
        self.tb_dims = (self.num_strands, self.num_sites_per_strand)

        self.complementary_base_dict = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "F": "G",
        }
        # Generate the DNA sequence
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

    # ------------------------------------------------------------------------

    def _create_dna_seq(self):
        """
        Create the DNA sequence based on the object's attributes.
        This method generates the DNA sequence considering whether it is double-stranded,
        methylated, and/or has a backbone. It uses the `upper_strand` attribute as the
        primary sequence and generates the `lower_strand` if the DNA is double-stranded.

        Returns
        -------
        tuple
            A tuple representing the DNA sequence. The contents of the tuple vary based
            on the attributes:
              If double-stranded and has a backbone: (backbone_strand, upper_strand, lower_strand, backbone_strand)
              If double-stranded and no backbone: (upper_strand, lower_strand)
              If single-stranded and has a backbone: (backbone_strand, upper_strand, backbone_strand)
              If single-stranded and no backbone: (upper_strand,)
        """

        # Generate the lower strand if it is not provided and the model is double-stranded
        if self.double_stranded:
            if not self.lower_strand:
                self.lower_strand = "".join(
                    self.complementary_base_dict[dna_base]
                    for dna_base in self.upper_strand
                )
                if self.methylated:
                    self._add_methylation()

        # Generate the backbone strand if the model includes a backbone
        if self.backbone:
            self.backbone_strand = "B" * len(self.upper_strand)

        # Return the DNA sequence based on the model properties
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
        Adds methylation to the lower DNA strand according to the fragile X syndrome.
        This method searches for occurrences of the sequence "FG" in the upper DNA strand.
        For each occurrence, it modifies the corresponding position in the lower DNA strand
        by changing the character following the match to "F".
        """
        # Find all occurrences of "cG" in the upper DNA strand
        matches = [match.start() for match in re.finditer("FG", self.upper_strand)]

        # Modify the lower DNA strand based on the matches
        lower_strand_list = list(self.lower_strand)
        for match in matches:
            lower_strand_list[match + 1] = "F"
        self.lower_strand = "".join(lower_strand_list)


def create_upper_strands(num_dna_bases, dna_bases):
    """
    Generate all possible upper DNA strands of a given length using specified DNA bases.

    Parameters
    ----------
    num_dna_bases : int
        The number of DNA bases in each strand.
    dna_bases : list of str
        A list of DNA bases to use for generating the strands.

    Returns
    -------
    list of str
        A list containing all possible upper DNA strands of the specified length.

    Raises
    ------
    AssertionError
        If any element of `dna_bases` is not in the configured DNA bases.
    """

    # Check that the DNA bases are valid
    assert all(
        [dna_base in DNA_BASES for dna_base in dna_bases]
    ), f"Elements of dna_bases must be in {DNA_BASES}"

    # Generate all possible upper DNA strands
    return list(
        "".join(upper_strand)
        for upper_strand in product(dna_bases, repeat=num_dna_bases)
    )
