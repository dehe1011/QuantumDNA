from typing import List, Tuple, Type
from itertools import product
import re

from utils import get_config

# Shortcuts:
# seq: sequence
# props: properties
# A: adenine
# C: cytosine
# G: guanine
# T: thymine
# c: 5'-methylcytosine

__all__ = ['DNA_Seq', 'DNASeqType', 'create_upper_strands', 'TB_MODELS_PROPS']

TB_MODELS_PROPS = {
    'WM': {
        'backbone': False,
        'double_stranded': False,
        'num_strands': 1,
        'diagonal_hopping': False,
        },
    'LM': { 
        'backbone': False,
        'double_stranded': True,
        'num_strands': 2,
        'diagonal_hopping': False,
        },
    'ELM': { 
        'backbone': False,
        'double_stranded': True,
        'num_strands': 2,
        'diagonal_hopping': True,
        },
    'FWM': {
        'backbone': True,
        'double_stranded': False,
        'num_strands': 3,
        'diagonal_hopping': False,
        },
    'FLM': { 
        'backbone': True,
        'double_stranded': True,
        'num_strands': 4,
        'diagonal_hopping': False,
        },
    'FELM': {
        'backbone': True,
        'double_stranded': True,
        'num_strands': 4,
        'diagonal_hopping': True,
        },
    'FC': {
        'backbone': True,
        'double_stranded': True,
        'num_strands': 4,
        'diagonal_hopping': True,
        },
}

# -------------------------------------------------------------------------

class DNA_Seq:
    def __init__(self, upper_strand: str, tb_model_name: str, methylated: bool = True):
        self.upper_strand = upper_strand
        self.methylated = methylated
        self.tb_model_name = tb_model_name
        
        self.tb_model_props = TB_MODELS_PROPS[self.tb_model_name]
        self.backbone = self.tb_model_props['backbone']
        self.double_stranded = self.tb_model_props['double_stranded']
        self.num_strands = self.tb_model_props['num_strands']
        
        self.num_sites_per_strand = len(self.upper_strand)
        self.tb_dims = (self.num_strands, self.num_sites_per_strand)
        
        self.complementary_base_dict = {"A": "T", "T": "A", "G": "C", "C": "G", "c": "G"}
        self.dna_seq = self._create_dna_seq()

    def __vars__(self) -> dict:
        return vars(self)

    def __repr__(self) -> str:
        return f"DNA_Seq({self.upper_strand}, {self.tb_model_name}, methylated = {self.methylated})"
        
    def __eq__(self, other):
        return self.__repr__() == other.__repr__()

    # ------------------------------------------------------------------
    
    def _create_dna_seq(self) -> Tuple[str]:
        if self.double_stranded:
            self.lower_strand = ''.join(list(self.complementary_base_dict[dna_base] for dna_base in self.upper_strand))
            if self.methylated:
                self._add_methylation()
        if self.backbone: 
            self.backbone_strand = "B" * len(self.upper_strand)
        if self.double_stranded and self.backbone:
            return (self.backbone_strand, self.upper_strand, self.lower_strand, self.backbone_strand)
        if self.double_stranded and not self.backbone:
            return (self.upper_strand, self.lower_strand)
        if not self.double_stranded and self.backbone:
            return (self.backbone_strand, self.upper_strand, self.backbone_strand)
        if not self.double_stranded and not self.backbone:
            return (self.upper_strand, )
            
    def _add_methylation(self):
        matches = [match.start() for match in re.finditer('cG', self.upper_strand)]
        lower_strand_list = list(self.lower_strand)
        for match in matches:
            lower_strand_list[match + 1] = 'c'
        self.lower_strand = ''.join(lower_strand_list)

DNASeqType = Type[DNA_Seq]

# -----------------------------------------------------------------------------

def create_upper_strands(num_dna_bases: int, dna_bases: List[str]) -> List[str]:
    """
    Creates all possible permutations of the given 
    """
    DNA_BASES = get_config()['DNA_BASES']
    assert all([dna_base in DNA_BASES for dna_base in dna_bases]), f"elements of dna_bases must be in {DNA_BASES}"
    return list(''.join(upper_strand) for upper_strand in product(dna_bases, repeat=num_dna_bases))
