from itertools import product

__all__ = ['DNA_BASES', 'add_complementary_strand', 'add_backbone', 'create_DNA_sequences']

# ----------------------------------- basic DNA stuff ------------------------------------------------------

DNA_BASES = ["A", "T", "G", "C"]

def add_complementary_strand(DNAstring):
    """
    Adds the second DNA strand determined by the nucleobases on the first strand.
    
    Example:
        >>> add_complementary_strand('GAC')
        ['GAC', 'CTG']
    """
    
    if isinstance(DNAstring, list): 
        DNAstring = DNAstring[0]
    complementary_base = {"A": "T", "T": "A", "G": "C", "C": "G", "c": "g", "g": "c"}
    complementary_strand = [complementary_base[base] for base in DNAstring]
    complementary_strand = "".join(complementary_strand)
    return [DNAstring, complementary_strand]

def add_backbone(DNAstring):
    """
    Adds the sugar-phosphate backbone to a given DNA sequence.
    
    Example:
        >>> add_backbone('GAC')
        ['BBB', 'GAC', 'BBB']
    """
        
    if isinstance(DNAstring, str): 
        DNAstring = [DNAstring]
    return ["B" * len(DNAstring[0])] + DNAstring + ["B" * len(DNAstring[0])]

def create_DNA_sequences(num_bases, dna_bases=DNA_BASES, double_stranded = True, backbone = False):
    """
    Creates a list of all sequences with a given number of DNA bases/ base pairs.
    
    Examples:
        >>> create_DNA_sequences(2, dna_bases=['A','T'], double_stranded=False,backbone=True)
        [['BB', 'AA', 'BB'],
         ['BB', 'AT', 'BB'],
         ['BB', 'TA', 'BB'],
         ['BB', 'TT', 'BB']]        
    """
    
    DNA_sequences = [["".join(p)] for p in product(dna_bases, repeat=num_bases)]
    if double_stranded: 
        DNA_sequences = [add_complementary_strand(DNAstring) for DNAstring in DNA_sequences]
    if backbone: 
        DNA_sequences = [add_backbone(DNAstring) for DNAstring in DNA_sequences]
    return DNA_sequences
