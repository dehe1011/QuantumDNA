from itertools import product

# ----------------------------------- basic DNA stuff ------------------------------------------------------

__all__ = ['DNA_BASES', 'PARTICLES', 'complementary_strand', 'DNA_sequences']

DNA_BASES = ['A','T','G','C']
PARTICLES = ['electron', 'hole', 'exciton']

def complementary_strand(DNAstring):
    # returns the second DNA strand determined by the nucleobases on the first strand 
    complementary_base = {'A':'T','T':'A','G':'C','C':'G','c':'g','g':'c'}
    complementary_strand = ''
    for i in range(len(DNAstring)):
        complementary_strand += complementary_base[DNAstring[i]]
    return complementary_strand

def DNA_sequences(N):
    # returns a list of all permutations of DNA bases of length N
    # e.g., ['AA','AT','AG','AC','TA','TT','TG','TC','GA','GT','GG','GC','CA','CT','CG','CC']
    DNA_sequences = [''.join(p) for p in product(DNA_BASES, repeat=N)]
    return DNA_sequences 
