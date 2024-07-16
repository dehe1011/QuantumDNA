from typing import Dict
import numpy as np 
import matplotlib.pyplot as plt

def dna_base_counter(dna_base: str, dna_dict: Dict[str, float]) -> np.ndarray:
    base_frequency = []
    dna_base_counter = 0 
    for sequence_num, (dna_sequence, val) in enumerate(dna_dict.items()):
        num_sites_per_strand = len(dna_sequence)
        dna_base_counter += dna_sequence.count(dna_base) 
        base_frequency.append( dna_base_counter / (num_sites_per_strand*(sequence_num+1)) )
    return np.array(base_frequency)

def plot_dna_base_frequency(lifetime_dict: Dict[str, float], cutoff_num: int = 10):
    fig,ax=plt.subplots(1,1,figsize=(25,5))
    
    # lifetimes given in fs in the dictionary, but they should be plotted in ps
    lifetimes = np.array( list(lifetime_dict.values())[10:] )*1e-3
    frequency_A = dna_base_counter('A', lifetime_dict)[cutoff_num:]
    frequency_T = dna_base_counter('T', lifetime_dict)[cutoff_num:]
    frequency_G = dna_base_counter('G', lifetime_dict)[cutoff_num:]
    frequency_C = dna_base_counter('C', lifetime_dict)[cutoff_num:]
    
    ax.plot(lifetimes, (frequency_A+frequency_T)*100, 'o--' )
    ax.plot(lifetimes, (frequency_G+frequency_C)*100, 'o--' )
    
    ax.set_ylabel('Percentage',fontsize=25)
    ax.set_xlabel('Exciton lifetime (ps)',fontsize=25) 
    ax.invert_xaxis()
    ax.legend(['A-T','G-C'],prop={'size': 30}, loc='upper right')
    ax.axhline(50,linestyle='--',color='black',alpha=0.8,lw=2.2)
    plt.xticks(fontsize=20), plt.yticks(fontsize=20)