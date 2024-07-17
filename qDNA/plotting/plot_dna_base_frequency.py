"""
This module provides functions to calculate the frequency of DNA bases in a given set of DNA sequences 
and to plot the frequency of these bases against the exciton lifetime.
"""

import numpy as np
import matplotlib.pyplot as plt


def dna_base_counter(dna_base, dna_dict):
    """
    Calculates the frequency of a specific DNA base in a given dictionary of DNA sequences.

    Parameters
    ----------
    dna_base : str
        The DNA base to count (e.g., 'A', 'T', 'G', 'C').
    dna_dict : Dict[str, float]
        Dictionary where keys are DNA sequences and values are associated float values (e.g., lifetimes).

    Returns
    -------
    np.ndarray
        Array of base frequencies.
    """
    base_frequency = []
    counter = 0
    for sequence_num, (dna_sequence, _) in enumerate(dna_dict.items()):
        num_sites_per_strand = len(dna_sequence)
        counter += dna_sequence.count(dna_base)
        base_frequency.append(counter / (num_sites_per_strand * (sequence_num + 1)))
    return np.array(base_frequency)


def plot_dna_base_frequency(lifetime_dict, cutoff_num=10):
    """
    Plots the frequency of A-T and G-C base pairs against the exciton lifetime.

    Parameters
    ----------
    lifetime_dict : Dict[str, float]
        Dictionary where keys are DNA sequences and values are lifetimes in femtoseconds.
    cutoff_num : int, optional
        The number of initial sequences to exclude from the plot, by default 10.
    """
    fig, ax = plt.subplots(1, 1, figsize=(25, 5))

    lifetimes = (
        np.array(list(lifetime_dict.values())[cutoff_num:]) * 1e-3
    )  # Convert fs to ps
    freq_a = dna_base_counter("A", lifetime_dict)[cutoff_num:]
    freq_t = dna_base_counter("T", lifetime_dict)[cutoff_num:]
    freq_g = dna_base_counter("G", lifetime_dict)[cutoff_num:]
    freq_c = dna_base_counter("C", lifetime_dict)[cutoff_num:]

    ax.plot(lifetimes, (freq_a + freq_t) * 100, "o--")
    ax.plot(lifetimes, (freq_g + freq_c) * 100, "o--")

    ax.set_ylabel("Percentage", fontsize=25)
    ax.set_xlabel("Exciton lifetime (ps)", fontsize=25)
    ax.invert_xaxis()
    ax.legend(["A-T", "G-C"], prop={"size": 30}, loc="upper right")
    ax.axhline(50, linestyle="--", color="black", alpha=0.8, lw=2.2)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    return fig, ax
