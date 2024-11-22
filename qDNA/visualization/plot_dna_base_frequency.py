"""
This module provides functions to calculate the frequency of DNA bases in a given set of DNA sequences
and to plot the frequency of these bases against the exciton lifetime.
"""

import numpy as np
import matplotlib.pyplot as plt

from . import COLORS_DNA_BASES

# ------------------------------------------------------------------------------


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
    base_counter = 0
    for sequence_num, (dna_sequence, _) in enumerate(dna_dict.items()):
        num_sites_per_strand = len(dna_sequence)
        num_bases = num_sites_per_strand * (sequence_num + 1)
        # count how many times the base appears in the sequence
        base_counter += dna_sequence.count(dna_base)
        base_frequency.append(base_counter / num_bases)
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

    # calculation
    lifetimes = np.array(list(lifetime_dict.values())[cutoff_num:])  # in fs
    lifetimes *= 1e-3  # Convert fs to ps

    base_freq_dict = {}
    for base in ["A", "T", "G", "C"]:
        base_freq_dict[base] = dna_base_counter(base, lifetime_dict)[cutoff_num:]

    # plotting
    fig, ax = plt.subplots(1, 1, figsize=(12, 4))
    ax.plot(
        lifetimes,
        (base_freq_dict["G"] + base_freq_dict["C"]) * 100,
        "o--",
        markersize=4,
        color=COLORS_DNA_BASES["C"],
    )
    ax.plot(
        lifetimes,
        (base_freq_dict["A"] + base_freq_dict["T"]) * 100,
        "o--",
        markersize=4,
        color=COLORS_DNA_BASES["T"],
    )

    # plot settings
    ax.set_ylabel("Percentage")
    ax.set_xlabel("Exciton lifetime [ps]")
    ax.invert_xaxis()
    ax.legend(["G-C", "A-T"], loc="upper right")
    ax.axhline(50, linestyle="--", color="black", alpha=0.8, lw=2.2)
    plt.xticks()
    plt.yticks()

    return fig, ax
