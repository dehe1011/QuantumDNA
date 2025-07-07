import os
import numpy as np
from ..io import load_json

# ----------------------------------------------------------------------


def sorted_dict(dictionary, reverse=True):
    """Sorts the dictionary by its values.

    Parameters
    ----------
    dictionary : dict
        The dictionary to be sorted.
    reverse : bool, optional
        If True, the dictionary is sorted in descending order. If False, in ascending order. Default is True.

    Returns
    -------
    dict
        A new dictionary sorted by its values.
    """
    return dict(sorted(dictionary.items(), key=lambda item: item[1], reverse=reverse))


def get_dominant_dict(dominant_filename, directory):
    """Load and sort a JSON file into a dictionary.

    Parameters
    ----------
    dominant_filename : str
        The name of the JSON file to be loaded.
    directory : str
        The directory where the JSON file is located.
    Returns
    -------
    dict
        A sorted dictionary containing the data from the JSON file.
    """

    filepath = os.path.join(directory, dominant_filename)
    dominant_dict = sorted_dict(load_json(filepath))
    return dominant_dict


def get_sorted_dict(dominant_filename, filename, directory):
    """Sorts a dictionary based on the order of keys in a dominant dictionary.

    Parameters
    ----------
    dominant_filename : str
        The filename of the dominant dictionary JSON file.
    filename : str
        The filename of the dictionary JSON file to be sorted.
    directory : str
        The directory where the JSON files are located.
    Returns
    -------
    dict
        A dictionary sorted based on the order of keys in the dominant dictionary.
    """

    dominant_dict = get_dominant_dict(dominant_filename, directory)

    # If the file contains the dominant dictionary, return the dominant dictionary
    if dominant_filename == filename:
        return dominant_dict

    # Sort the dictionary based on the order of keys in the dominant dictionary
    sequence_ordering = list(dominant_dict.keys())
    filepath = os.path.join(directory, filename)
    lifetime_dict = load_json(filepath)
    sorted_keys = sorted(lifetime_dict.keys(), key=sequence_ordering.index)
    lifetime_dict = {key: lifetime_dict[key] for key in sorted_keys}
    return lifetime_dict


def get_correlation(dominant_filename, filename, directory):
    """Calculate the correlation coefficient between two datasets.

    Parameters
    ----------
    dominant_filename : str
        The filename of the dominant dataset.
    filename : str
        The filename of the dataset to compare against the dominant dataset.
    directory : str
        The directory where the files are located.
    Returns
    -------
    float
        The correlation coefficient between the two datasets.
    Notes
    -----
    This function assumes that `get_dominant_dict` and `get_sorted_dict` are
    defined elsewhere and return dictionaries where the values are the data
    points to be correlated.
    """

    # Get the values from the dominant and sorted dictionaries
    dominant_vals = get_dominant_dict(dominant_filename, directory).values()
    sorted_vals = get_sorted_dict(dominant_filename, filename, directory).values()

    # Calculate the correlation coefficient
    correlation = np.corrcoef(list(dominant_vals), list(sorted_vals))[0, 1]
    return correlation


# ----------------------------------------------------------------------
