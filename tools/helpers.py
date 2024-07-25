import time
from functools import wraps
import os
import numpy as np

from .save_load import my_load

__all_ = [
    "timeit",
    "sorted_dict",
    "get_dominant_dict",
    "get_sorted_dict",
    "get_correlation",
]

# ------------------------------------ wrappers ----------------------------------------


def timeit(f):
    """
    This decorator measures the execution time of a function. Use it by adding '@timeit' above the definition of the function.
    """

    @wraps(f)
    def timed(*args, **kw):
        start = time.time()
        result = f(*args, **kw)
        end = time.time()
        print(f"Time: {round(end-start,2)} s")
        return result

    return timed


# ---------------------------------------------------------------------------------------------------------


def sorted_dict(dictionary, reverse=True):
    """
    Sorts the dictionary by value.
    """
    return dict(sorted(dictionary.items(), key=lambda item: item[1], reverse=reverse))


def get_dominant_dict(dominant_filename, directory):
    """returns soreted dictionary from a given filepath"""
    dominant_dict = sorted_dict(my_load(dominant_filename, directory=directory))
    return dominant_dict


def get_sorted_dict(dominant_filename, filename, directory):
    """sorts the dict according to the ordering given by the dict in the dominant filepath"""
    dominant_dict = get_dominant_dict(dominant_filename, directory)
    if dominant_filename == filename:
        return dominant_dict
    sequence_ordering = list(dominant_dict.keys())
    lifetime_dict = my_load(filename, directory=directory)
    sorted_keys = sorted(lifetime_dict.keys(), key=lambda x: sequence_ordering.index(x))
    lifetime_dict = {key: lifetime_dict[key] for key in sorted_keys}
    return lifetime_dict


def get_correlation(dominant_filename, filename, directory):
    """returns the correlation coefficient between dicts saved in the given filpaths when ordered
    according to the dict saved in the first path"""
    A = get_dominant_dict(dominant_filename, directory)
    B = get_sorted_dict(dominant_filename, filename, directory)
    correlation = np.corrcoef(list(B.values()), list(A.values()))[0, 1]
    return correlation
