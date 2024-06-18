from utils import sorted_dictionary, my_load
import numpy as np

def get_dominant_dict(dominant_filepath):
    ''' returns soreted dictionary from a given filepath '''
    dominant_dict = sorted_dictionary( my_load(dominant_filepath) )
    return dominant_dict

def get_sorted_dict(dominant_filepath, filepath):
    ''' sorts the dict according to the ordering given by the dict in the dominant filepath '''
    dominant_dict = get_dominant_dict(dominant_filepath)
    if dominant_filepath == filepath: return dominant_dict
    sequence_ordering = list(dominant_dict.keys())
    lifetime_dict = my_load( filepath )
    sorted_keys = sorted(lifetime_dict.keys(), key=lambda x: sequence_ordering.index(x))
    lifetime_dict = {key: lifetime_dict[key] for key in sorted_keys}
    return lifetime_dict

def get_correlation(dominant_filepath, filepath):
    ''' returns the correlation coefficient between dicts saved in the given filpaths when ordered
    according to the dict saved in the first path '''
    A = get_dominant_dict(dominant_filepath)
    B = get_sorted_dict(dominant_filepath, filepath)
    correlation = np.corrcoef( list(B.values()),list(A.values()) )[0, 1]
    return correlation