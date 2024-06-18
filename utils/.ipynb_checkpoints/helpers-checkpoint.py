import time 
from functools import wraps 

__all_ = ['timeit', 'sorted_dictionary']

# ------------------------------------ wrappers ----------------------------------------

def timeit(f):
    # this decorator measures the execution time of a function
    # use it by adding '@timeit' above the definition of the function
    @wraps(f)
    def timed(*args, **kw):
        start = time.time()
        result = f(*args, **kw)
        end = time.time()
        print( f'Time: {round(end-start,2)} s' )
        return result
    return timed

# ---------------------------------------------------------------------------------------------------------

def sorted_dictionary(dictionary, reverse = True):
    # returns the dictionary sorted by value 
    dictionary = dict(sorted(dictionary.items(), key = lambda item: item[1], reverse=reverse))
    return dictionary