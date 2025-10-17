import time
from functools import wraps

# ----------------------------------------------------------------------


def timeit(f):
    """
    Decorator to measure and print the execution time of a function.
    """

    @wraps(f)
    def timed(*args, **kw):
        start = time.time()
        result = f(*args, **kw)
        end = time.time()
        print(f"Time: {round(end-start,2)} s")
        return result

    return timed


# ----------------------------------------------------------------------
