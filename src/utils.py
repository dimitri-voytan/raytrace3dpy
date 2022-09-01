from time import time
from functools import wraps


def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        print('func:%r took: %2.4f sec' %
              (f.__name__, te - ts))
        return result
    return wrap


def check_complex(init_conds):
    import numpy as np
    np.csingle
    np.cdouble
    np.clongdouble
