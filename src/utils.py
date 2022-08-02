def terminalevent(f):
    f.terminal = True
    return f

class Stopper:
    '''
    Helper class for stopping when the ray hits the boundary
    '''
    def __init__(self,
                 idx : int,
                 bound : float):
        '''
        Bound is the boundary value.
        '''
        self.idx = idx
        self.bound = bound

    @terminalevent
    def stopping_criteria(self, t, y_bar):
        '''
        Stop if a ray exits the domain (i.e. the return of this function changes sign)
        '''
        return y_bar[self.idx]-self.bound

def check_complex(init_conds):
    import numpy as np
    np.csingle
    np.cdouble
    np.clongdouble