from dataclasses import dataclass
from typing import List, Tuple
import numpy as np


class Domain():
    def __init__(self,
                 x_coords: np.ndarray,
                 y_coords: np.ndarray,
                 z_coords: np.ndarray):

        self.x_coords = x_coords
        self.y_coords = y_coords
        self.z_coords = z_coords

        self.x_min, self.x_max = np.min(x_coords), np.max(x_coords)
        self.y_min, self.y_max = np.min(y_coords), np.max(y_coords)
        self.z_min, self.z_max = np.min(z_coords), np.max(z_coords)

        self.dx = x_coords[1] - x_coords[0]
        self.dy = y_coords[1] - y_coords[0]
        self.dz = z_coords[1] - z_coords[0]

        self.nx = len(x_coords)
        self.ny = len(y_coords)
        self.nz = len(z_coords)


def terminalevent(f):
    '''
    Creates an attribute of a function, f
    indicating that the function is a terminal event.

    This is used by scipy's ODE integrator
    '''
    f.terminal = True
    return f


class Stopper:
    '''
    Helper class for stopping when the ray hits the boundary.
    '''

    def __init__(self,
                 idx: int,
                 bound: float):

        self.idx = idx
        self.bound = bound

    @terminalevent
    def stopping_criteria(self, t, y_bar):
        '''
        Stop if a ray exits the domain which occurs when the return of this 
        function changes sign
        '''
        return y_bar[self.idx] - self.bound


class Boundaries:
    '''
    Domain boundaries.
    '''

    def __init__(self,
                 domain: Domain):

        self.bounds = [(0, domain.x_min + domain.dx),
                       (0, domain.x_max - domain.dx),
                       (1, domain.y_min + domain.dy),
                       (1, domain.y_max - domain.dy),
                       (2, domain.z_min + domain.dz),
                       (2, domain.z_max - domain.dz)]

        self.event_list = [Stopper(bound[0], bound[1])
                           for bound in self.bounds]
        self.events = [event.stopping_criteria for event in self.event_list]

    def add_bound(self,
                  bounds: List[Tuple[int, int]]):
        '''
        Adds a new boundary to the domain.
        '''
        new_events = [Stopper(bound[0], bound[1]) for bound in bounds]
        self.events.append(event.stopping_criteria for event in new_events)
        self.event_list.append(new_events)


@dataclass
class Slowness:
    vals: np.ndarray
    grad: np.ndarray

    def __getitem__(self, idx):
        return self.vals[idx]
