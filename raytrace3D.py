from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator
import numpy as np
from typing import Tuple

'''
Solves one-point ray tracing
Requires src location, takeoff angle, and velocity
'''

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
        

class OnePointTrace3D():
    def __init__(self,
                 src_coord : Tuple,
                 takeoff_angle : Tuple,
                 velocity : np.ndarray,
                 x_coords : np.ndarray,
                 y_coords : np.ndarray,
                 z_coords : np.ndarray,
                 lf : float,
                 n_pad : int = 10):
        '''
        takeoff angle: Tuple of (inclination [0:pi), azimuth [0:2pi) )
        '''
        # Regularly sampled grid
        self.x_coords = x_coords
        self.y_coords = y_coords
        self.z_coords = z_coords

        self.x_min, self.x_max = np.min(x_coords), np.max(x_coords)
        self.y_min, self.y_max = np.min(y_coords), np.max(y_coords)
        self.z_min, self.z_max = np.min(z_coords), np.max(z_coords)

        self.dx = x_coords[1]-x_coords[0]
        self.dy = y_coords[1]-y_coords[0]
        self.dz = z_coords[1]-z_coords[0]

        self.velocity = velocity

        self.n_pad = 5

        # Create interpolator objects for common terms
        self.slowness = self.initialize_slowness()

        # Initial conditions
        self.src_coord = src_coord
        self.takeoff_angle = takeoff_angle
        self.y0 = np.array([*self.src_coord, *self.initialize_p(), 0])

        # Stopping conditions. Ray will run until it hits the boundary or tf.
        # tf has the physical meaning of length along the ray
        self.lf = lf

    def initialize_p(self):
        s_0 = self.slowness(self.src_coord)
        alpha, beta = [self.deg2rad(item) for item in self.takeoff_angle]
        p_1 = s_0*np.sin(alpha)*np.cos(beta)
        p_2 = s_0*np.sin(alpha)*np.sin(beta)
        p_3 = s_0*np.cos(alpha)
        return (p_1, p_2, p_3)

    def initialize_slowness(self):
        '''
        Add padding so that if the slowness is computed outside of the domain (when the ray escapes)
        this will be properly handled
        '''
        x_coords = np.pad(self.x_coords, self.n_pad, mode='reflect', reflect_type='odd')
        y_coords = np.pad(self.y_coords, self.n_pad, mode='reflect', reflect_type='odd')
        z_coords = np.pad(self.z_coords, self.n_pad, mode='reflect', reflect_type='odd')
        vel = np.pad(self.velocity, self.n_pad, mode='edge')

        return RegularGridInterpolator((x_coords,
                                        y_coords,
                                        z_coords),
                                        1.0/vel)

    def deg2rad(self, theta):
        return theta*(np.pi/180)

    def get_grad_s(self, s, coord):
        '''
        Simple central difference to approximate the gradient.
        '''
        x, y, z = [int(item) for item in coord]
        # Could replace with a higher order scheme?
        ds_dx = (s((x+self.dx,y,z))-s((x-self.dx,y,z)))/(2*self.dx)
        ds_dy = (s((x,y+self.dy,z))-s((x,y-self.dy,z)))/(2*self.dy)
        ds_dz = (s((x,y,z+self.dz))-s((x,y,z-self.dz)))/(2*self.dz)
        return ds_dx, ds_dy, ds_dz

    def rhs(self, t, y_bar):
        '''
        Right hand side of the ODE
        y_bar has elements x, y, z, p_1, p_2, p_3, T
        '''
        x, y, z, p1, p2, p3, T = y_bar

        s = self.slowness((x, y, z))
        one_over_s = 1/s

        ds_dx, ds_dy, ds_dz = self.get_grad_s(self.slowness, (x, y, z))

        return (one_over_s)*p1, (one_over_s)*p2, (one_over_s)*p3, ds_dx, ds_dy, ds_dz, s

    def run(self, **kwargs):
        '''
        Can pass kwargs to the rk solver
        '''
        # List of index, val which tells the solver to stop 
        # if it encounters a model boundary
        bounds = [[0, self.x_min],
                  [0, self.x_max],
                  [1, self.y_min],
                  [1, self.y_max],
                  [2, self.z_min],
                  [2, self.z_max]]

        event_list = [Stopper(bound[0], bound[1]) for bound in bounds]
        events = [event.stopping_criteria for event in event_list]

        out = solve_ivp(self.rhs, (0, self.lf), self.y0, events=events, **kwargs)
        return out
