from scipy.integrate import solve_ivp
import numpy as np
from typing import Tuple
from domain import Domain, Boundaries, Slowness
from utils import timing


class OnePointTrace3D():

    def __init__(self,
                 src_coords: Tuple,
                 takeoff_angles: Tuple,
                 domain: Domain,
                 slowness: Slowness,
                 lambda_final: float
                 ):

        self.lf = lambda_final
        self.domain = domain
        self.s = slowness

        self.y0 = []

        # If a single tuple is passed (i.e. one source and takeoff)
        # wrap in a list so that zip expands correctly.
        if np.ndim(src_coords) <= 1:
            src_coords = [src_coords]

        if np.ndim(takeoff_angles) <= 1:
            takeoff_angles = [takeoff_angles]

        for src, angle in zip(src_coords, takeoff_angles):
            self.y0.append(np.array([*src, *self.initialize_p(src, angle), 0]))

    def initialize_p(self, src, angle):
        alpha, beta = [OnePointTrace3D.deg2rad(item) for item in angle]
        x_idx, y_idx, z_idx = self.find_nearest_idx(*src)
        s_0 = self.s[x_idx, y_idx, z_idx]
        p_1 = s_0 * np.sin(alpha) * np.cos(beta)
        p_2 = s_0 * np.sin(alpha) * np.sin(beta)
        p_3 = s_0 * np.cos(alpha)
        return (p_1, p_2, p_3)

    def find_nearest_idx(self, x, y, z):
        idx_x = int(np.round((x - self.domain.x_min) / self.domain.dx))
        idx_y = int(np.round((y - self.domain.y_min) / self.domain.dy))
        idx_z = int(np.round((z - self.domain.z_min) / self.domain.dz))
        return idx_x, idx_y, idx_z

    @staticmethod
    def deg2rad(theta):
        return theta * (np.pi / 180)

    def rhs(self, t, y_bar):
        '''
        Right hand side of the ODE
        '''
        x, y, z, p1, p2, p3, T = y_bar

        x_idx, y_idx, z_idx = self.find_nearest_idx(x, y, z)

        try:
            s_i = self.s[x_idx, y_idx, z_idx]
            one_over_s = 1 / s_i
            ds_dx, ds_dy, ds_dz = self.s.grad[:, x_idx, y_idx, z_idx]
        except IndexError:
            # Gracefully handles domain exits if an index error occurs
            s_i = 0
            one_over_s = 0
            ds_dx, ds_dy, ds_dz = (0, 0, 0)

        return (one_over_s) * p1, (one_over_s) * p2, (one_over_s) * p3, \
            ds_dx, ds_dy, ds_dz, s_i

    def trace_single_ray(self, ray_init, events, **kwargs):
        out = solve_ivp(self.rhs,
                        (0, self.lf),
                        ray_init,
                        events=events,
                        **kwargs)
        # Store the initial conditions if they are needed later
        out['init_conds'] = ray_init
        return out

    @timing
    def run(self,
            bounds: Boundaries,
            parallel: bool = False,
            **kwargs):
        '''
        Can pass kwargs to the rk solver to change order etc.
        See scipy.integrate.solveivp for details
        '''

        out = []

        if parallel:
            raise NotImplementedError
            # It will usually be more expensive to run in parallel than to run serially
            # due to pickling of the class. 
            # from multiprocessing import Pool
            # from functools import partial

            # n_procs = kwargs.pop('n_procs', os.cpu_count())
            # with Pool(n_procs) as p:
            #     f = partial(self.trace_single_ray, events=bounds.events, **kwargs)
            #     out = p.map(f, self.y0)
        else:
            for ray_init in self.y0:
                out.append(
                    self.trace_single_ray(
                        ray_init,
                        bounds.events,
                        **kwargs))
        return out
