# raytrace3dpy
3D one-point ray tracing code

Solves the ODE ray system

$$
\frac{d\mathbf{x}}{d\lambda} = \frac{1}{S(\mathbf{x})}\mathbf{p}
$$

$$
\frac{d\mathbf{p}}{d\lambda} = \nabla S(\mathbf{x})
$$

$$
\frac{d T}{d\lambda} = S(\mathbf{x})
$$

Given a takeoff direction specified by the inclination angle, $\alpha$, azimuth angle, $\beta$, and a velocity model solve the ODE ray system above. The system is integrated by `scipy.integrate.solve_ivp` which by default uses an explicit Runge-Kutta method of order 5(4). Note that the parameter $\lambda$ has a physical meaning of length along the ray; any other parameter that increasses monotonically along the ray could be used. It can be set arbitrarily large.


# Requirements

The only dependancies are scipy, numpy, and matplotlib.

# Getting Started

The tutorial notebook `tutorial.ipynb` shows a demo example on a constant gradient of velocity model.
