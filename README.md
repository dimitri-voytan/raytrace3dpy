# raytrace3dpy

This code is designed to solve one-point ray tracing problems in 3D heterogenous media.
The code solves the ODE ray system

$$
\frac{d\mathbf{x}}{d\lambda} = \frac{1}{S(\mathbf{x})}\mathbf{p}
$$

$$
\frac{d\mathbf{p}}{d\lambda} = \nabla S(\mathbf{x})
$$

$$
\frac{d T}{d\lambda} = S(\mathbf{x})
$$

Where $S(\mathbf{x}) = \frac{1}{V(\mathbf{x})}$ is the slowness or reciprocal of velocity, $\mathbf{x}$ is the spatial coordinate of the ray, $\mathbf{p}$ is the ray vecotor, $T$ is the traveltime along the ray, and $\lambda$ is a paramater that increases monotonically along the ray; it has the physical meaning of length along the ray

The user provides a velocity model, an initial source location $(x_0, y_0, z_0)$, and a takeoff direction specified by the inclination angle, $\alpha$ and azimuth angle, $\beta$.  The system is integrated by `scipy.integrate.solve_ivp` which by default uses an explicit Runge-Kutta method of order 5(4). For running the solver, the stopping "time" (alluding to the dependent variable of ODEs often being time) of the integration $\lambda_f$ can be set arbitrarily large, because the solver stops when a ray exits the domain.

## Takeoff angle

The takeoff direction is specified with a pair of angles $(\alpha, \beta)$. $\alpha$ is the inclination angle and $\beta$ is the azimuth. $\alpha=0$ when $\mathbf{p}$ is aligned with the $z$ direction (downward) and increases counterclockwise. $\beta=0$ when $\mathbf{p}$ is aligned with the $x$ axis and increases clockwise.

## Tracing multiple rays

Multiple rays can be traced by passing a list of tuples for the source coordinates and takeoff angles e.g. `srcs=[(0,0,0), (0,0,0)]` `angles=[(0,45), (0,30)]`. Tracing can be run in parallel by passing `parallel=True` and `n_procs` to the runner. For example

```
out = tracer.run(parallel=True, n_procs=6)
```

# Requirements

The only dependancies are Scipy and Numpy. An installation of Scipy includes numpy so a working environemt can be built by

```
conda install -c conda-forge scipy 
```

# Getting Started

The tutorial notebook `tutorial.ipynb` shows a demo example on a constant gradient of velocity model.
