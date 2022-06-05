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

given a takeoff direction specified by the inclination angle, $\alpha$, azimuth angle, $\beta$, and a velocity model. The system is integrated by `scipy.integrate.solve_ivp` which by default uses an explicit Runge-Kutta method of order 5(4). Note that the parameter $\lambda$ has a physical meaning of length along the ray; any other parameter that increases monotonically along the ray could be used. For running the solver, the stopping "time" of the integration $\lambda_f$ can be set arbitrarily large.

## Takeoff angle

The takeoff direction is specified with a pair of angles $(\alpha, \beta)$. $\alpha$ is the inclination angle and $\beta$ is the azimuth. $\alpha=0$ when $\mathbf{p}$ is aligned with the $z$ direction (downward) and increases counterclockwise. $\beta=0$ when $\mathbf{p}$ is aligned with the $x$ axis and increases clockwise.

## Tracing multiple rays

Multiple rays can be traced by passing a list of tuples for the source coordinates and takeoff angles e.g. `srcs=[(0,0,0), (0,0,0)]` `angles=[(0,45), (0,30)]`. Tracing can be run in parallel by passing `parallel=True` and `n_procs` to the runner. For example

```
out = tracer.run(parallel=True, n_procs=6)
```

# Requirements

The only dependancies are scipy and numpy.

# Getting Started

The tutorial notebook `tutorial.ipynb` shows a demo example on a constant gradient of velocity model.
