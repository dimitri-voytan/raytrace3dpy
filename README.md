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

# Requirements

The only dependancies are Scipy and Numpy. An installation of Scipy includes Numpy so a working environemt can be built by

```
conda create -n raytrace3d python=3.x #Change x to any scipy compatible python
conda install -c conda-forge scipy 
```

# Getting Started

The tutorial notebook `tutorial.ipynb` shows an example of setting up the solver and tracing rays in a simple velocity model. Below I list a few key considerations for a user 

1.The output of a run is a list of [solution objects](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp) for each ray. It has attributes which can be accesed in a dictionary like manner to view the solution and the evaluation points. The solution is stored in the attribute `y` and the evaluation $\lambda$'s are stored in `t`. For example, if you traced a single ray and stored the solution in `out` you can access the $\mathbf{x}$ coordinates of the ray, $\mathbf{p}$ values, and traveltime $T$ by:

```
x, y, z, = out['y'][0], out['y'][1], out['y'][2]
p_0, p_1, p_2, = out['y'][3], out['y'][4], out['y'][5]
T = out['y'][-1]
lambdas = out['t']
```
The initial conditions are also available under the attribute 

```
init_conds = out['init_conds']
```
Refer to the Scipy documentation for more details about the solver and the output object.

2. The properties of the solver can be modified by passing `kwargs` through to the scipy API.
3. The code may throw an error if, during the solving phase, a value of $\mathbf{x}$ outside of the slowness domain is queried (this happens as the ray is exiting the domiain). This can be avoided by passing as an argument to the solver `max_step=n*max_dx` where $n$ can be up to 5. 
4. The parameter $\lambda_f$ can be set arbitrarily high as the solution will terminate if the ray leaves the domain. 


## Specifiying the takeoff angle

The takeoff direction is specified with a pair of angles $(\alpha, \beta)$. $\alpha$ is the inclination angle and $\beta$ is the azimuth. $\alpha=0$ when $\mathbf{p}$ is aligned with the $z$ direction (downward) and increases counterclockwise. $\beta=0$ when $\mathbf{p}$ is aligned with the $x$ axis and increases clockwise.

## Tracing multiple rays

Multiple rays can be traced by passing a list of tuples for the source coordinates and takeoff angles e.g. `srcs=[(0,0,0), (0,0,0)]` `angles=[(0,45), (0,30)]`. Tracing can be run in parallel by passing `parallel=True` and `n_procs` to the runner. For example

```
out = tracer.run(parallel=True, n_procs=6)
```
