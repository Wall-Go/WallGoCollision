# Configuring the integrator

Collision integrations in **WallGoCollision** are performed using the Vegas Monte Carlo algorithm (G. P. Lepage, J. Comput. Phys. 27 (1978) 192) which is an adaptive probabilistic integrator. We use the [GSL implementation](https://www.gnu.org/software/gsl/doc/html/montecarlo.html#vegas) of the algorithm. **WallGoCollision** applies it to 5-dimensional integrals of form described in the [physics documentation](./physics.md). The set of integration variables consists of four angles (or cosines of thereof) and one radial variable, corresponding to magnitude of the incoming $\mathbf{p}_2$ momentum.

[WIP: the defaults and other details here are subject to change, and may be moved to dedicated API documentation]

The `IntegrationOptions` class collects parameters that can be adjusted for granular control over the integrator:
- `maxIntegrationMomentum` specifies the upper bound on the radial integration variable. Must be given in units of the temperature. A value of order 10 is often sufficient because of Boltzmann suppression at high momenta. Default is 20.

- `calls` sets the number of function evaluations between convergence checks. After `calls` evaluations, we check the statistical error and finish if it is within the specified `relativeErrorGoal` or `absoluteErrorGoal`. If not, we repeat the process with another `calls` evaluations. We also implement a check based on chi-squared of the Vegas algorithm, requiring that $\chi^2$ per degrees of freedom is comparable to 1 (threshold: $|\chi^2 / d - 1| \leq 0.5$). If the $\chi^2$ check fails, the current error estimate is deemed unreliable and we perform another `calls` function evaluations. Default is 50000.

- `maxTries` sets the maximum number of time the process described above is repeated before giving up. Therefore there will be a maximum of `calls * maxTries` function evaluations. Default is 50.

- `relativeErrorGoal` and `absoluteErrorGoal` set error tolerances used in the convergence check.

These options are unique for each `CollisionTensor` object. You can update them by filling in a `IntegrationOptions` object and passing it as `CollisionTensor.setIntegrationOptions(myOptions)` (in the Python API). The C++ API additionally exposes overloads of `CollisionTensor::computeIntegralsAll()` that accept an `IntegrationOptions` object as an input.
