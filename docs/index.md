# WallGoCollision documentation

**WallGoCollision** is a scientific library for computing Boltzmann collision integrals for use with the Python package [**WallGo**](https://github.com/Wall-Go/WallGo). The library is written in C++17 and can be compiled as a Python module for seamless interoperation with WallGo. Collision integrations are typically by far the most computationally intensive part of the **WallGo** wall velocity pipeline. **WallGoCollision** performs these integrations using native C++, allowing for performance optimizations that would not be possible in a pure Python package.

> [!IMPORTANT]
> **WallGoCollision** is still in beta and may undergo large changes without notice. Use with care.

## What it does

The purpose of **WallGoCollision** is to compute integrals of form
```math
\mathcal{C}_a[\delta f] \propto \sum_{bcd} \int \frac{d^3\mathbf{p}_2 d^3\mathbf{p}_3 d^3\mathbf{p}_4}{E_2 E_3 E_4} \delta^4(p_1 + p_2 - p_3 - p_4) |M_{ab \rightarrow cd}(p_1, p_2; p_3, p_4)|^2 \mathcal{P}_{ab \rightarrow cd}[\delta f],
```
appearing in the collision operator of the relativistic Boltzmann equation for particle species $a$. **WallGo**, and **WallGoCollision** by extension, uses the spectral approach of Cline & Laurent ([paper](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.106.023501)) to Boltzmann integration, meaning the integrals must be evaluated on a 4D grid of momenta and suitable basis polynomials. The integrations are performed using importance-sampled Monte Carlo. Integrals on different grid points are independent and trivially computable in parallel. We support OpenMP for parallel evaluation. Currently **WallGoCollision** is restricted to the linearized version of collision integrals. See [the physics page](./physics.md) for more details on the physics.

The user must provide relevant particle content for their model, matrix elements $M_{ab \rightarrow cd}$ of their choice in a symbolic form, and numerical values for model parameters appearing in the matrix elements. To get started, see [the quickstart page](./quickstart.md), [code examples in C++](../examples) or [code examples using the Python interface on **WallGo** repo](https://github.com/Wall-Go/WallGo/tree/main/Models).


## Limitations
