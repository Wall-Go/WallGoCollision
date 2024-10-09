# WallGoCollision documentation

**WallGoCollision** is a scientific library for computing Boltzmann collision integrals for use with the Python package [**WallGo**](https://github.com/Wall-Go/WallGo). The library is written in C++17 and can be compiled as a Python module for seamless interoperation with WallGo. Collision integrations are typically by far the most computationally intensive part of the **WallGo** wall velocity pipeline. **WallGoCollision** performs these integrations using native C++, allowing for performance optimizations that would not be possible in a pure Python package.

> [!IMPORTANT]
> **WallGoCollision** is still in beta and may undergo large changes without notice. Use with care.

## What it does

The purpose of **WallGoCollision** is to compute integrals of form
```math
\mathcal{C}_a[\delta f] \propto \sum_{bcd} \int \frac{d^3\mathbf{p}_2 d^3\mathbf{p}_3 d^3\mathbf{p}_4}{E_2 E_3 E_4} \delta^4(p_1 + p_2 - p_3 - p_4) |M_{ab \rightarrow cd}(p_1, p_2; p_3, p_4)|^2 \mathcal{P}_{ab \rightarrow cd}[\delta f],
```
appearing in the collision operator of the relativistic Boltzmann equation of particle $a$. 
