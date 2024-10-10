# WallGoCollision quickstart

In order to compute collision integrals with **WallGoCollision** you have to
1) Create a `PhysicsModel` object (particle and model parameter definitions).
2) Load in symbolic matrix elements to the model.
3) Use the model to create a `CollisionTensor` object and pass it the size of your momentum/polynomial grid.
4) Calling `computeIntegralsAll()` from your `CollisionTensor` computes all required integrals for your particle content and grid size. The results can be stored in binary `.hdf5` format and later loaded into the **WallGo** Boltzmann solver.

This page demonstrates **WallGoCollision** use up for a Boltzmann system of one fermion and one boson species, and contains also an arbitrary number of identical "light" fermions that are assumed to remain in thermal equilibrium but appear in matrix elements. The following is written using the Python API. See [the examples directory](../examples) for code examples using the C++ API. 

The following assumes you have the **WallGoCollision** Python module [installed](../README.md).

