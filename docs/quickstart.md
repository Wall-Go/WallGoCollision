# WallGoCollision quickstart

In order to compute collision integrals with **WallGoCollision** you have to
1) Create a `PhysicsModel` object (particle and model parameter definitions).
2) Load in symbolic matrix elements to the model.
3) Use the model to create a `CollisionTensor` object and pass it the size of your momentum/polynomial grid.
4) Calling `computeIntegralsAll()` from your `CollisionTensor` computes all required integrals for your particle content and grid size. The results can be stored in binary `.hdf5` format and later loaded into the **WallGo** Boltzmann solver.

This page demonstrates **WallGoCollision** use up for a Boltzmann system of one fermion and one boson species, and contains also an arbitrary number of identical "light" fermions that are assumed to remain in thermal equilibrium but appear in matrix elements. The following is written using the Python API. See [the examples directory](../examples) for code examples using the C++ API. 

The following assumes you have the **WallGoCollision** Python module [installed](../README.md). Load it in your Python program with `import WallGoCollision`.

## Defining the model

Model definition is done by filling in a ModelDefinition helper object and passing it to the PhysicsModel constructor. Start by defining your model parameters:
```
modelDef = WallGoCollision.ModelDefinition()
gs = 1.228 # Corresponds to the QCD coupling
modelDef.defineParameter("gs", gs) 
```
Parameters must be given as (name, value) pairs. In the above, "gs" is the parameter name, and any appearance of "gs" in the symbolic matrix elements will be replaced with the numeric value during collision integration.

Next we define our particle content using the ParticleDefinition class. Each particle species must have a unique name (string) as well as a unique integer identifier ("particle index"). The index is used to associate loaded matrix elements with the correct particles. Additionally, you must specify the particle statistics type (boson or fermion), and whether the species is assumed to remain in thermal equilibrium. The latter can be used to reduce the number of collision integrations for models containing particle species for which deviations from equilibrium are negligible. Here we define a "top quark" and a "gluon" as out-of-equilibrium particles, and a generic "light quark" that is kept in equilibrium.
```
topQuark = WallGoCollision.ParticleDescription()
topQuark.name = "top"
topQuark.index = 0
topQuark.type = WallGoCollision.EParticleType.eFermion
topQuark.bInEquilibrium = False
modelDef.defineParticleSpecies(topQuark)

gluon = WallGoCollision.ParticleDescription()
gluon.name = "gluon"
gluon.index = 1
gluon.type = WallGoCollision.EParticleType.eBoson
gluon.bInEquilibrium = False
modelDef.defineParticleSpecies(gluon)

# Generic light quark, assumed to remain in thermal equilibrium
lightQuark = WallGoCollision.ParticleDescription()
lightQuark.name = "lightQuark"
lightQuark.index = 2
lightQuark.type = WallGoCollision.EParticleType.eFermion
lightQuark.bInEquilibrium = True
modelDef.defineParticleSpecies(lightQuark)
```
The above particle species are all defined as "ultrarelativistic" (the default behavior). In **WallGoCollision**, an ultrarelativistic particle species $a$ means that its energy-momentum dispersion relation is approximated as $E_a \approx |p_a|$ in collision integrals. This allows for heavy optimizations and should be preferred whenever the approximation makes sense for your particles. For relaxing this approximation, see [here (TODO! empty link for now)]().

Finally, we should define numerical values for masses that appear in propagators of matrix elements. In **WallGoCollision**, these masses are treated as model parameters analogous to the "gs" parameter defined above, ie. a propagator mass is just another user-defined symbol that can appear in matrix elements and must be given a numberical value. In the matrix elements for this example, the mass-square of a quark propagator is denoted by "mq2" and the mass-square of a gluon propagator is "mg2". For this model we approximate them as their asymptotic QCD thermal masses (same for all quark species):
```
modelDef.defineParameter("mq2", gs**2 / 3.0)
modelDef.defineParameter("mg2", gs**2)
```

