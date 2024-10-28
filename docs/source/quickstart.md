# Quickstart

In order to compute collision integrals with **WallGoCollision** you have to
1) Create a `PhysicsModel` object (particle and model parameter definitions).
2) Load in symbolic matrix elements to the model.
3) Use the model to create a `CollisionTensor` object and pass it the size of your momentum/polynomial grid.
4) Calling `computeIntegralsAll()` from your `CollisionTensor` computes all required integrals for your particle content and grid size. The results can be stored in binary `.hdf5` format and later loaded into the **WallGo** Boltzmann solver.

This page demonstrates **WallGoCollision** use for a Boltzmann system of one fermion and one boson species, and contains also an arbitrary number of identical "light" fermions that are assumed to remain in thermal equilibrium but appear in matrix elements. The following is written using the Python API. See [the examples directory](https://github.com/Wall-Go/WallGoCollision/tree/main/examples) for code examples using the C++ API. 

The following assumes you have the **WallGoCollision** Python module [installed](install.md). Load it in your Python program with `import WallGoCollision`.

## Defining the model

Model definition is done by filling in a ModelDefinition helper object and passing it to the PhysicsModel constructor. Start by defining your model parameters:
```
modelDef = WallGoCollision.ModelDefinition()
gs = 1.228 # Corresponds to the QCD coupling
modelDef.defineParameter("gs", gs) 
```
Parameters must be given as (name, value) pairs. The value must be floating-point type or convertible to such, eg. complex numbers are not allowed. In the above, "gs" is the parameter name, and any appearance of "gs" in the symbolic matrix elements will be replaced with the numeric value during collision integration.

Next we define our particle content using the ParticleDefinition class. Each particle species must have a unique name (string) as well as a unique integer identifier ("particle index"). The index is used to associate loaded matrix elements with the correct particles. Additionally, you must specify the particle statistics type (boson or fermion), and whether the species is assumed to remain in thermal equilibrium. The latter can be used to reduce the number of collision integrations for models containing particle species for which deviations from equilibrium are negligible. Here we define a "top quark" and a "gluon" as out-of-equilibrium particles, and a generic "light quark" that is kept in equilibrium.
```
topQuark = WallGoCollision.ParticleDescription()
topQuark.name = "Top"
topQuark.index = 0
topQuark.type = WallGoCollision.EParticleType.eFermion
topQuark.bInEquilibrium = False
modelDef.defineParticleSpecies(topQuark)

gluon = WallGoCollision.ParticleDescription()
gluon.name = "Gluon"
gluon.index = 1
gluon.type = WallGoCollision.EParticleType.eBoson
gluon.bInEquilibrium = False
modelDef.defineParticleSpecies(gluon)

# Generic light quark, assumed to remain in thermal equilibrium
lightQuark = WallGoCollision.ParticleDescription()
lightQuark.name = "LightQuark"
lightQuark.index = 2
lightQuark.type = WallGoCollision.EParticleType.eFermion
lightQuark.bInEquilibrium = True
modelDef.defineParticleSpecies(lightQuark)
```
The above particle species are all defined as "ultrarelativistic" (the default behavior). In **WallGoCollision**, an ultrarelativistic particle species $a$ means that its energy-momentum dispersion relation is approximated as $E_a \approx |p_a|$ in collision integrals. This allows for heavy optimizations and should be preferred whenever the approximation makes sense for your particles. For going beyond the ultrarelativistic approximation, see [here (TODO! empty link for now)]().

Finally, we should define numerical values for masses that appear in propagators of matrix elements. In **WallGoCollision**, these masses are treated as model parameters analogous to the "gs" parameter defined above, ie. a propagator mass is just another user-defined symbol that can appear in matrix elements and must be given a numberical value. In the matrix elements for this example, the mass-square of a quark propagator is denoted by "mq2" and the mass-square of a gluon propagator is "mg2". For this model we approximate them as their asymptotic QCD thermal masses (same for all quark species):
```
modelDef.defineParameter("mq2", gs**2 / 3.0)
modelDef.defineParameter("mg2", gs**2)
```
> [!NOTE]
> Any dimensionful parameters must be given in units of the temperature. Therefore the above corresponds to eg. $\frac13 g_s^2 T^2$ for the quark mass-square.

To finalize our model definition we create a concrete `PhysicsModel` based on the information in our `modelDef` object:
```
collisionModel = WallGoCollision.PhysicsModel(modelDef)
```
Once created, the particle and parameter content of a `PhysicsModel` is fixed and cannot be changed. You are still allowed to update values of existing parameters by calling `model.updateParameter(name, newValue)`. Runtime changes to a PhysicsModel will automatically propagate to `CollisionTensor` objects created from it.

## Loading matrix elements

The next step is to specify what collision processes are allowed in the model. This is done by loading symbolic (ie. math expressions) matrix elements for the relevant processes from an external file. Finding matrix elements for a given particle physics model is an exercise in field theory and not within capabilities of **WallGoCollision**. The companion Mathematica package [**WallGoMatrix**](https://github.com/Wall-Go/WallGoMatrix) can be used to generate matrix elements for arbitrary models, directly in the format required to **WallGoCollision**.

**WallGoCollision** supports matrix element parsing in JSON format (needs `.json` file extension) or from generic text files (legacy option). The JSON format is generally recommended for better validation; the `.json` matrix elements relevant for this example are available [here](https://github.com/Wall-Go/WallGoCollision/tree/main/examples/MatrixElements/MatrixElements_QCD.json). A legacy `.txt` version is available [here](https://github.com/Wall-Go/WallGoCollision/tree/main/examples/MatrixElements/MatrixElements_QCD.txt), note that the formatting must be exactly as shown in the file for parsing to work. Each matrix element must specify indices of external particles participating in the collision process, and a symbolic expression that is typically a function of Mandelstam variables (denoted by reserved symbols "_s", "_t", "_u") and of any user-specified symbols (such as "gs" in this example).

Matrix elements are loaded to a PhysicsModel as follows:
```
# Replace with your own path. The second argument is optional and specifies whether the loaded matrix elements should be printed for logging purposes
bSuccess: bool = collisionModel.loadMatrixElements("MatrixElements/MatrixElements_QCD.json", True)
```
This parses the file and finds matrix elements for particles that the model is aware of (the particle index attribute is used for this purpose). Any processes involving undefined particle indices are ignored (currently we also print a warning for each ignored matrix element). At this stage the program also attempts to numerically evaluate each loaded matrix element in order to verify that all required symbols were included in the model definition, if this is not the case an error message is printed, load is aborted and `loadMatrixElements()` returns False.

[TODO: should explain how we group off-eq particles to pairs and require specific indexing in matrix elements]

## `CollisionTensor` creation and evaluation

The `CollisionTensor` class contains collision integrals in an unevaluated but otherwise ready form. Once matrix elements have been loaded we can create a `CollisionTensor` object as follows:
```
# Modify according to your needs. Using a trivially small grid to make this example run fast
gridSize = 5
# Type hinting included for completeness
collisionTensor: WallGoCollision.CollisionTensor = collisionModel.createCollisionTensor(gridSize)
```
`gridSize` is the number of basis polynomials on your momentum/polynomial grid as used by the **WallGo** pipeline, and must be integer. It can be freely changed later with `collisionTensor.changePolynomialBasisSize(newSize)`.

Before starting integrations it can be useful to configure the integrator to best suit your needs. A handful of settings is available in the `IntegrationOptions` class, including error tolerances for the Monte Carlo integration and the upper limit on momentum integration. You can then pass your modified `IntegrationOptions` object to `CollisionTensor` as `collisionTensor.setIntegrationOptions(yourOptionsObject)`. Here we skip this step and use the default settings. Similarly, the `CollisionTensorVerbosity` class can be used to configure verbosity settings, such as progress reporting and time estimates. Set it as `collisionTensor.setIntegrationVerbosity(yourVerbosityObject)`.

To perform the actual integrations, use
```
results: WallGoCollision.CollisionTensorResult = collisionTensor.computeIntegralsAll()
```
This is typically a very long running functions for nontrivial models and grid sizes. For grid size $N$ the number of required integrals scales as $(N-1)^4$. The return type CollisionTensorResult is a wrapper around 4D numerical array that contains statistical error estimates of the Monte Carlo integration.  Once finished, you can save the results in `.hdf5` format as follows:
```
# Replace with your output directory
results.writeToIndividualHDF5("CollisionOutDir/")
```
This creates a separate `.hdf5` for each pair of out-of-equilibrium particles in the model in the specified output directory. The data can then be readily loaded into [**WallGo**](https://wallgo.readthedocs.io).
