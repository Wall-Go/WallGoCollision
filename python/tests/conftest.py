import pytest
import pathlib

import WallGoCollision

fileDir = pathlib.Path(__file__).parent.resolve()


@pytest.fixture()
def collisionModelQCD() -> WallGoCollision.PhysicsModel:
    """SU(3) model with 1 off-eq quark, off-eq gluon and 5 in-eq light quarks.
    All ultrarelativistic. Also loads matrix elements"""

    defs = WallGoCollision.ModelDefinition()

    gs = 1.2279920495357861

    defs.defineParameter("gs", 1.2279920495357861)
    defs.defineParameter("mq2", gs**2 / 6.0)
    defs.defineParameter("mg2", 2.0 * gs**2)

    topQuark = WallGoCollision.ParticleDescription()
    topQuark.name = "top"
    topQuark.index = 0
    topQuark.type = WallGoCollision.EParticleType.eFermion
    topQuark.bInEquilibrium = False
    topQuark.bUltrarelativistic = True

    defs.defineParticleSpecies(topQuark)

    gluon = WallGoCollision.ParticleDescription()
    gluon.name = "gluon"
    gluon.index = 1
    gluon.type = WallGoCollision.EParticleType.eBoson
    gluon.bInEquilibrium = False
    gluon.bUltrarelativistic = True
    defs.defineParticleSpecies(gluon)

    lightQuark = topQuark
    lightQuark.bInEquilibrium = True
    lightQuark.name = "light quark"
    lightQuark.index = 2
    defs.defineParticleSpecies(lightQuark)

    model = WallGoCollision.PhysicsModel(defs)

    model.loadMatrixElements(
        str(fileDir / "MatrixElements/MatrixElements_QCD.json"), False
    )

    yield model

@pytest.fixture()
def collisionTensorQCD(
    collisionModelQCD: WallGoCollision.PhysicsModel
) -> WallGoCollision.CollisionTensor:
    """Collision tensor for QCD model, top/gluon are off eq.
    Initialized to N=11."""

    collisionTensor = collisionModelQCD.createCollisionTensor(11)
    yield collisionTensor

@pytest.fixture()
def resultsQCD_N3(
    collisionTensorQCD: WallGoCollision.CollisionTensor,
) -> WallGoCollision.CollisionTensorResult:
    """Collision integration results at N=3 for QCD model"""

    WallGoCollision.setSeed(0)

    collisionTensorQCD.changePolynomialBasisSize(3)

    options = WallGoCollision.IntegrationOptions()
    options.maxIntegrationMomentum = 20
    options.absoluteErrorGoal = 1e-8
    options.relativeErrorGoal = 1e-1
    options.maxTries = 50
    options.calls = 50000
    options.bIncludeStatisticalErrors = True

    collisionTensorQCD.setIntegrationOptions(options)

    resultTensor = collisionTensorQCD.computeIntegralsAll()

    yield resultTensor
