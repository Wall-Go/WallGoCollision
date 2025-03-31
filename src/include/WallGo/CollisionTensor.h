#pragma once

#include <vector>
#include <map>
#include <string>
#include <chrono>
#include <filesystem>
#include <utility> // std::pair
#include <memory>
#include <cstdint>

#include "Common.h"
#include "ModelParameters.h"
#include "CollisionElement.h"
#include "ParticleSpecies.h"
#include "CollisionIntegral.h"
#include "ResultContainers.h"
#include "ModelObserver.h"


namespace wallgo
{

class PhysicsModel;


/* Class CollisionTensor -- Contains collision integrals in an unevaluated form
and provides an interface for evaluating them.
Uses an observer pattern to automatically detect changes in the PhysicsModel instance that created it.
*/
class CollisionTensor : public IModelObserver
{

public:

    CollisionTensor(PhysicsModel* creator);
    CollisionTensor(PhysicsModel* creator, size_t basisSize);
    // Destructor that calls unregisterObserver() on the creator
    virtual ~CollisionTensor();

    /* Rule of three: need custom copy/assignment operators because of the custom dtor that handles unregistration from model.
    * Not ideal, but these ensure that a copied CollisionTensor is also registered as a model observer.
    */
    CollisionTensor(const CollisionTensor& other);
    CollisionTensor& operator=(const CollisionTensor& other);

    /* Configures default integration options that are used by computeIntegralsForPair() and related functions
    if no IntegrationOptions object is passed when calling them. */
    void setDefaultIntegrationOptions(const IntegrationOptions& options);

    /* Configures default integration verbosity that is used by computeIntegralsForPair() and related functions
    if no CollisionIntegralVerbosity object is passed when calling them. */
    void setDefaultIntegrationVerbosity(const CollisionTensorVerbosity& verbosity);

    // Change basis size used by the polynomial grid. This is a fast operation, does not require rebuild of stored integral objects
    void changePolynomialBasisSize(size_t newBasisSize);

    /* Adds a new collision integral object for an off - eq particle pair.
    Currently we always have just one CollisionIntegral4 for each pair, so this will override any existing integral for that pair.
    NB: does NOT perform sensibility checks on the input integral.*/
    void addCollisionIntegral(const ParticleNamePair& particleNames, const CollisionIntegral4& inIntegral);

    // Get pointer to specified CollisionIntegral4. Will be null if the integral is not found.
    CollisionIntegral4* getIntegralForPair(const ParticleNamePair& particleNames);

    /* Calculate collision integrals for particle pair (particle1, particle2) everywhere on the grid. */
    CollisionResultsGrid computeIntegralsForPair(
        const std::string& particle1,
        const std::string& particle2,
        const IntegrationOptions& options,
        const CollisionTensorVerbosity& verbosity);

    /* Calculate collision integrals for particle pair (particle1, particle2) everywhere on the grid.
    Uses default integration options. */
    CollisionResultsGrid computeIntegralsForPair(
        const std::string& particle1,
        const std::string& particle2,
        const CollisionTensorVerbosity& verbosity);

    /* Calculate collision integrals for particle pair (particle1, particle2) everywhere on the grid.
    Uses default verbosity settings */
    CollisionResultsGrid computeIntegralsForPair(
        const std::string& particle1,
        const std::string& particle2,
        const IntegrationOptions& options);


    /* Calculate collision integrals for particle pair (particle1, particle2) everywhere on the grid.
    Uses default integration options and verbosity settings. */
    CollisionResultsGrid computeIntegralsForPair(
        const std::string& particle1,
        const std::string& particle2);

    /* Calculates all integrals previously initialized with setupCollisionIntegrals(). */
    CollisionTensorResult computeIntegralsAll(const IntegrationOptions& options, const CollisionTensorVerbosity& verbosity);

    /* Calculates all integrals previously initialized with setupCollisionIntegrals().
    Uses default verbosity. */
    CollisionTensorResult computeIntegralsAll(const IntegrationOptions& options);

    /* Calculates all integrals previously initialized with setupCollisionIntegrals().
    Uses default options. */
    CollisionTensorResult computeIntegralsAll(const CollisionTensorVerbosity& verbosity);

    /* Calculates all integrals previously initialized with setupCollisionIntegrals().
    Uses default options and verbosity. */
    CollisionTensorResult computeIntegralsAll();

    IntegrationResult computeSingleIntegral(
        const std::string& particle1,
        const std::string& particle2,
        const GridPoint& gridPoint,
        const IntegrationOptions& options);

    IntegrationResult computeSingleIntegral(
        const std::string& particle1,
        const std::string& particle2,
        const GridPoint& gridPoint);

    // Count how many independent collision integrals we have. Scales as N^4 * M^2, N = grid size, M = number of off-eq particles
    size_t countIndependentIntegrals() const;

    // IModelObserver interface
    virtual void handleModelChange(const ModelChangeContext& context) override;
    virtual void handleModelDestruction() override;
    // ~IModelObserverInterface

private:

    PhysicsModel* mObservingModel = nullptr;

    size_t mBasisSize;

    IntegrationOptions mDefaultIntegrationOptions;
    CollisionTensorVerbosity mDefaultVerbosity;

    // Holds collision integrals so that they can be reused
    std::map<ParticleNamePair, CollisionIntegral4> mCachedIntegrals;

};

} // namespace
