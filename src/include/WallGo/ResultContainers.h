#pragma once

#include <filesystem>
#include <unordered_map>

#include "Common.h"
#include "hdf5Interface.h"
#include "IntegrationOptions.h"
#include "ModelParameters.h"

namespace wallgo
{

/** The linearized Boltzmann equations for deltaF (deviation from equilibrium) is of form
    (L_ab + T^2 C_ab) deltaF_b = S_a
where a,b label particle species, L is the Liouville operator and S is a source term.
C_ab is the (dimensionless) collision operator which we focus on here.

WallGo follows the spectral approach of 2204.13120 in which deltaF is expanded as a series of "basis polynomials" T_m.
The collision operator is rank-2 in particle indices (a,b), rank-2 in polynomial indices (m,n)
and rank-2 in momentum indices (j, k) (momentum components of particle 'a' that are parallel and perpendicular to the wall).
In total C is a rank-6 tensor, assuming the collision integrals are not explicitly position-dependent.

Note that the last assumption fails if the (dimensionless) C_ab contain explicit T-dependence,
which occurs eg. if vacuum masses are included in matrix elements. We do not consider this more complicated case.

Here we describe this as a 2D array of 4D objects so that particle indices are separated from the (polynomial, momentum) "grid" indices.
*/


// Holds metadata about collision integrations
struct CollisionMetadata
{
    size_t basisSize = 1;
    bool bStatisticalErrors = false;
    std::string basisName = "Unknown";
    std::string integrator = "Unknown";
    uint64_t seed = 0;
    // How many threads were used. This is useful info because multithreading affects order of RNG
    uint32_t numThreads = 1;

    // Copy of the integration options used when producing this collision data
    IntegrationOptions usedIntegrationOptions;

    // Evaluation time in seconds
    uint32_t timeSpent = 0;

    // Copy of model parameters that were used
    ModelParameters modelParameters;
};

/* Rank 4 tensor that holds collision integration results on the grid for (particle1, particle2) pair
Range of indices is: [2, N] for m, [1, N-1] for njk.
Attempting to access grid points with invalid indices results in failed assert (Debug builds)
or a crash (Release builds). */
class CollisionResultsGrid
{
public:

    CollisionResultsGrid() : mParticlePair("Unknown", "Unknown"), mElementsPerDimension(0) {}
    CollisionResultsGrid(const ParticleNamePair& particlePair, const CollisionMetadata& metadata);

    inline bool hasStatisticalErrors() const { return mMetadata.bStatisticalErrors; }
    
    size_t getBasisSize() const { return mMetadata.basisSize; }
    ParticleNamePair getParticleNamePair() const { return mParticlePair; }

    // ---- Read-write accessors

    double valueAt(const GridPoint& gridPoint) const;
    double& valueAt(const GridPoint& gridPoint);

    double errorAt(const GridPoint& gridPoint) const;
    double& errorAt(const GridPoint& gridPoint);

    // Write (value, error) pair. Error will only be written if the CollisionResultsGrid object was created with bIncludeStatisticalErrors = true
    void updateValue(const GridPoint& gridPoint, double newValue, double newError);

    /* Write array contents to a HDF5 file. This always overrides file if it exists.
    Dataset name for the integration results will be "particle1, particle2" with particle names read from mDescription.
    Dataset name for integration errors will be "particle1, particle2 errors"
    Return value is false if something goes wrong. */
    bool writeToHDF5(
        const std::filesystem::path& filePath,
        bool bWriteErrors = true) const;

    // Sets all elements to a specified value
    void fillWithValue(double newValue, double newError);

private:

    ParticleNamePair mParticlePair;
    
    size_t mElementsPerDimension;  // equals basisSize - 1

    /* The actual data, stored as 4D array with index order : (momentum1, momentum2, poly1, poly2).
    This ordering matches that used in the Python-side WallGo */
    Array4D mData;
    // Statistical errors of integrations. Has large memory cost, so we keep this empty if not needed
    Array4D mErrors;

    CollisionMetadata mMetadata;

    void initData();

    bool validateGridPoint(const GridPoint& gridPoint) const;

    // Allow CollisionIntegral4 to update metadata on the go
    friend class CollisionIntegral4;
};

struct CollisionTensorResult
{
public:

    CollisionTensorResult() {}
    CollisionTensorResult(size_t numParticlePairs) : mData(numParticlePairs) {}

    /* Writes the contents to disk, so that each particle pair is written into a separate HDF5 file.
    * In practice this just calls CollisionResultsGrid::writeToHDF5 for each particle pair.
    * Return value is false if something went wrong. */
    bool writeToIndividualHDF5(
        const std::filesystem::path& outDirectory,
        bool bWriteErrors = true) const;

    /* Writes the contents to current working directory, so that each particle pair is written into a separate HDF5 file.
    * In practice this just calls CollisionResultsGrid::writeToHDF5 for each particle pair.
    * Return value is false if something went wrong. */
    bool writeToIndividualHDF5(bool bWriteErrors = true) const;

    // Returns pointer to collision integration results of the specified particle pair. Can be nullptr if the pair is not found
    CollisionResultsGrid* getResultsForParticlePair(const ParticleNamePair& particlePair);

    // Returns pointer to collision integration results of the specified particle pair. Can be nullptr if the pair is not found
    CollisionResultsGrid* getResultsForParticlePair(const std::string& particleName1, const std::string& particleName2);

    // CollisionResultsGrid is not copyable, so use std::vector
    std::vector<CollisionResultsGrid> mData;
};

} // namespace
