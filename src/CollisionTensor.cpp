#include "CollisionTensor.h"

#include <array>
#include <fstream>
#include <sstream>
#include <regex> // Reading matrix elements from file
#include <algorithm> // std::remove_if
#include <chrono>
#include <filesystem>

namespace wallgo
{

// Function for this file only. Processes string of form "M[a,b,c,d] -> some funct" and stores in the arguments
void interpretMatrixElement(const std::string &inputString, std::vector<size_t> &indices, std::string &mathExpression)
{
    // First split the string by "->""
    std::vector<std::string> tokens(2);
    
    std::string delimiter = "->";
    std::string lhs = inputString.substr(0, inputString.find(delimiter));

    // RHS
    mathExpression = inputString.substr(lhs.length() + delimiter.length());

    // remove whitespaces from lhs to avoid weirdness
    lhs.erase(std::remove_if(lhs.begin(), lhs.end(), isspace), lhs.end());

    // ---- Extract the abcd indices from M[a,b,c,d]
    std::size_t start = lhs.find('[');
    std::size_t end = lhs.find(']');

    // Ensure '[' and ']' are found and the start position is before the end position
    if (start != std::string::npos && end != std::string::npos && start < end) 
    {
        std::string values = lhs.substr(start + 1, end - start - 1);

        // Use stringstream to tokenize and extract integers
        std::istringstream ss(values);
        indices.clear();
        indices.reserve(4);
        size_t num;

        while (ss >> num) 
        {
            indices.push_back(num);

            // Check for the ',' separator and ignore it
            if (ss.peek() == ',')
            {
                ss.ignore();
            }
        }
    }
}

CollisionTensor::CollisionTensor()
{
    // Set default options

    mBasisSize = 1;

    mDefaultIntegrationOptions.calls = 50000;
    mDefaultIntegrationOptions.maxIntegrationMomentum = 20;
    mDefaultIntegrationOptions.maxTries = 50;
    mDefaultIntegrationOptions.relativeErrorGoal = 1e-1;
    mDefaultIntegrationOptions.absoluteErrorGoal = 1e-8;
    mDefaultIntegrationOptions.bOptimizeUltrarelativistic = true;
    mDefaultIntegrationOptions.bIncludeStatisticalErrors = true;

    mDefaultVerbosity = CollisionTensorVerbosity();
}

CollisionTensor::CollisionTensor(size_t basisSize)
    : CollisionTensor()
{
    changePolynomialBasisSize(basisSize);
}

void CollisionTensor::defineParticle(const ParticleSpecies &particle)
{
    const std::string name = particle.getName();
    
    if (isParticleRegistered(particle))
    {
        std::cout << "Particle " << name << " already registered with CollisionTensor, doing nothing\n";
        return;
    }

    // TODO let the user specify the index themselves?
    particles.push_back(std::make_shared<ParticleSpecies>(particle));

    particleIndex[name] = particles.size() - 1;

    if (!particle.isInEquilibrium())
    {
        outOfEqParticles.push_back(particles.back());
    }
}


void CollisionTensor::setDefaultIntegrationOptions(const IntegrationOptions& options)
{
    mDefaultIntegrationOptions = options;
}

void CollisionTensor::setDefaultIntegrationVerbosity(const CollisionTensorVerbosity& verbosity)
{
    mDefaultVerbosity = verbosity;
}


void CollisionTensor::updateParticleMasses(const std::map<std::string, double> &msqVacuum, const std::map<std::string, double> &msqThermal)
{
    // TODO proper validation. The .at() will crash if the key is not found 
    for (const auto &[name, msq] : msqVacuum)
    {
        const size_t idx = particleIndex.at(name);
        particles[idx]->setVacuumMassSquared(msq);
    }

    for (const auto &[name, msq] : msqThermal)
    {
        const size_t idx = particleIndex.at(name);
        particles[idx]->setThermalMassSquared(msq);
    }
}

void CollisionTensor::changePolynomialBasisSize(size_t newBasisSize)
{
    mBasisSize = newBasisSize;
    for (auto& [key, integral] : mCachedIntegrals)
    {
        integral.changePolynomialBasis(newBasisSize);
    }
}


void CollisionTensor::defineVariable(const std::string &varName, double value)
{
    if (mModelParameters.contains(varName))
    {
        std::cerr << "Error: variable " << varName << " has already been defined for this CollisionTensor instance" << std::endl;
    }
    else
    {
        mModelParameters.addOrModifyParameter(varName, value);
    }
}

void CollisionTensor::defineVariables(const std::map<std::string, double> &variables)
{
    for (const auto& [name, value] : variables)
    {
        defineVariable(name, value);
    }
}

void CollisionTensor::setVariable(const std::string &varName, double value)
{
    if (!mModelParameters.contains(varName))
    {
        std::cerr << "Error: can't change value of undefined variable " << varName << "\n";
        return;
    }
    mModelParameters.addOrModifyParameter(varName, value);

    // Sync collision elements    
    for (auto& [_, integral] : mCachedIntegrals)
    {
        integral.updateModelParameter(varName, value);
    }
}

void CollisionTensor::setVariables(const std::map<std::string, double> &newValues)
{
    for (const auto &[key, val] : newValues)
    {
        setVariable(key, val);
    }
}

CollisionIntegral4 CollisionTensor::setupCollisionIntegral(const std::shared_ptr<ParticleSpecies>& particle1, const std::shared_ptr<ParticleSpecies>& particle2, 
    const std::string &matrixElementFile, size_t inBasisSize, bool bVerbose)
{
    std::vector<CollElem<4>> collisionElements = parseMatrixElements(particle1->getName(), particle2->getName(), matrixElementFile, bVerbose);

    ParticleNamePair namePair(particle1->name, particle2->name);
    CollisionIntegral4 collisionIntegral(inBasisSize, namePair);

    for (const CollElem<4> &elem : collisionElements)
    {
        collisionIntegral.addCollisionElement(elem);
    }

    return collisionIntegral;
}


void CollisionTensor::setupCollisionIntegrals(const std::filesystem::path& matrixElementFile, bool bVerbose)
{
    clearIntegralCache();

    for (const std::shared_ptr<ParticleSpecies>& particle1 : outOfEqParticles)
    for (const std::shared_ptr<ParticleSpecies>& particle2 : outOfEqParticles)
    {
        const std::string name1 = particle1->getName();
        const std::string name2 = particle2->getName();

        const auto namePair = std::make_pair(name1, name2);

        CollisionIntegral4 newIntegral = setupCollisionIntegral(particle1, particle2, matrixElementFile.string(), mBasisSize, bVerbose);

        mCachedIntegrals.insert({ namePair, newIntegral });
    }
}

void CollisionTensor::clearIntegralCache()
{
    mCachedIntegrals.clear();
}

CollisionResultsGrid CollisionTensor::computeIntegralsForPair(
    const std::string& particle1,
    const std::string& particle2,
    const IntegrationOptions& options,
    const CollisionTensorVerbosity& verbosity)
{
    const auto pairName = std::make_pair(particle1, particle2);
    if (mCachedIntegrals.count(pairName) < 1)
    {
        std::cerr << "No cached collision integral found for particle pair: ("
            << particle1 << ", " << particle2 << ")" << std::endl;
        
        assert(false && "Particle pair not found");
        
        return CollisionResultsGrid(ParticleNamePair("Unknown", "Unknown"), CollisionMetadata());
    }

    return mCachedIntegrals.at(pairName).evaluateOnGrid(options, verbosity);
}

CollisionResultsGrid CollisionTensor::computeIntegralsForPair(
    const std::string& particle1,
    const std::string& particle2,
    const CollisionTensorVerbosity& verbosity)
{
    return computeIntegralsForPair(particle1, particle2, mDefaultIntegrationOptions, verbosity);
}

CollisionResultsGrid CollisionTensor::computeIntegralsForPair(
    const std::string& particle1,
    const std::string& particle2,
    const IntegrationOptions& options)
{
    return computeIntegralsForPair(particle1, particle2, options, mDefaultVerbosity);
}

CollisionResultsGrid CollisionTensor::computeIntegralsForPair(
    const std::string& particle1,
    const std::string& particle2)
{
    return computeIntegralsForPair(particle1, particle2, mDefaultIntegrationOptions, mDefaultVerbosity);
}

CollisionTensorResult CollisionTensor::computeIntegralsAll(const IntegrationOptions& options, const CollisionTensorVerbosity& verbosity)
{
    if (mCachedIntegrals.size() < 1)
    {
        std::cerr << "Warning: computeIntegralsAll() called on empty CollisionTensor object" << std::endl;
        // return empty result
        return CollisionTensorResult();
    }

    CollisionTensorResult result(mCachedIntegrals.size());

    std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
    if (verbosity.bPrintElapsedTime)
    {
        startTime = std::chrono::steady_clock::now();
    }

    size_t i = 0;
    for (auto& [namePair, integral] : mCachedIntegrals)
    {
        result.mData[i] = integral.evaluateOnGrid(options, verbosity);

        std::cout << particlePairToString(namePair) << " done\n";

        ++i;
    }

    if (verbosity.bPrintElapsedTime)
    {
        auto endTime = std::chrono::steady_clock::now();
        auto elapsedTime = endTime - startTime;
        auto elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds>(elapsedTime).count();

        std::cout << "\nAll done, took " << elapsedSeconds << " seconds.\n" << std::endl;
    }

    return result;
}

CollisionTensorResult CollisionTensor::computeIntegralsAll(const IntegrationOptions& options)
{
    return computeIntegralsAll(options, mDefaultVerbosity);
}

CollisionTensorResult CollisionTensor::computeIntegralsAll(const CollisionTensorVerbosity& verbosity)
{
    return computeIntegralsAll(mDefaultIntegrationOptions, verbosity);
}

CollisionTensorResult CollisionTensor::computeIntegralsAll()
{
    return computeIntegralsAll(mDefaultIntegrationOptions, mDefaultVerbosity);
}


/*
CollisionTensorResult CollisionTensor::calculateAllIntegrals(bool bVerbose)
{
    CollisionTensorResult result(mCachedIntegrals.size());

    // Initialize progress tracking 
    totalIntegralCount = countIndependentIntegrals(mBasisSize, outOfEqParticles.size());
    computedIntegralCount = 0;
    bFinishedInitialProgressCheck = false;
    startTime = std::chrono::steady_clock::now();


    // make rank 2 tensor that mixes out-of-eq particles (each element is a collision integral, so actually rank 6, but the grid indices are irrelevant here)

    size_t i = 0;
    for (auto & [namePair, integral] : mCachedIntegrals)
    {
        std::chrono::steady_clock::time_point pairStartTime = std::chrono::steady_clock::now();

        const std::string name1 = namePair.first;
        const std::string name2 = namePair.second;

        result.mData[i] = integral.evaluateOnGrid(mDefaultIntegrationOptions, mDefaultVerbosity);

        i++;

        std::cout << particlePairToString(namePair) << " done\n";

        // Create a new HDF5 file. H5F_ACC_TRUNC means we overwrite the file if it exists
        const std::string fileNameBase = "collisions_" + name1 + "_" + name2 + ".hdf5";
        std::filesystem::path outputPath(outputDirectory);
        
        outputPath = outputPath / fileNameBase;

        // How long did this all take
        std::chrono::duration<double> duration = std::chrono::steady_clock::now() - pairStartTime;
        double seconds = duration.count();
        int hours = static_cast<int>(seconds / 3600);
        // leftover mins
        int minutes = static_cast<int>(seconds - hours * 3600 / 60);
        std::cout << "[" << name1 << ", " << name2 << "] done in " << hours << "h " << minutes << "min." << std::endl;
    }

    return result;
}
*/

size_t CollisionTensor::countIndependentIntegrals() const
{
    size_t res = 0;
    for (const auto& [_, integral] : mCachedIntegrals)
    {
        res += integral.countIndependentIntegrals();
    }
    return res;
}


CollElem<4> CollisionTensor::makeCollisionElement(const std::string &particleName2, const std::vector<size_t> &indices, 
    const std::string &expr, const std::vector<std::string> &symbols)
{
    assert(indices.size() == 4);

    const std::shared_ptr<ParticleSpecies> p1 = particles[indices[0]];
    const std::shared_ptr<ParticleSpecies> p2 = particles[indices[1]];
    const std::shared_ptr<ParticleSpecies> p3 = particles[indices[2]];
    const std::shared_ptr<ParticleSpecies> p4 = particles[indices[3]];
    
    CollElem<4> collisionElement( { p1, p2, p3, p4} );

    std::map<std::string, double> variables;
    for (const std::string& s : symbols)
    {
        if (!mModelParameters.contains(s))
        {
            std::cerr << "CollisionTensor error: symbol " << s << " not found in modelParameters map, can't parse matrix element " << expr << "\n";
            return collisionElement;
        }
        variables[s] = mModelParameters.getParameterValue(s);
    }

    collisionElement.matrixElement.initSymbols(variables);
    // Parse expr as a math function
    collisionElement.matrixElement.setExpression(expr);

    // Set deltaF flags: in general there can be 4 deltaF terms but we only take the ones with deltaF of particle2
    for (size_t i = 0; i < collisionElement.bDeltaF.size(); ++i) 
    {
        collisionElement.bDeltaF[i] = (collisionElement.particles[i]->getName() == particleName2);
    }
    
    return collisionElement;
}


std::vector<CollElem<4>> CollisionTensor::parseMatrixElements(
    const std::string &particleName1,
    const std::string &particleName2,
    const std::string &inMatrixElementFile,
    bool bVerbose)
{
    // Just for logging 
    const std::string pairName = "[" + particleName1 + ", " + particleName2 + "]"; 

    // @todo should actually check if these are in outOfEquilibriumParticles vector
    if (particleIndex.count(particleName1) < 1 || particleIndex.count(particleName2) < 1 ) 
    {
        std::cerr << "Error: particles missing from list! Was looking for out-of-eq particles " << pairName << "\n";
        exit(9);
    }

    if (bVerbose) std::cout << "\n" <<"Parsing matrix elements for off-equilibrium pair " << pairName << "\n";
    

    /* Since we don't have per-expression symbol lists implemented in the files,
    here I just define each symbol in our modelParameters list as a variable. This wastes memory.
    Alternatively, could try to implicitly find free symbols with the parser, but this seems complicated. */
    std::vector<std::string> symbols = mModelParameters.getParameterNames();

    std::ifstream file(inMatrixElementFile);

    /*M_ab->cd, with a = particle1 and at least one of bcd is particle2.Suppose that for each out - of - eq pair, Mathematica gives these in form
    M[a, b, c, d] -> (some symbolic expression), where abcd are integer indices that need to match our ordering in particleIndex map
    Here we parse the lhs to extract indices, then parse the rhs as a math expression (function of s,t,u and other symbols that we infer).
    For each of these we make a CollElem<4> object with correct deltaF structure which we infer from the indices */

    std::vector<CollElem<4>> collisionElements;

    if (!file.is_open()) {
        std::cerr << "!!! Error: Failed to open matrix element file " << inMatrixElementFile << std::endl;
        return collisionElements;
    }

    /* Now use regex to read all lines of form M[...] -> expr
    For each line we check if the first particle index matches that of particleName1
    and require that at least one other index matches that of particleName2.
    This is not optimal because we end up reading the full file for each off-eq pair.
    */
    std::string line;
    std::string expr;
    std::vector<size_t> indices;
    indices.resize(4);
    
    while (std::getline(file, line)) {
        if (std::regex_search(line, std::regex("M\\[.*\\] -> (.*)"))) {

            /* Big TODO. Change matrix element file format so that 
            1. It's easier to parse without regex hacks. Eg: JSON format
            2. Each matrix element could be associated with a list of symbols needed to evaluate it.
            This would make it possible to safely define just enough symbols needed for each matrix element.
            */
            
            interpretMatrixElement(line, indices, expr);
            
            if (indices[0] != particleIndex[particleName1]) continue;
            if ( std::find(indices.begin(), indices.end(), particleIndex[particleName2]) == indices.end() ) continue;
 
            CollElem<4> newElem = makeCollisionElement(particleName2, indices, expr, symbols);
            collisionElements.push_back(newElem);

            if (bVerbose)
            {
                std::cout << "Loaded matrix element:\n";
                std::cout << line << "\n";
            }
        }
    }

    file.close();

    if (bVerbose) std::cout << std::endl;

    return collisionElements;
}

/*
void CollisionTensor::reportProgress()
{
    if (totalIntegralCount > 0)
    {
        double elapsedSeconds = elapsedTime.count();
        double timePerIntegral = elapsedSeconds / computedIntegralCount;
        double timeRemaining = timePerIntegral * (totalIntegralCount - computedIntegralCount);

        int percentage = int(100 * static_cast<double>(computedIntegralCount) / totalIntegralCount);
        std::cout << "Integral progress: " << computedIntegralCount << " / " << totalIntegralCount << " (" << percentage << "%). "; 
        std::cout << "Estimated time remaining: " << std::floor(timeRemaining / 3600) << "h " << (int(timeRemaining) % 3600 ) / 60 << "min" << std::endl;
    }
    bFinishedInitialProgressCheck = true;
}
*/

bool CollisionTensor::isParticleRegistered(const ParticleSpecies& particle)
{
    return particleIndex.count(particle.getName()) > 0;
}


} // namespace