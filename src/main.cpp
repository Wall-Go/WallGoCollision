/*******/

#include <iostream>
#include <cmath>
#include <chrono>

#include <cstring>
#include <tclap/CmdLine.h> // command line arguments

#include "Common.h"
#include "CollisionManager.h"
#include "CollElem.h"
#include "CollisionIntegral.h"
#include "hdf5Interface.h"
#include "gslWrapper.h"
#include "MatrixElement.h"
#include "ConfigParser.h"

/* Test/example routine, calculates all collision integrals with QCD interactions. 
The structure here illustrates how the same could be done from Python with arbitrary inputs */ 
void collisionsQCD(uint N) {

	const double gs = 1.2279920495357861;

    //**** Masses squared. These need to be in units of temperature, ie. (m/T)^2 **//
	// Thermal
	const double mq2 = 0.251327; // quark
	const double mg2 = 3.01593; // SU(3) gluon
	// Vacuum, @TODO if needed
	const double msqVacuum = 0.0;

	const bool bUltraRelativistic = true;

	// For now take just the top to be out-of-eq
    ParticleSpecies topQuark("top", EParticleType::FERMION, false, msqVacuum, mq2, bUltraRelativistic);
    ParticleSpecies gluon("gluon", EParticleType::BOSON, false, msqVacuum, mg2, bUltraRelativistic);
	ParticleSpecies lightQuark("quark", EParticleType::FERMION, true, msqVacuum, mq2, bUltraRelativistic);

	// Main control object
	CollisionManager collisionManager;
	collisionManager.addParticle(topQuark);
	collisionManager.addParticle(gluon);
	collisionManager.addParticle(lightQuark);

	collisionManager.addCoupling(gs);

	collisionManager.calculateCollisionIntegrals(N);

}


int main(int argc, char** argv) {

	//--------------- How this works

	/* class Collision: This is a control class that collects all functionality in one place, 
	* basically a replacement for main function. */

	/* class CollElem : Describes a matrix element with fixed external particles (ordering matters!).
	* This is the object that calculates |M|^2 and the statistical 'population factor' P once the external momenta are fixed.
	* We need a separate CollElem object for each scattering process that contributes to the collision integral (tt->gg, tg->tg, tq->tq are separate CollElems)
	* Currently the matrix elements are just hard coded, in a more realistic setting they would probably be read from elsewhere. */

	/* class ParticleSpecies : Quite self-explanatory. Contains info about particle statistics, masses and whether the particle species stays in equilibrium, etc.
	* These are given as inputs to CollElem when constructing CollElem objects.
	* The particle name property is important as it is used to read in the correct matrix element (this needs improvement in the future). */

	/* class CollisionIntegral4 : This describes the whole 2-by-2 collision integral for a given particle type (top quark in this case).
	* IE: this object calculates eq. (A1) in 2204.13120, with delta f replace with Chebyshev polynomials.
	* Therefore the class it needs to know the size of our polynomial basis (N) and CollElem objects that make up the integral.
	* The class calculates 5D integrals with same integration variables as Benoit had. See notes in the private repo. */

	/* Currently the interface between CollisionIntegral4 and CollElem is not optimal and there is some redundancy
	* in how the CollisionIntegral4 obtains particle masses etc. This needs to be improved in next version before we
	* get started with generic matrix elements */

	//---------------

	// basis size, default value
	int basisSizeN = 20;
    // config file, default name
    std::string configFileName = "CollisionDefaults.ini";


	// ---- Setup and parse command line args
	try 
	{
		TCLAP::CmdLine cmd("WallGo Collision program (standalone)", ' ', /*version*/"1.0");

		TCLAP::ValueArg<int> basisSizeArg("n", "basisSize", "Polynomial basis size", /*required*/false, /*default*/11, "Positive integer");
		cmd.add(basisSizeArg);

		TCLAP::ValueArg<std::string> configFileArg("c", "configFile", "Path to config file", false, "CollisionDefaults.ini", "String");
		cmd.add(configFileArg);

		cmd.parse(argc, argv);

		basisSizeN = basisSizeArg.getValue();
		configFileName = configFileArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
		exit(1);
 	}


	if (basisSizeN < 1) {
		std::cerr << "Invalid basis size N = " << basisSizeN << "\n";
		return 2;
	}

	// Load config
	ConfigParser& config = ConfigParser::get();

	if (config.load(configFileName)) {
		std::cout << "Read config:\n";
		config.printContents();
		std::cout << std::endl;
	} else {
		return 3;
	}

	std::cout << "Running with basis size "<< basisSizeN << "\n";


	gslWrapper::initializeRNG();

	collisionsQCD(basisSizeN);


/*
	//-------------------- Measure wall clock time

	std::cout << "Running speed test: integral C[2,1,1,1]\n";
	auto startTime = std::chrono::steady_clock::now();
	collInt.evaluate(2, 1, 1, 1);

	auto endTime = std::chrono::steady_clock::now();

	auto elapsedTime = endTime - startTime;
	// Convert the elapsed time to milliseconds
	auto elapsedTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(elapsedTime).count();

	// How long for all collision integrals
	auto totalTime = elapsedTime * nCollisionTerms;

	auto hours = std::chrono::duration_cast<std::chrono::hours>(totalTime).count();
	auto minutes = std::chrono::duration_cast<std::chrono::minutes>(totalTime).count() % 60;

	std::cout << "Test done, took " << elapsedTimeMs << "ms\n";

	std::cout << "Estimated time-per-thread for all " << nCollisionTerms << " collision integrals: "
			<< hours << " hours " << minutes << " minutes\n";

	//--------------------
*/

	// Cleanup
	gslWrapper::clearRNG();
	return 0;
}
