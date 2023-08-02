
#include <iostream>

#include "CollisionIntegral.h"
#include "CollElem.h"

// Monte Carlo integration
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte_vegas.h>

namespace gslWrapper {
     // Helpers for GSL integration routines. note that we cannot pass a member function by reference,
     // so we dodge this in the wrapper by passing a reference to the object whose function we want to integrate

     struct gslFunctionParams {
          // Pointer to the object whose member function we are integrating
          CollisionIntegral4* pointerToObject;
          // Additional (non-integration variable) parameters needed to evaluate the integrand
          CollisionIntegral4::IntegrandParameters integrandParameters;
     };

     // pp should be of gslFunctionParams type
     inline double integrandWrapper(double* intVars, size_t dim, void* pp) {
          
          // GSL requires size_t as input: the logic is that dim should be used to check that intVars is correct size.
          // But I'm a badass and just do this:
          (void)dim;

          gslFunctionParams* params = static_cast<gslFunctionParams*>(pp);

          double p2 = intVars[0];
          double phi2 = intVars[1];
          double phi3 = intVars[2];
          double cosTheta2 = intVars[3];
          double cosTheta3 = intVars[4];
          return params->pointerToObject->calculateIntegrand(p2, phi2, phi3, cosTheta2, cosTheta3, params->integrandParameters);
     } 
}


// This calculates the full collision integral C[m,n; j,k]. NOTE: has to be thread safe!!
std::array<double, 2> CollisionIntegral4::evaluate(int m, int n, int j, int k) {

     IntegrandParameters integrandParameters = initializeIntegrandParameters(m, n, j, k); 

     // Integral dimensions
     constexpr size_t dim = 5;
     // Define the integration limits for each variable: {p2, phi2, phi3, cosTheta2, cosTheta3}
     double integralLowerLimits[dim] = {0.0, 0.0, 0.0, -1., -1.}; // Lower limits
     double integralUpperLimits[dim] = {maxIntegrationMomentum, 2.0*constants::pi, 2.0*constants::pi, 1., 1.}; // Upper limits

     // Construct parameter wrapper struct
     gslWrapper::gslFunctionParams gslWrapper;
     gslWrapper.pointerToObject = this;
     gslWrapper.integrandParameters = integrandParameters;

     gsl_monte_function G;
     G.f = &gslWrapper::integrandWrapper;
     G.dim = dim;
     G.params = &gslWrapper;

     // How many Monte Carlo iterations
     size_t calls = 50000;
     size_t warmupCalls = 0.1*calls;
     double mean = 0.0;
     double error = 0.0;

     // Unique Monte Carlo state for this integration (NB: keep this here: thread safety)
     gsl_monte_vegas_state* gslState = gsl_monte_vegas_alloc(dim);

     // Warmup?!?
     gsl_monte_vegas_integrate(&G, integralLowerLimits, integralUpperLimits, dim, warmupCalls, this->gslRNG, gslState, &mean, &error);

     /* TODO: GSL instructions on the Vegas routine did a bunch of "warmup" iterations before the actual "convergence" runs. 
     Need to understand what the optimal usage of Vegas routines is. */

     gsl_monte_vegas_integrate(&G, integralLowerLimits, integralUpperLimits, dim, calls, this->gslRNG, gslState, &mean, &error);
     
     gsl_monte_vegas_free(gslState);

     return std::array<double, 2>( {mean, error} );
}



double CollisionIntegral4::calculateIntegrand(double p2, double phi2, double phi3, double cosTheta2, double cosTheta3, 
        const IntegrandParameters &integrandParameters) const
     {
     
     double fullIntegrand = 0.0;

     // Variables from the input structure, redefine locally here for convenience
     const int m = integrandParameters.m;
     const int n = integrandParameters.n;

     const double p1 = integrandParameters.p1;
     const double pZ1 = integrandParameters.pZ1;
     const double pPar1 = integrandParameters.pPar1;

     // Sines
     const double sinTheta2 = std::sin(std::acos(cosTheta2));
     const double sinTheta3 = std::sin(std::acos(cosTheta3));
     const double sinPhi2 = std::sin(phi2);
     const double sinPhi3 = std::sin(phi3);
     
     // Cosines
     const double cosPhi2 = std::cos(phi2);
     const double cosPhi3 = std::cos(phi3);

     
     // SLOPPY: Create 3-vectors, but I just use FourVectors with vanishing 0-component
     FourVector FV1dummy(0.0, pPar1, 0.0, pZ1);
     FourVector FV2dummy(0.0, p2*sinTheta2*cosPhi2, p2*sinTheta2*sinPhi2, p2*cosTheta2);
     // 'p3Hat': like p3, but normalized to 1. We will fix its magnitude using a delta(FV4^2 - msq4)
     const FourVector FV3Hat(0.0, sinTheta3*cosPhi3, sinTheta3*sinPhi3, cosTheta3);

     // dot products of 3-vectors. Need minus sign here since I'm hacking this with 4-vectors with (+1 -1 -1 -1) metric
     const double p1p2Dot = -1.0 * FV1dummy*FV2dummy;
     const double p1p3HatDot = -1.0 * FV1dummy*FV3Hat;
     const double p2p3HatDot = -1.0 * FV2dummy*FV3Hat;

     // Kinematics differs for each collision process since the masses are generally different
     // TODO optimize so that if everything is ultrarelativistic, calculate kinematic factors only once
     for (CollElem<4> collElem : collisionElements) {

          // Vacuum masses squared
          std::array<double, 4> massSquared;
          for (int i=0; i<4; ++i) {
               massSquared[i] = collElem.particles[i].getVacuumMassSquared();
          }

          // Energies: Since p3 is not fixed yet we only know E1, E2 
          double E1 = std::sqrt(p1*p1 + massSquared[0]);
          double E2 = std::sqrt(p2*p2 + massSquared[1]);


          //------------------------------- TODO move this bit elsewhere?
          
          // Now def. function g(p3) = FV4(p3) - msq4 and express the delta function in 
          // terms of roots of g(p3) and delta(p3 - p3root); p3 integral becomes trivial.
          // Will need to solve a quadratic equation; some helper variables for it:
          double Q = massSquared[0] + massSquared[1] + massSquared[2] - massSquared[3];
          double kappa = Q + 2.0 * (E1*E2 - p1p2Dot);
          double eps = 2.0 * (E1 + E2);
          double delta = 2.0 * (p1p3HatDot + p2p3HatDot);
          
          auto funcG = [&](double p3) {
               double m3sq = massSquared[2];
               return kappa + delta*p3 - eps * sqrt(p3*p3 + m3sq);
          };

          // Quadratic eq. A p3^2 + B p3 + C = 0, take positive root(s)
          double A = delta*delta - eps*eps;
          double B = 2.0 * kappa * delta;
          double C = kappa*kappa - eps*eps*massSquared[2];
          // Roots of g(p3):
          double discriminant = B*B - 4.0*A*C;
          double root1 = 0.5 * (-B - sqrt(discriminant)) / A;
          double root2 = 0.5 * (-B + sqrt(discriminant)) / A;

          // TODO:
          // 1) Check for possible singularities
          // 2) Check that the root(s) satisfy the original eq. with square root
          // 3) Is it possible to have 2 positive roots??

          //-------------------------------
          std::array<double, 2> rootp3({root1, root2});

          // Now proceed to fix remaining 4-momenta
          for (double p3 : rootp3) if (p3 >= 0.0) {

               if (std::abs(funcG(p3)) > 1e-8) {
                    std::cerr << "! Invalid root in CollisionIntegral4::calculateIntegrand \n";
               }

               double E3 = std::sqrt(p3*p3 + massSquared[2]);

               // Fix 4-momenta for real this time
               FourVector FV1 = FV1dummy;
               FV1[0] = E1;
               FourVector FV2 = FV2dummy;
               FV2[0] = E2;
               FourVector FV3 = p3*FV3Hat;
               FV3[0] = E3;

               // momentum conservation fixes P4 
               const FourVector FV4 = FV1 + FV2 - FV3;

               if (FV4.energy() < 0.0) {
                    std::cerr << "! Negative energy E4: " << FV4.energy() << "\n";
                    continue;
               }

               // Kinematic prefactor, ie. everything from integration measure
               double kinPrefac = 1.0;
               // Avoid spurious singularity at zero momenta, zero mass
               if (std::abs(massSquared[1]) < massSquaredLowerBound) kinPrefac *= p2;
               else kinPrefac *= p2*p2 / E2;

               if (std::abs(massSquared[2]) < massSquaredLowerBound) kinPrefac *= p3;
               else kinPrefac *= p3*p3 / E3;
               // additional factor from delta(g(p3))
               double gDer = 0.0;
               if (std::abs(massSquared[2]) < massSquaredLowerBound) gDer = delta - eps;
               else gDer = delta - eps * p3 / E3;
               
               kinPrefac /= std::abs(gDer);

               // Calculate deltaF's for all momenta.
               // In our spectral approach this means that we replace deltaF with 
               // Tm(rhoZ) Tn(rhoPar), where Tm, Tn are appropriate basis polynomials
               const double TmTn_p1 = integrandParameters.TmTn_p1;
               const double TmTn_p2 = polynomialBasis.TmTn(m, n, FV2);
               const double TmTn_p3 = polynomialBasis.TmTn(m, n, FV3);
               const double TmTn_p4 = polynomialBasis.TmTn(m, n, FV4);
               const std::array<double, 4> TmTn_p({TmTn_p1, TmTn_p2, TmTn_p3, TmTn_p4});

               std::array<double, 4> deltaF;
                for (int i=0; i<4; ++i) {
                    if (collElem.particles[i].isInEquilibrium()) {
                         // particle species in equilibrium, so deltaF = 0
                         deltaF[i] = 0.0;
                    } else {
                         deltaF[i] = TmTn_p[i];
                    }
                }

               // Now evaluate |M|^2 P[f] for this CollElem
               std::array<FourVector, 4> fourMomenta({FV1, FV2, FV3, FV4});
               double integrand = collElem.evaluate( fourMomenta, deltaF );

               integrand *= kinPrefac;

               fullIntegrand += integrand;
          } // end p3 : rootp3

     } // end collElem : collisionElements
     
     // Common numerical prefactor
     constexpr double PI = constants::pi;
     constexpr double pi2Pow5 = (2.0*PI) * (2.0*PI) * (2.0*PI) * (2.0*PI) * (2.0*PI);
     return fullIntegrand / pi2Pow5 / 8.0;
}
