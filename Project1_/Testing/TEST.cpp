#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <iostream>
#include "../sampler.h"
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include "../WaveFunctions/simplegaussian.h"
#include "../WaveFunctions/simplegaussian_numerical.h"
#include "../WaveFunctions/ellipticalgaussian.h"
#include "../WaveFunctions/ellipticalgaussiananalytical.h"
#include "../Hamiltonians/hamiltonian.h"
#include "../Hamiltonians/harmonicoscillator.h"
#include "../Hamiltonians/harmonicoscillatorinteracting.h"
#include "../InitialStates/initialstate.h"
#include "../InitialStates/randomuniform.h"
#include "../Math/random.h"
#include <math.h>


class System* setUpSystemSimpleGaussian(double alph, int numParticles, int numDims)
{
    int numberOfDimensions  = numDims;
    int numberOfParticles   = numParticles;
    int numberOfSteps       = (int) pow(2.0, 19.0);
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = alph;          // Variational parameter.
    double stepLength       = 0.5;          // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used for equilibration.
    double D                = 0.5;          // Diffusion constant. Will use D = 0.5.
    double timestep         = 0.3;          // Importance sampling time step. 

    // INITIALIZE A NEW SYSTEM

    System* system = new System();

    // CHOOSE HAMILTONIAN:

    system->setHamiltonian              (new HarmonicOscillator(system, omega)); double a = 0.0;

    system->setWaveFunction             (new SimpleGaussian(system, alpha)); 

    // SET PARAMS:

    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setDiffusionConstant        (D);

    return system;
}

class System* setUpSystemNumericalSimpleGaussian(double alph , int numParticles, int numberOfDims)
{
    int numberOfDimensions  = numberOfDims;
    int numberOfParticles   = numParticles;
    int numberOfSteps       = (int) pow(2.0, 19.0);
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = alph;          // Variational parameter.
    double stepLength       = 0.5;          // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used for equilibration.
    double D                = 0.5;          // Diffusion constant. Will use D = 0.5.
    double timestep         = 0.3;          // Importance sampling time step. 

    // INITIALIZE A NEW SYSTEM

    System* system = new System();

    // CHOOSE HAMILTONIAN:

    system->setHamiltonian              (new HarmonicOscillator(system, omega)); double a = 0.0;

    // CHOOSE WAVEFUNCTION:

    system->setWaveFunction             (new SimpleGaussianNumerical(system, alpha)); 

    // SET SYSTEM PARAMS:

    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setDiffusionConstant        (D);

    return system;
}

double testGetAlphaSimpleGaussian(double alpha)
{

    class System* system = setUpSystemSimpleGaussian(alpha, 10, 10);

    std::vector<double> params = system->getWaveFunction()->getParameters();

    return params[0];
}

double testEnergyCalculationMetropolisSimpleGaussian(double alpha,  int numParticles, int numDims)
{
    class System* system = setUpSystemSimpleGaussian(alpha, numParticles, numDims);

    system->runMetropolisSteps(0.5*1e6);

    double energy = system->getSampler()->getEnergy();

    return energy;

}

double testEnergyCalculationImportanceSamplingSimpleGaussian(double alpha, double stepLength, int numParticles, int numDims)
{
    class System* system = setUpSystemSimpleGaussian(alpha, numParticles, numDims);

    system->runImportanceSamplingSteps(pow(2,19), stepLength);

    double energy = system->getSampler()->getEnergy();

    return energy;
}

double testEnergyCalculationImportanceSamplingNumericalSimpleGaussian(double alpha, double stepLength, int numParticles, int numDims)
{
    class System* system = setUpSystemNumericalSimpleGaussian(alpha, numParticles, numDims);

    system->runImportanceSamplingSteps(pow(2,19), stepLength);

    double energy = system->getSampler()->getEnergy();

    return energy;
}

double testEnergyCalculationMetropolisNumericalSimpleGaussian(double alpha, int numParticles, int numDims)
{
    class System* system = setUpSystemNumericalSimpleGaussian(alpha, numParticles, numDims);

    system->runMetropolisSteps(pow(2,19));

    double energy = system->getSampler()->getEnergy();

    return energy;
}

double calculateExactEnergySimpleGaussian(double alpha,  int numParticles, int numDims)
{
    return (alpha/2 + 1/(8*alpha))*numDims*numParticles;
}

double testGradientDescentSimpleGaussian(double initialValue, double stepLength, int numParticles, int numDims)
{
    class System* system = setUpSystemSimpleGaussian(0.5,numParticles, numDims);

    double alpha_minimum = system->runGradientDescent(stepLength, initialValue, 60, 1e5, 5*1e-7);

    return alpha_minimum;
}


TEST_CASE("Test simple gaussian for returning the alpha value", "[Simple Gaussian Get Alpha]")
{
    REQUIRE(testGetAlphaSimpleGaussian(0.5) == 0.5);
    REQUIRE(testGetAlphaSimpleGaussian(0.1) == 0.1);
    REQUIRE(testGetAlphaSimpleGaussian(10)  == 10);

}

TEST_CASE("Test simple gaussian energy calculation with Metropolis brute-force", "[Simple Gaussian Energy Calculation Metropolis Brute-Force]")
{
    double tol = 0.05;

    REQUIRE(abs(testEnergyCalculationMetropolisSimpleGaussian(0.4, 1, 1)  - calculateExactEnergySimpleGaussian(0.4, 1,1)) < tol);
    REQUIRE(abs(testEnergyCalculationMetropolisSimpleGaussian(0.5, 1, 1)  - calculateExactEnergySimpleGaussian(0.5, 1,1)) < tol);
    REQUIRE(abs(testEnergyCalculationMetropolisSimpleGaussian(0.6, 1, 1)  - calculateExactEnergySimpleGaussian(0.6, 1,1)) < tol);

    REQUIRE(abs(testEnergyCalculationMetropolisSimpleGaussian(0.4, 10, 1) - calculateExactEnergySimpleGaussian(0.4, 10,1)) < tol);
    REQUIRE(abs(testEnergyCalculationMetropolisSimpleGaussian(0.5, 10, 1) - calculateExactEnergySimpleGaussian(0.5, 10,1)) < tol);
    REQUIRE(abs(testEnergyCalculationMetropolisSimpleGaussian(0.6, 10, 1) - calculateExactEnergySimpleGaussian(0.6, 10,1)) < tol);

    REQUIRE(abs(testEnergyCalculationMetropolisSimpleGaussian(0.4, 10, 3) - calculateExactEnergySimpleGaussian(0.4, 10,3)) < tol);
    REQUIRE(abs(testEnergyCalculationMetropolisSimpleGaussian(0.5, 10, 3) - calculateExactEnergySimpleGaussian(0.5, 10,3)) < tol);
    REQUIRE(abs(testEnergyCalculationMetropolisSimpleGaussian(0.6, 10, 3) - calculateExactEnergySimpleGaussian(0.6, 10,3)) < tol); 

}

TEST_CASE("Test simple gaussian energy calculation with Metropolis Importance Sampling", "[Simple Gaussian Energy Calculation Metropolis Importance Sampling]")
{
    double tol = 0.005;

    REQUIRE(abs(testEnergyCalculationImportanceSamplingSimpleGaussian(0.4, 0.4, 1, 1)  - calculateExactEnergySimpleGaussian(0.4, 1,1)) < tol);
    REQUIRE(abs(testEnergyCalculationImportanceSamplingSimpleGaussian(0.5, 0.4, 1, 1)  - calculateExactEnergySimpleGaussian(0.5, 1,1)) < tol);
    REQUIRE(abs(testEnergyCalculationImportanceSamplingSimpleGaussian(0.6, 0.4, 1, 1)  - calculateExactEnergySimpleGaussian(0.6, 1,1)) < tol);

    REQUIRE(abs(testEnergyCalculationImportanceSamplingSimpleGaussian(0.4, 0.4, 10, 1) - calculateExactEnergySimpleGaussian(0.4, 10,1)) < tol);
    REQUIRE(abs(testEnergyCalculationImportanceSamplingSimpleGaussian(0.5, 0.4, 10, 1) - calculateExactEnergySimpleGaussian(0.5, 10,1)) < tol);
    REQUIRE(abs(testEnergyCalculationImportanceSamplingSimpleGaussian(0.6, 0.4, 10, 1) - calculateExactEnergySimpleGaussian(0.6, 10,1)) < tol);

    REQUIRE(abs(testEnergyCalculationImportanceSamplingSimpleGaussian(0.4, 0.4, 10, 3) - calculateExactEnergySimpleGaussian(0.4, 10,3)) < tol);
    REQUIRE(abs(testEnergyCalculationImportanceSamplingSimpleGaussian(0.5, 0.4, 10, 3) - calculateExactEnergySimpleGaussian(0.5, 10,3)) < tol);
    REQUIRE(abs(testEnergyCalculationImportanceSamplingSimpleGaussian(0.6, 0.4, 10, 3) - calculateExactEnergySimpleGaussian(0.6, 10,3)) < tol); 

}

TEST_CASE("Test numerical simple Gaussian energy calculation with Metropolis brute-force", "[Numerical Simple Gaussian Energy Calculation Metropolis Brute Force]")
{
    double tol = 0.05; 

    REQUIRE(abs(testEnergyCalculationMetropolisNumericalSimpleGaussian(0.4, 1, 1)  - calculateExactEnergySimpleGaussian(0.4, 1,1)) < tol);
    REQUIRE(abs(testEnergyCalculationMetropolisNumericalSimpleGaussian(0.5, 1, 1)  - calculateExactEnergySimpleGaussian(0.5, 1,1)) < tol);
    REQUIRE(abs(testEnergyCalculationMetropolisNumericalSimpleGaussian(0.6, 1, 1)  - calculateExactEnergySimpleGaussian(0.6, 1,1)) < tol);

    REQUIRE(abs(testEnergyCalculationMetropolisNumericalSimpleGaussian(0.4, 10, 1) - calculateExactEnergySimpleGaussian(0.4, 10,1)) < tol);
    REQUIRE(abs(testEnergyCalculationMetropolisNumericalSimpleGaussian(0.5, 10, 1) - calculateExactEnergySimpleGaussian(0.5, 10,1)) < tol);
    REQUIRE(abs(testEnergyCalculationMetropolisNumericalSimpleGaussian(0.6, 10, 1) - calculateExactEnergySimpleGaussian(0.6, 10,1)) < tol);

    REQUIRE(abs(testEnergyCalculationMetropolisNumericalSimpleGaussian(0.4, 10, 3) - calculateExactEnergySimpleGaussian(0.4, 10,3)) < tol);
    REQUIRE(abs(testEnergyCalculationMetropolisNumericalSimpleGaussian(0.5, 10, 3) - calculateExactEnergySimpleGaussian(0.5, 10,3)) < tol);
    REQUIRE(abs(testEnergyCalculationMetropolisNumericalSimpleGaussian(0.6, 10, 3) - calculateExactEnergySimpleGaussian(0.6, 10,3)) < tol); 
}

TEST_CASE("Test numerical simple Gaussian energy calculation with Metropolis Importance Sampling", "[Numerical Simple Gaussian Energy Calculation Metropolis Importance Sampling]")
{
    double tol = 0.05; 

    REQUIRE(abs(testEnergyCalculationImportanceSamplingNumericalSimpleGaussian(0.4, 0.4, 1, 1)  - calculateExactEnergySimpleGaussian(0.4, 1,1)) < tol);
    REQUIRE(abs(testEnergyCalculationImportanceSamplingNumericalSimpleGaussian(0.5, 0.4, 1, 1)  - calculateExactEnergySimpleGaussian(0.5, 1,1)) < tol);
    REQUIRE(abs(testEnergyCalculationImportanceSamplingNumericalSimpleGaussian(0.6, 0.4, 1, 1)  - calculateExactEnergySimpleGaussian(0.6, 1,1)) < tol);

    REQUIRE(abs(testEnergyCalculationImportanceSamplingNumericalSimpleGaussian(0.4, 0.4, 10, 1) - calculateExactEnergySimpleGaussian(0.4, 10,1)) < tol);
    REQUIRE(abs(testEnergyCalculationImportanceSamplingNumericalSimpleGaussian(0.5, 0.4, 10, 1) - calculateExactEnergySimpleGaussian(0.5, 10,1)) < tol);
    REQUIRE(abs(testEnergyCalculationImportanceSamplingNumericalSimpleGaussian(0.6, 0.4, 10, 1) - calculateExactEnergySimpleGaussian(0.6, 10,1)) < tol);

    REQUIRE(abs(testEnergyCalculationImportanceSamplingNumericalSimpleGaussian(0.4, 0.4, 10, 3) - calculateExactEnergySimpleGaussian(0.4, 10,3)) < tol);
    REQUIRE(abs(testEnergyCalculationImportanceSamplingNumericalSimpleGaussian(0.5, 0.4, 10, 3) - calculateExactEnergySimpleGaussian(0.5, 10,3)) < tol);
    REQUIRE(abs(testEnergyCalculationImportanceSamplingNumericalSimpleGaussian(0.6, 0.4, 10, 3) - calculateExactEnergySimpleGaussian(0.6, 10,3)) < tol); 
}

TEST_CASE("Test Gradient Descent Simple Gaussian", "[Gradient Descent Simple Gaussian]")
{
    double tol = 0.005;

    REQUIRE(abs(testGradientDescentSimpleGaussian(0.5, 0.1, 1, 1) - 0.5) < tol);
    REQUIRE(abs(testGradientDescentSimpleGaussian(0.5, 0.1, 1, 3) - 0.5) < tol);
    REQUIRE(abs(testGradientDescentSimpleGaussian(0.5, 0.01, 10, 1) - 0.5) < tol);
    REQUIRE(abs(testGradientDescentSimpleGaussian(0.5, 0.01, 10, 3) - 0.5) < tol);

    REQUIRE(abs(testGradientDescentSimpleGaussian(0.45, 0.1, 1, 1) - 0.5) < tol);
    REQUIRE(abs(testGradientDescentSimpleGaussian(0.45, 0.1, 1, 3) - 0.5) < tol);
    REQUIRE(abs(testGradientDescentSimpleGaussian(0.45, 0.01, 10, 1) - 0.5) < tol);
    REQUIRE(abs(testGradientDescentSimpleGaussian(0.45, 0.01, 10, 3) - 0.5) < tol);

    REQUIRE(abs(testGradientDescentSimpleGaussian(0.55, 0.1, 1, 1) - 0.5) < tol);
    REQUIRE(abs(testGradientDescentSimpleGaussian(0.55, 0.1, 1, 3) - 0.5) < tol);
    REQUIRE(abs(testGradientDescentSimpleGaussian(0.55, 0.01, 10, 1) - 0.5) < tol);
    REQUIRE(abs(testGradientDescentSimpleGaussian(0.55, 0.01, 10, 3) - 0.5) < tol);
}


