#include <iostream>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/simplegaussian_numerical.h"
#include "WaveFunctions/ellipticalgaussian.h"
#include "WaveFunctions/ellipticalgaussiananalytical.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/harmonicoscillatorinteracting.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include <math.h>

using namespace std;

 
int main() {
    int numberOfDimensions  = 3;
    int numberOfParticles   = 10;
    int numberOfSteps       = (int) pow(2.0, 20.0);
    double omega            = 1.0;          // Oscillator frequency.
    double alpha            = 0.6;          // Variational parameter.
    double stepLength       = 0.5;          // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used for equilibration.
    double D                = 0.5;          // Diffusion constant. Will use D = 0.5.
    double timestep         = 0.01;          // Importance sampling time step. 
    double a                = 0.0043;       
    double beta             = 2.828;
    double gamma            = 2.828;

    // INITIALIZE A NEW SYSTEM

    System* system = new System();

    // CHOOSE HAMILTONIAN:

    // system->setHamiltonian              (new HarmonicOscillator(system, omega));
    system->setHamiltonian              (new HarmonicOscillatorInteracting(system, omega));

    // CHOOSE WAVEFUNCTION: 

    // system->setWaveFunction             (new SimpleGaussian(system, alpha));
    // system->setWaveFunction             (new SimpleGaussianNumerical(system, alpha));
    system->setWaveFunction             (new EllipticalGaussianAnalytical(system, alpha));
    // system->setWaveFunction             (new EllipticalGaussian(system, alpha));

    // SET PARAMS:

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setA                        (a);
    system->setBeta                     (beta);
    system->setGamma                    (gamma);
    system->setDiffusionConstant        (D);

    // CHOOSE SIMULATION: 

    system->runMetropolisSteps          (numberOfSteps);
    // system->runImportanceSamplingSteps  (numberOfSteps, timestep);
    // system->runGradientDescent          (0.001, 0.498);

    return 0;
}