#include <iostream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep) {
    //Initialize the parameters if it is the start of simulation. 
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_cumulativeEnergySquared = 0;
        m_cumulativeAlphaDerivative = 0;
        m_cumulativeLocalEnergyTimesAlphaDerivative = 0;
    } 

    //Sampling of the interesting quantities. 
    double localEnergy = m_system->getHamiltonian()-> 
                         computeLocalEnergy(m_system->getParticles());

    double localAlphaDerivative = m_system->getWaveFunction()->computeAlphaDerivative(m_system->getParticles());


    m_localEnergy = localEnergy;
    m_cumulativeEnergy  += localEnergy;
    m_cumulativeEnergySquared += localEnergy*localEnergy;
    m_cumulativeAlphaDerivative += localAlphaDerivative;
    m_cumulativeLocalEnergyTimesAlphaDerivative += localEnergy*localAlphaDerivative;
    m_stepNumber++;
}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Reults -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << endl;
}

void Sampler::computeAverages() {
    /* Computes the final estimates of the quantities of interest. 
     */
    m_energy = m_cumulativeEnergy / m_system->getNumberOfMetropolisSteps();
    m_energySquared = m_cumulativeEnergySquared / m_system->getNumberOfMetropolisSteps();
    m_alphaDerivative = m_cumulativeAlphaDerivative / m_system->getNumberOfMetropolisSteps();
    m_energyTimesAlphaDerivative = m_cumulativeLocalEnergyTimesAlphaDerivative / m_system->getNumberOfMetropolisSteps();
}
