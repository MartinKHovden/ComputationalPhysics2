#include "randomuniform.h"
#include <iostream>
#include <cassert>
// #include "Math/random.h"
#include "../particle.h"
#include "../system.h"
#include <random>

using std::cout;
using std::endl;

std::random_device rd2;
std::mt19937_64 gen2(rd2());
std::normal_distribution<double> Normaldistribution2(0.0, 1.0);

RandomUniform::RandomUniform(System*    system,
                             int        numberOfDimensions,
                             int        numberOfParticles)  :
        InitialState(system) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;

    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    setupInitialState();
}

void RandomUniform::setupInitialState() {
    /* Function for initializing the particle positions.
     */
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (int j=0; j < m_numberOfDimensions; j++) {
            //The particles are placed according to a gaussian function, since this 
            //resembles the wavefunction of the system. 
            position.push_back(Normaldistribution2(gen2));//*m_system->getStepLength());
        }

        m_particles.push_back(new Particle());
        m_particles.at(i)->setNumberOfDimensions(m_numberOfDimensions);
        m_particles.at(i)->setPosition(position);
    } 

    // for(int i=0; i < m_numberOfParticles; i++)
    // {
    //     for(int j = )
    // }
}
