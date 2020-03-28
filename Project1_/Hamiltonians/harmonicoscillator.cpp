#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Function for computing the local energy.
     *
     * Params:
     * -------
     * std::vector<Particle*> particles: vector containing the particles in the system. 
     * 
     * Returns:
     * --------
     * double localEnergy
     */

    class WaveFunction* temp = m_system->getWaveFunction();

    std::vector<double> system_params = temp->getParameters();

    double alpha = system_params[0];
    double mass = 1;

    int num_particles = m_system->getNumberOfParticles();
    int num_dimensions = m_system->getNumberOfDimensions();

    //Calculating the kinetic energy:
    double kinetic_energy = m_system->getWaveFunction()->computeDoubleDerivative(particles);
    kinetic_energy *= -0.5;

    //Calculating the potential energy:
    double square_sum = 0;

    for(int i = 0; i < num_particles; i++)
    {
        double temp = 0;

        class Particle current_particle = *particles[i];

        std::vector<double> current_particle_position = current_particle.getPosition();

        for(int j = 0; j < num_dimensions; j++)
        {
            square_sum += current_particle_position[j]*current_particle_position[j];
        }
    } 

    double potentialEnergy = 0.5*mass*m_omega*m_omega*square_sum;


    return kinetic_energy + potentialEnergy;
}


