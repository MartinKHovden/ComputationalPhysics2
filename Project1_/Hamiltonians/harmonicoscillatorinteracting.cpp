#include "harmonicoscillatorinteracting.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"
#include "math.h"

using std::cout;
using std::endl;

HarmonicOscillatorInteracting::HarmonicOscillatorInteracting(System* system, double omega) :
        Hamiltonian(system) {
    assert(omega > 0);
    m_omega  = omega;
}

double HarmonicOscillatorInteracting::computeLocalEnergy(std::vector<Particle*> particles) {
    /* Function for calculating the local energy of the system.
     *
     * Params:
     * -------
     * std::vector<Particle*> particles: vector containing the particles in the system. 
     * 
     * Returns:
     * --------
     * double: returns the local energy
     */

    double mass = 1;
    double gamma = m_system->getGamma();

    int num_particles = m_system->getNumberOfParticles();
    int num_dimensions = m_system->getNumberOfDimensions();

    //Calculating the kinetic energy:
    double kinetic_energy = m_system->getWaveFunction()->computeDoubleDerivative(particles);
    kinetic_energy *= -0.5;

    //Calculating the potential energy:
    double term2 = 0;

    for(int i = 0; i < num_particles; i++)
    {
        std::vector<double> particle_position = particles[i]->getPosition();

        for(int j = 0; j < num_dimensions; j++)
        {

            if(j == 2)
            {
                term2 += gamma*gamma*particle_position[j]*particle_position[j];
            }
            else
            {
                term2 += particle_position[j]*particle_position[j];
            }
        }
    }

    term2 = 0.5*term2;

    //Calculating the last term of the Hamiltonianm which can be set to zero, since we will never go to a position where anything else can happen:
    double term3=0;

    double local_energy = kinetic_energy + term2 + term3; 

    return local_energy;
} 

// double HarmonicOscillatorInteracting::computeTotalPotentialInt(std::vector<Particle*> particles)
// {
//     int num_particles = m_system->getNumberOfParticles();
//     int num_dims = m_system->getNumberOfDimensions();

//     double interactingPotential = 0;
//     double a = m_system->getA(); 

//     for(int i = 0; i < num_particles; i++)
//     {
//         for(int j = i+1; j < num_particles; j++)
//         {
//             interactingPotential += computePotentialInt(particles, i, j, a);
//         }
//     }
// }

// double HarmonicOscillatorInteracting::computePotentialInt(std::vector<Particle*> particles, int i , int j, double a)
// {
//     std::vector<double> position_particle1 = particles[i]->getPosition();
//     std::vector<double> position_particle2 = particles[j]->getPosition();
//     double return_value;
//     if(computeDistance(position_particle1, position_particle2) <= a)
//     {
//         return_value = 10000000000000000;
//     }
//     else 
//     {
//         return_value = 0;
//     }
// }

double HarmonicOscillatorInteracting::computeDistance(std::vector<double> p1, std::vector<double> p2)
{
    /* Function for computing the distance between two particles given their position. 
     * 
     * Params:
     * ------
     * std::vector<double> p1: vector with particle coordinates.
     * 
     * Returns:
     * -------
     * double: returns the distance between the two particles. 
     */
    double distance = 0;
    for(int i = 0; i < m_system->getNumberOfDimensions(); i++)
    {
        distance += (p1[i] - p2[i])*(p1[i] - p2[i]);
    }
    return sqrt(distance);
}




