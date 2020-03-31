#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include "stdio.h"
#include <iostream>

using namespace std;

SimpleGaussian::SimpleGaussian(class System* system, double alpha) :
        WaveFunction(system, alpha) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    /* Function for evaluating the wave function value for the particles in the system.
     * 
     * Params:
     * ------
     * particles: vector that contains each particle in the system.
     * 
     * Returns:
     * --------
     * double g_product: The wave function value for the given system.  
     */

    double exp_argument = 0;

    // Computes the argument of the exponential in the wave function. 
    for(int i = 0; i < particles.size(); i++)
    {
        double temp = 0;

        class Particle particle = *particles[i]; 
        vector<double> particle_position = particle.getPosition();

        for(int dim = 0; dim < particle_position.size(); dim ++)
        {
            temp += particle_position[dim]*particle_position[dim];
        }

        exp_argument += temp;
    }
    double alpha = m_parameters[0];
    double g_product = exp(-alpha*exp_argument);
    return g_product;
}


double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* Computes the double derivative of the wavefunction divided by the wavefunction value. 
     * This quantity is needed to compute the energy.
     * 
     * Params:
     * -------
     * particles: Vector that contains the particles in the system. 
     * 
     * Returns: 
     * Returns the laplacian of the wave function divided by the wavefunction value. 
     */

    int num_part = m_system->getNumberOfParticles();
    int num_dims = m_system->getNumberOfDimensions();
    double alpha = m_parameters[0];
    
    double kinetic_energy = -2*alpha*num_part*num_dims;

    double square_sum = 0;

    //Loops over all the particles and computed the laplacian. 
    for(int i = 0; i < num_part; i++)
    {
        double temp = 0;

        class Particle current_particle = *particles[i];

        std::vector<double> current_particle_position = current_particle.getPosition();

        for(int j = 0; j < num_dims; j++)
        {
            square_sum += current_particle_position[j]*current_particle_position[j];
        }
    } 

    double mass = 1;

    kinetic_energy += 4*alpha*alpha*square_sum;

    return kinetic_energy;
}

void SimpleGaussian::computeDerivative(double *derivative, std::vector<class Particle*> particles, int particle_number)
{
    /* Computes the derivative of the wavefunction with respect to one of the particles divided by the wavefunction value to be used 
     * in importance sampling. 
     * 
     * Params:
     * -------
     * *derivative: pointer to a vector that will contain the gradient of the wavefunction. 
     * particles: vector of particles in system. 
     * particle_number: Computes the gradient with respect to this particle. 
     */

    class Particle particle = *particles[particle_number];
    vector<double> particle_coordinates = particle.getPosition();

    int number_of_particles = particles.size();
    int number_of_dimensions = particle_coordinates.size();

    double alpha = m_parameters[0];

    for(int i = 0; i < number_of_dimensions; i ++)
    {
        derivative[i] = -2*alpha*particle_coordinates[i];
    }
}

void SimpleGaussian::computeDriftForce(double *drift_force, double * gradient, int particle_number)
{
    /* Computes the driftforce of the wavefunction with respect to one particle to be used in importance
     * sampling. 
     * 
     * Params:
     * -------
     * drift_force: pointer to a vector that will contain the drift force. 
     * gradient: pointer to a vector with the gradient in it. 
     * particle_number: Computes the drift force with respect to this particle. 
     */ 

    int number_of_dimensions = m_system->getNumberOfDimensions();

    for(int j = 0; j < number_of_dimensions; j ++)
    {
        drift_force[j] = 2*gradient[j];
    }
}
 
double SimpleGaussian::computeAlphaDerivative(std::vector<class Particle*> particles)
{
    /* Computes the derivative of the wave function with respect to alpha.
     * 
     * Params:
     * -------
     * particles: vector of particles in the system. 
     * 
     * Returns:
     * --------
     * double derivative: The derivative of the wavefunction with respect to alpha divided by the wavefunction value. 
     */

    int num_particles = m_system->getNumberOfParticles();
    int num_dims = m_system->getNumberOfDimensions();

    double beta = m_system->getBeta();

    double derivative = 0;

    // Loops over all the particles and adds its contribution. 
    for(int i = 0; i < num_particles; i++)
    {

        std::vector<double> particle_coordinates = particles[i]->getPosition();

        for(int dim = 0; dim < num_dims; dim++)
        {
            if(dim == 2)
            {
                derivative += particle_coordinates[dim]*particle_coordinates[dim];
            }
            else 
            {
                derivative += particle_coordinates[dim]*particle_coordinates[dim];
            }
        }
    }

    return -derivative;

}
