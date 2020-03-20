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
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
    double exp_argument = 0;

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
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    // int num_particles = particles.size();

    // double term_2 = 0;

    // for(int i = 0; i < num_particles; i++)
    // {
    //     double temp_2 = 0;

    //     class Particle particle = *particles[i];
    //     vector<double> particle_coordinates = particle.getPosition();
    //     int num_dimensions = particle_coordinates.size();

    //     for(int dim = 0; dim < num_dimensions; dim++)
    //     {
    //         temp_2 += particle_coordinates[dim]*particle_coordinates[dim];
    //     }

    //     term_2 += temp_2;
    // }

    // return term_2*evaluate(particles);

    int num_part = m_system->getNumberOfParticles();
    int num_dims = m_system->getNumberOfDimensions();
    double alpha = m_parameters[0];
    
    double kinetic_energy = -2*alpha*num_part*num_dims;

    double square_sum = 0;

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
    /* Computes the derivative of the wavefunction with respect to one of the particles to be used 
     * in importance sampling. 
     */
    class Particle particle = *particles[particle_number];
    vector<double> particle_coordinates = particle.getPosition();

    int number_of_particles = particles.size();
    int number_of_dimensions = particle_coordinates.size();

    // vector<double> derivative = vector<double>();

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
     */ 

    int number_of_dimensions = m_system->getNumberOfDimensions();

    for(int j = 0; j < number_of_dimensions; j ++)
    {
        drift_force[j] = 2*gradient[j];
    }
}

double SimpleGaussian::computeAlphaDerivative(std::vector<class Particle*> particles)
{

}
