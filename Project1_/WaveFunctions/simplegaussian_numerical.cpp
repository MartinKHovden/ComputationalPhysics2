#include "simplegaussian_numerical.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include "stdio.h"
#include <iostream>
 
using namespace std;

SimpleGaussianNumerical::SimpleGaussianNumerical(class System* system, double alpha) :
        WaveFunction(system, alpha) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussianNumerical::evaluate(std::vector<class Particle*> particles) {
    /* Function for evaluating the wave function value. 
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

double SimpleGaussianNumerical::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* Computes the Laplacian of the wavefunction numerically. 
     * This quantity is needed to compute the (local) energy.
     */

    int num_particles = m_system->getNumberOfParticles();

    double double_derivative = 0;

    double wavefunction_value = evaluate(particles);

    double h = 0.001;  //Step-length in the numerical derivative.
    double alpha = m_parameters[0];

    for(int i = 0; i < num_particles; i++)
    {
        class Particle particle = *particles[i];

        vector<double> particle_coordinates = particle.getPosition();

        vector<double> particle_coordinates_plus_h = particle_coordinates;
        vector<double> particle_coordinates_minus_h = particle_coordinates;

        // Computes the numerical double derivative.
        for(int k = 0; k < m_system->getNumberOfDimensions(); k++)
        {
            particle_coordinates_minus_h[k] = (particle_coordinates[k] - h);
            particle_coordinates_plus_h[k] = (particle_coordinates[k] + h);

            particles[i]->setPosition(particle_coordinates_plus_h);

            double wavefunction_value_plus_h = evaluate(particles);

            particles[i]->setPosition(particle_coordinates_minus_h);
            double wavefunction_value_minus_h = evaluate(particles);

            double_derivative += (wavefunction_value_plus_h - 2*wavefunction_value + wavefunction_value_minus_h)/(h*h);

            particle_coordinates_minus_h[k] = particle_coordinates[k];
            particle_coordinates_plus_h[k] = particle_coordinates[k];            

            particles[i]->setPosition(particle_coordinates);
        }
    }

    double_derivative /= (wavefunction_value);

    return double_derivative;
}

void SimpleGaussianNumerical::computeDerivative(double *derivative, std::vector<class Particle*> particles, int particle_number)
{
    /* Computes the derivative of the wavefunction with respect to one of the particles to be used 
     * in importance sampling. 
     */ 

    class Particle particle = *particles[particle_number];
    vector<double> particle_coordinates = particle.getPosition();
    vector<double> particle_coordinates_plus_h = particle_coordinates;
    double h = 0.0001;

    double wavefunction_value = evaluate(particles);


    for(int i = 0; i < particle_coordinates.size(); i++)
    {
        particle_coordinates_plus_h[i] += h;
        particles[particle_number]->setPosition(particle_coordinates_plus_h);
        double wavefunction_value_plus_h = evaluate(particles);
        particles[particle_number]->setPosition(particle_coordinates);
        particle_coordinates_plus_h[i] -= h;
        derivative[i] = (wavefunction_value_plus_h - wavefunction_value)/(h*wavefunction_value);
    }
}

void SimpleGaussianNumerical::computeDriftForce(double *drift_force, double * gradient, int particle_number)
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

double SimpleGaussianNumerical::computeAlphaDerivative(std::vector<Particle*> particles)
{
    ;
}
