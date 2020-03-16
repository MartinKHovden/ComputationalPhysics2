#include "ellipticalgaussian.h"
#include <math.h>
#include <stdio.h>
#include "../system.h"
#include "../particle.h"
#include <cassert>
#include <iostream>
#include "wavefunction.h"

using namespace std;

EllipticalGaussian::EllipticalGaussian(class System *system, double alpha): WaveFunction(system, alpha)
{
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}


double EllipticalGaussian::evaluate(std::vector<class Particle*> particles)
{
    double return_value = evaluate_non_interacting_part(particles)*evaluate_correlation_part(particles);
    // cout << "WaveFunctionVal = " << return_value << endl;
    return return_value;
}

double EllipticalGaussian::evaluate_non_interacting_part(std::vector<class Particle*> particles)
{
    double exp_argument = 0;

    double beta = 2.67;

    for(int i = 0; i < particles.size(); i++)
    {
        double temp = 0;

        class Particle particle = *particles[i]; 
        std::vector<double> particle_position = particle.getPosition();

        for(int dim = 0; dim < particle_position.size(); dim ++)
        {
            if(dim == 2)
            {
                temp += beta*particle_position[dim]*particle_position[dim];            
            }
            else 
            {
                temp += particle_position[dim]*particle_position[dim];
            }
        }

        exp_argument += temp;
    }
    double alpha = m_parameters[0];
    double g_product = exp(-alpha*exp_argument);

    return g_product; 

}

double EllipticalGaussian::evaluate_correlation_part(std::vector<class Particle*> particles)
{
    double num_particles = m_system->getNumberOfParticles();
    double a = 0.0;

    double result = 1;

    for(int j = 0; j < num_particles; j++)
    {
        for(int k = j +1; k < num_particles; k++)
        {
            std::vector<double> r1 = particles[j]->getPosition();
            std::vector<double> r2 = particles[k]->getPosition();

            double distance = calculate_distance(r1, r2);

            if(distance < a)
            {
                result *= 0.0;
            }
            else 
            {
                result *= (1 - a/distance);
            }
        }
    }

    return result;
    
}

double EllipticalGaussian::calculate_distance(std::vector<double> r1_position, std::vector<double> r2_position)
{
    double temp = 0;

    double num_dims = m_system->getNumberOfDimensions();

    for(int i = 0; i < num_dims; i++)
    {
        temp += (r1_position[i] - r2_position[i])*(r1_position[i] - r2_position[i]);
    }

    return sqrt(temp);
}

double EllipticalGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    int num_particles = m_system->getNumberOfParticles();

    double double_derivative = 0;

    double wavefunction_value = evaluate(particles);

    double h = 0.00001;
    
    for(int i = 0; i < num_particles; i++)
    {
        class Particle particle = *particles[i];

        vector<double> particle_coordinates = particle.getPosition();

        vector<double> particle_coordinates_plus_h = particle_coordinates;
        vector<double> particle_coordinates_minus_h = particle_coordinates;


        // for(int j = 0; j < particle_coordinates.size(); j++)
        // {
        //     particle_coordinates_minus_h.push_back(particle_coordinates[j]);
        //     particle_coordinates_plus_h.push_back(particle_coordinates[j]);
        // }

        for(int k = 0; k < m_system->getNumberOfDimensions(); k++)
        {
            particle_coordinates_minus_h[k] = (particle_coordinates[k] - h);
            particle_coordinates_plus_h[k] = (particle_coordinates[k] + h);

            // cout << "R: " << particle_coordinates_minus_h[k] << ", " << particle_coordinates_plus_h[k] << endl;

            particle.setPosition(particle_coordinates_plus_h);
            particles[k]->setPosition(particle_coordinates_plus_h);

            double wavefunction_value_plus_h = evaluate(particles);

            particle.setPosition(particle_coordinates_minus_h);
            particles[k]->setPosition(particle_coordinates_minus_h);
            double wavefunction_value_minus_h = evaluate(particles);

            // cout << wavefunction_value << ", " << wavefunction_value_minus_h << ", " << wavefunction_value_plus_h << endl;

            double_derivative += (wavefunction_value_plus_h - 2*wavefunction_value + wavefunction_value_minus_h);

            particle_coordinates_minus_h[k] = particle_coordinates[k];
            particle_coordinates_plus_h[k] = particle_coordinates[k];            

            particle.setPosition(particle_coordinates);
            particles[k]->setPosition(particle_coordinates);

        }

    }

    double_derivative /= (wavefunction_value*h*h);

    return double_derivative;
}

void EllipticalGaussian::computeDerivative(double *derivative, std::vector<class Particle*> particles, int particle_number )
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
        // cout << wavefunction_value << ", " << wavefunction_value_plus_h << endl;
        derivative[i] = (wavefunction_value_plus_h - wavefunction_value)/(h*wavefunction_value);
        // cout << "Derivative:     " << derivative[i] << endl;
    }
}

void EllipticalGaussian::computeDriftForce(double *drift_force, double * gradient, int particle_number)
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

double EllipticalGaussian::computeAlphaDerivative(std::vector<class Particle*> particles)
{
    ;
}
