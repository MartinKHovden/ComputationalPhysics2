#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include "stdio.h"
#include <iostream>
#include "ellipticalgaussiananalytical.h"


using namespace std;

EllipticalGaussianAnalytical::EllipticalGaussianAnalytical(class System* system, double alpha) : WaveFunction(system, alpha)
{
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double EllipticalGaussianAnalytical::evaluate(std::vector<class Particle*> particles)
{
    /* Function for evaluating the wave function value. 
     */

    double return_value = evaluate_non_interacting_part(particles)*evaluate_correlation_part(particles);

    return return_value;
}

double EllipticalGaussianAnalytical::evaluate_non_interacting_part(std::vector<class Particle*> particles)
{
    /* Function for evaluating the gaussian part of the wave function. 
     */

    double exp_argument = 0;

    double beta = m_system->getBeta();

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

double EllipticalGaussianAnalytical::evaluate_correlation_part(std::vector<class Particle*> particles)
{
    /* Function for evaluating the Jastrow factor (correlation part). 
     */

    double num_particles = m_system->getNumberOfParticles();
    double a = m_system->getA();

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

double EllipticalGaussianAnalytical::calculate_distance(std::vector<double> r1_position, std::vector<double> r2_position)
{
    /* Function for calcualting the distance between two particles. 
     */ 

    double temp = 0;

    double num_dims = m_system->getNumberOfDimensions();

    for(int i = 0; i < num_dims; i++)
    {
        temp += (r1_position[i] - r2_position[i])*(r1_position[i] - r2_position[i]);
    }

    return sqrt(temp);
}


double EllipticalGaussianAnalytical::computeDoubleDerivative(std::vector<class Particle*> particles)
{
    double double_derivative = 0;
    double alpha = m_parameters[0];

    double beta = m_system->getBeta();
    double a = m_system->getA();

    double num_particles = m_system->getNumberOfParticles(); 

    for(int k = 0; k < num_particles; k++)
    {
        std::vector<double> particle_k_position = particles[k]->getPosition();

        //Calculating term 1: laplacian of phi divided by phi
        double term1 = 0;
        int num_dims = m_system->getNumberOfDimensions();

        for(int dim = 0; dim < num_dims; dim++)
        {
            if(dim == 2)
            {
                term1 += 4*alpha*alpha*beta*beta*particle_k_position[dim]*particle_k_position[dim];
                term1 -= 2*alpha*beta;
            }
            else 
            {
                term1 += 4*alpha*alpha*particle_k_position[dim]*particle_k_position[dim];
                term1 -= 2*alpha;
            }
        }

        // std::cout <<  "Temp1 : " << temp1 << std::endl;

        double_derivative += term1;

        //Calculating term 2 and 3

        double term2 = 0;
        double term3 = 0;

        std::vector<double> temp4;
        temp4.resize(num_dims);

        for(int dim = 0; dim < num_dims; dim++)
        {
            temp4[dim] = 0;
        }

        for(int j = 0; j < num_particles; j++)
        {
            std::vector<double> particle_j_position = particles[j]->getPosition();

            double temp3 = 0;

            if(j != k)
            {
                double r_jk = calculate_distance(particle_j_position, particle_k_position);
                double u_derivative_r_jk = computeUDerivative(particle_j_position, particle_k_position, a);

                for (int dim = 0; dim < num_dims; dim++)
                {
                    if(dim == 2)
                    {
                        temp3 += -2*alpha*particle_k_position[dim]*(particle_k_position[dim] - particle_j_position[dim]);

                    }
                    else
                    {
                        temp3 += -2*alpha*beta*particle_k_position[dim]*(particle_k_position[dim] - particle_j_position[dim]);
                    }
                    
                    temp4[dim] += (particle_k_position[dim] - particle_j_position[dim])*u_derivative_r_jk/r_jk;
                }

                temp3 = temp3*computeUDerivative(particle_j_position, particle_k_position, a)/r_jk;

                term2 += temp3;

            }
            else 
            {
                term2 += 0;
            }

            
        }


        for(int dim = 0; dim < num_dims; dim++)
        {
            term3 += temp4[dim]*temp4[dim];
        }

        double_derivative += term3;

        // std::cout <<  "Temp2 : " << temp2 << std::endl;

        double_derivative += term2;
   
        //Calculating term 5
        double temp6 = 0;

        for(int j = 0; j < num_particles; j ++)
        {
            std::vector<double> particle_j_position = particles[j]->getPosition();
            if(j != k)
            {
                temp6 += computeUDoubleDerivative(particle_k_position, particle_j_position, a) ;
                temp6 += 2*computeUDerivative(particle_k_position, particle_j_position, a)/calculate_distance(particle_k_position, particle_j_position);
                
            }
            else
            {
                temp6 += 0;
            } 
 
        }

        double_derivative += temp6;

    }    

    return double_derivative;
}

void EllipticalGaussianAnalytical::computeDerivative(double *derivative, std::vector<class Particle*> particles, int particle_number)
{
    /* Computes the derivative of the wavefunction with respect to particle k. 
     */

    std::vector<double> particle_k_coordinates = particles[particle_number]->getPosition();
    
    double alpha = m_parameters[0];
    int num_dims = m_system->getNumberOfDimensions();
    int num_particles = m_system->getNumberOfParticles();
    double a = m_system->getA();

    for(int dim = 0; dim < num_dims; dim ++)
    {
        derivative[dim] = -2*alpha*particle_k_coordinates[dim];

        for(int m = 0; m < num_particles; m++)
        {

            std::vector<double> particle_m_coordinates = particles[m]->getPosition();

            if(m != particle_number)
            {
                double distance_km = calculate_distance(particle_k_coordinates, particle_m_coordinates);
                derivative[dim] += computeUDerivative(particle_k_coordinates, particle_m_coordinates, a)*(particle_k_coordinates[dim] - particle_m_coordinates[dim])/distance_km;
            }
        }
    }

}

void EllipticalGaussianAnalytical::computeDriftForce(double *drift_force, double * gradient, int particle_number)
{ 
    int number_of_dimensions = m_system->getNumberOfDimensions();

    for(int j = 0; j < number_of_dimensions; j ++)
    {
        drift_force[j] = 2*gradient[j];
    }
}


double EllipticalGaussianAnalytical::computeUDerivative(std::vector<double> r1_position, std::vector<double> r2_position, double a)
{
    return a/(calculate_distance(r1_position, r2_position)*(calculate_distance(r1_position, r2_position)) - a);
}

double EllipticalGaussianAnalytical::computeUDoubleDerivative(std::vector<double> r1_position, std::vector<double> r2_position, double a)
{
    return (a*a - 2*a*calculate_distance(r1_position, r2_position))/((calculate_distance(r1_position, r2_position) - a*calculate_distance(r1_position, r2_position))*(calculate_distance(r1_position, r2_position) - a*calculate_distance(r1_position, r2_position)));
}

double EllipticalGaussianAnalytical::computeAlphaDerivative(std::vector<Particle*> particles)
{
    int num_particles = m_system->getNumberOfParticles();
    int num_dims = m_system->getNumberOfDimensions();

    double beta = m_system->getBeta();

    double derivative = 0;

    for(int i = 0; i < num_particles; i++)
    {

        std::vector<double> particle_coordinates = particles[i]->getPosition();

        for(int dim = 0; dim < num_dims; dim++)
        {
            if(dim == 2)
            {
                derivative += beta*particle_coordinates[dim]*particle_coordinates[dim];
            }
            else 
            {
                derivative += particle_coordinates[dim]*particle_coordinates[dim];
            }
        }
    }

    return -derivative;
}