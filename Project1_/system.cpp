#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <random>
#include <fstream> 
#include <time.h>
#include "InitialStates/randomuniform.h"

using namespace std;

std::random_device rd;
std::mt19937_64 gen(rd());
std::normal_distribution<double> Normaldistribution(0.0, 1.0);
std::uniform_real_distribution<double> UniformGenerator(0.0, 1.0);

bool System::metropolisStep() {
    /* Function to perform on step in the Metropolis algorithm. Proposes a new 
     * position, calculates the probability ratio and then accepts/rejects the 
     * proposed step. The new position is propsed for one random particle. 
     * 
     * The function updates the system variables while running. 
     * 
     * Returns:
     * A boolean variable that indicates wheter the new step is accepted or rejected.
     */

    // Chooses a random particle number to update the position. 
    int random_particle_number       = rand() % m_particles.size(); 
    
    class Particle random_particle   = *m_particles[random_particle_number];

    vector<double> old_position      = random_particle.getPosition();

    class WaveFunction *temp         = getWaveFunction();

    // Evaluates the current wavefunction value. 
    double wave_function_old         =  temp->evaluate(m_particles);

    vector<double> new_position;

    // Proposes a new position for the randommly chosen particle. 
    for(int i = 0; i < old_position.size(); i++)
    {
        int temp_ran = rand() % 10;

        if(temp_ran < 4.5)
        {
            new_position.push_back(old_position[i] - 1*getStepLength());
        }
        else
        {
            new_position.push_back(old_position[i] + 1*getStepLength());
        }
    }

    // Sets the chosen particles position to the proposed position and calculates the new wave function value. 
    m_particles[random_particle_number]->setPosition(new_position);

    double wave_function_new = temp->evaluate(m_particles);

    double prob_limit = (double) rand() / RAND_MAX;

    // Does the Metropolis test by calculating the ratio between the old and new wave function value. 
    if(prob_limit < (wave_function_new*wave_function_new)/(wave_function_old*wave_function_old))
    {
        // If accepted, keeps the new position. 
        return true;
    }
    else
    {
        // If rejected, sets the position back to the old one. 
        m_particles[random_particle_number] -> setPosition(old_position);
        return false;
    }

}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    /* Function for running the Metropolis algorithm. Calls the metropolisStep function
     * multiple times and calculates the final Monte Carlo estimation of the ground state 
     * energy.
     * 
     * Params:
     * numberOfMetropolisSteps: int, Number of iterations in the Metropolis algorithm. 
     */

    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    ofstream outfile;
    outfile.open("Data/interacting_metropolis_local_energy_values_a_" + to_string(getA()) + 
                        "_num_particles_" + to_string(getNumberOfParticles()) + "_num_dims_" + 
                                        to_string(getNumberOfDimensions()) + "_alpha_" + 
                                                to_string(getWaveFunction()->getParameters()[0]) + ".txt"  );

  
    clock_t start_time = clock(); 

    // Sets up a vector to record the local energies for each iteration of the Metropolis algo. 
    std::vector<double> samples;
    samples.reserve(numberOfMetropolisSteps);

    // Does the Metropolis steps "numberOfMetropolisSteps" times and samples the energy for each step. 
    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();

        m_sampler->sample(acceptedStep);

        samples[i] = m_sampler->getLocalEnergy();

        if( i % 100000 ==0)
        {
            cout << i << endl;
        }

    }

    clock_t end_time = clock();

    // Writes the recorded samples to file. 
    for(int i= 0; i < numberOfMetropolisSteps; i++)
    {
        outfile << samples[i] << endl;
    }


    // Computes the values and writes to terminal. 
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();

    cout << "Time elapsed: " << ((float)(end_time - start_time))/CLOCKS_PER_SEC << " sec" << endl;

    outfile << "Time: " << ((float)(end_time - start_time))/CLOCKS_PER_SEC << endl;

    outfile.close();

}

bool System::importanceSamplingStep(double timestep)
{
    /* Function to perform on step in the Metropolis Importance Sampling algorithm. 
     * Proposes a new position, calculates the probability ratio and then accepts/rejects the 
     * proposed step. The new position is propsed for one random particle. 
     * 
     * The function updates the system variables while running. 
     * 
     * Params:
     * double timestep: The timestep in the position update. 
     * 
     * Returns:
     * A boolean variable that indicates wheter the new step is accepted or rejected.
     */

    int random_particle_number       = rand() % m_particles.size(); 
    
    class Particle random_particle   = *m_particles[random_particle_number];

    vector<double> old_position      = random_particle.getPosition();

    class WaveFunction *temp         = getWaveFunction();

    double wave_function_old         = temp->evaluate(m_particles);


    vector<double> new_position;

    //Sets up containers for the drift-forces and the gradients. 
    double prev_drift_force[m_numberOfDimensions]; 
    double new_drift_force[m_numberOfDimensions];
    double prev_gradient[m_numberOfDimensions];
    double new_gradient[m_numberOfDimensions];

    temp->computeDerivative(prev_gradient, m_particles, random_particle_number);
    
    temp->computeDriftForce(prev_drift_force, prev_gradient , random_particle_number);

    for(int i = 0; i < m_numberOfDimensions; i++)
    {
        double random_int_1 = Normaldistribution(gen);

        new_position.push_back(old_position[i]  + random_int_1*sqrt(timestep) + prev_drift_force[i]*m_D*timestep);
    }

    m_particles[random_particle_number]->setPosition(new_position);

    temp->computeDerivative(new_gradient, m_particles, random_particle_number);

    temp->computeDriftForce(new_drift_force, new_gradient, random_particle_number);

    double wave_function_new = temp->evaluate(m_particles);

    double greens_function_argument = 0.0;

    // Calculate the greens function argument of the ratio for this particle. 
    for(int j = 0; j < m_numberOfDimensions; j++)
    {
        greens_function_argument += -(old_position[j] - new_position[j] - m_D*timestep*new_drift_force[j])*(old_position[j] - new_position[j] - m_D*timestep*new_drift_force[j]);
        greens_function_argument -= -(new_position[j] - old_position[j] - m_D*timestep*prev_drift_force[j])*(new_position[j] - old_position[j] - m_D*timestep*prev_drift_force[j]);
    }

    greens_function_argument /= 4*m_D*timestep;

    double greens_function = exp(greens_function_argument);

    double prob_limit = UniformGenerator(gen); 

    if(prob_limit <= greens_function*(wave_function_new*wave_function_new)/(wave_function_old*wave_function_old))
    {
        return true;
    }
    else
    {
        m_particles[random_particle_number]->setPosition(old_position);
        return false;
    }
}

void System::runImportanceSamplingSteps(int numberOfImportanceSamplingSteps, double timestep)
{
    m_particles                 = m_initialState->getParticles();
    // bool temp = importanceSamplingStep(timestep);
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfImportanceSamplingSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfImportanceSamplingSteps);

    ofstream outfile;
    // outfile.open("importance_sampling_local_energy_values_a_" + to_string(getA()) + "num_particles_" + to_string(getNumberOfParticles()) + "num_dims_" + to_string(getNumberOfDimensions()) + "alpha_" + to_string(getWaveFunction()->getParameters()[0]) + ".txt" );
    outfile.open("Data/non_interacting_importance_sampling_metropolis_local_energy_values_a_" + to_string(getA()) + "_num_particles_" + to_string(getNumberOfParticles()) + "_num_dims_" + to_string(getNumberOfDimensions()) + "_alpha_" + to_string(getWaveFunction()->getParameters()[0]) + ".txt"  );
    std::vector<double> samples;
    samples.reserve(numberOfImportanceSamplingSteps);

    time_t start_time = clock();

    for (int i=0; i < numberOfImportanceSamplingSteps; i++) {
        bool acceptedStep = importanceSamplingStep(timestep);

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */
        m_sampler->sample(acceptedStep);
        samples[i] = m_sampler->getLocalEnergy();
        // outfile << m_sampler->getLocalEnergy() << endl;
        if(i % 100000 == 0)
        {
            cout << i << endl;
        }
    }

    time_t end_time = clock();

    for(int i = 0; i < numberOfImportanceSamplingSteps; i++)
    {
        outfile << samples[i] << endl;
    }

    outfile << "Time: " << ((float)(end_time - start_time))/CLOCKS_PER_SEC << endl;

    outfile.close();
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();

    cout << "Time elapsed: " << ((float)(end_time - start_time))/CLOCKS_PER_SEC << " sec" << endl;
}

double System::runGradientDescent(double stepLength, double initialAlphaValue, int numberOfGDSteps, int numMetropolisSteps, double tolerance)
{
    /* Returns the optimal alpha value for the given system. 
     *
     * 
     */
    int max_number_of_steps     = numberOfGDSteps;
    double tol                  = tolerance;
    int numberOfMetropolisSteps = numMetropolisSteps;

    double current_alpha_value  = initialAlphaValue;
    double gradient = 1;

    int step = 0;

    while(step < max_number_of_steps)
    {
        m_particles                 = m_initialState->getParticles();
        m_sampler                   = new Sampler(this);
        m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
        m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    

        for (int i=0; i < numberOfMetropolisSteps; i++) {
            // bool acceptedStep = metropolisStep(); 
            bool acceptedStep = importanceSamplingStep(0.3);
            m_sampler->sample(acceptedStep);
        }

        m_sampler->computeAverages();





        

        gradient = 2*(m_sampler->getEnergyTimesAlphaDerivative() - (m_sampler->getEnergy()*m_sampler->getAlphaDerivative()));
        cout << "EnergyTimesAlphaDerivative: " << m_sampler->getEnergyTimesAlphaDerivative() << endl;
        cout << "Energy: " << m_sampler->getEnergy() << endl;
        cout << "Alpha Derivative" << m_sampler->getAlphaDerivative() << endl;
        cout <<"Step-length: " << stepLength << endl;
        cout << "Gradient: " << gradient << endl; 
        cout << "Grad times step-length :   "<< stepLength*gradient << endl;
        cout << "Current alpha value       :   " << current_alpha_value << endl;

        m_waveFunction->setAlpha(current_alpha_value);

        current_alpha_value = current_alpha_value - stepLength*gradient;

        setInitialState(new RandomUniform(this, m_numberOfDimensions, m_numberOfParticles));

        if(abs(gradient) < tol && step > 10)
        {
            break;
        }

        step++;

    }

    cout << current_alpha_value << endl;
    return current_alpha_value;

}

void System::runComputeOneBodyDensity(int numberOfMetropolisSteps)
{
    m_particles                 = m_initialState->getParticles();
    // m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    // m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    ofstream outfile;
    outfile.open("NonInteractingOneBodyDensity_" + to_string(m_numberOfParticles) + "_2.txt");

  
    clock_t start_time = clock();


    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();
        // bool acceptedStep = importanceSamplingStep(0.3);

        for(int j = 0; j < m_numberOfParticles; j++)
        {
            double norm = 0;
            std::vector<double> particle_position = m_particles[j]->getPosition();

            for(int dim = 0; dim < m_numberOfDimensions; dim++)
            {
                norm += particle_position[dim]*particle_position[dim];
                // cout << norm << endl;
            }

            norm = sqrt(norm);

            outfile << norm << endl;


        }

        if( i % 1000 ==0)
        {
            cout << i << endl;
        }

    }

    clock_t end_time = clock();

    cout << "Time elapsed: " << ((float)(end_time - start_time))/CLOCKS_PER_SEC << " sec" << endl;

    // outfile << "Time: " << ((float)(end_time - start_time))/CLOCKS_PER_SEC << endl;

    outfile.close();
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setA(double a)
{
    m_a = a;
}

void System::setBeta(double beta)
{
    m_beta = beta;
}

void System::setGamma(double gamma)
{ 
    m_gamma = gamma;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}

void System::setDiffusionConstant(double D)
{
    m_D = D;
}


