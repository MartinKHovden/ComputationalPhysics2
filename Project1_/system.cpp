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

using namespace std;

std::random_device rd;
std::mt19937_64 gen(rd());
std::normal_distribution<double> Normaldistribution(0.0, 1.0);
std::uniform_real_distribution<double> UniformGenerator(0.0, 1.0);

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */

    int random_particle_number       = rand() % m_particles.size(); 
    
    class Particle random_particle   = *m_particles[random_particle_number];

    vector<double> old_position      = random_particle.getPosition();

    class WaveFunction *temp         = getWaveFunction();

    double wave_function_old         =  temp->evaluate(m_particles);

    vector<double> new_position;

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


    m_particles[random_particle_number]->setPosition(new_position);

    double wave_function_new = temp->evaluate(m_particles);

    // cout << "new          : " << wave_function_new << "      old: " << wave_function_old << endl;
    // cout << "position     : " << m_particles[random_particle_number]->getPosition()[0] << endl;

    double prob_limit = (double) rand() / RAND_MAX;

    // cout << prob_limit << endl;
    // cout << (wave_function_new*wave_function_new)/(wave_function_old*wave_function_old) << endl;

    if(prob_limit < (wave_function_new*wave_function_new)/(wave_function_old*wave_function_old))
    {
        return true;
    }
    else
    {
        m_particles[random_particle_number] -> setPosition(old_position);
        return false;
    }

}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    ofstream outfile;
    outfile.open("interacting_metropolis_local_energy_values_a_" + to_string(getA()) + "num_particles_" + to_string(getNumberOfParticles()) + "num_dims_" + to_string(getNumberOfDimensions()) + "alpha_" + to_string(getWaveFunction()->getParameters()[0]) + ".txt"  );

  
    clock_t start_time = clock();


    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */
        m_sampler->sample(acceptedStep);
        outfile << m_sampler->getLocalEnergy() << endl;

        if( i % 1000 ==0)
        {
            cout << i << endl;
        }

    }

    clock_t end_time = clock();


    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();

    cout << "Time elapsed: " << ((float)(end_time - start_time))/CLOCKS_PER_SEC << " sec" << endl;

    outfile << "Time: " << ((float)(end_time - start_time))/CLOCKS_PER_SEC << endl;

    outfile.close();

}

bool System::importanceSamplingStep(double timestep)
{
    int random_particle_number       = rand() % m_particles.size(); 
    
    class Particle random_particle   = *m_particles[random_particle_number];

    vector<double> old_position      = random_particle.getPosition();

    class WaveFunction *temp         = getWaveFunction();

    double wave_function_old         = temp->evaluate(m_particles);


    vector<double> new_position;

    double prev_drift_force[m_numberOfParticles]; 
    double new_drift_force[m_numberOfParticles];
    double gradient[m_numberOfParticles];

    temp->computeDerivative(gradient, m_particles, random_particle_number);
    
    temp->computeDriftForce(prev_drift_force, gradient , random_particle_number);

    for(int i = 0; i < m_numberOfDimensions; i++)
    {
        double random_int_1 = Normaldistribution(gen);

        new_position.push_back(old_position[i]  + random_int_1*sqrt(timestep) + prev_drift_force[i]*m_D*timestep);
    }

    m_particles[random_particle_number]->setPosition(new_position);

    temp->computeDriftForce(new_drift_force, gradient, random_particle_number);

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

    double prob_limit = UniformGenerator(gen); //rand()/RAND_MAX;

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
    outfile.open("importance_sampling_local_energy_values_a_" + to_string(getA()) + "num_particles_" + to_string(getNumberOfParticles()) + "num_dims_" + to_string(getNumberOfDimensions()) + "alpha_" + to_string(getWaveFunction()->getParameters()[0]) + ".txt" );


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
        outfile << m_sampler->getLocalEnergy() << endl;
        if(i % 100 == 0)
        {
            cout << i << endl;
        }
    }

    time_t end_time = clock();

    outfile << "Time elapsed: " << ((float)(end_time - start_time))/CLOCKS_PER_SEC << " sec" << endl;

    outfile.close();
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();

    cout << "Time elapsed: " << ((float)(end_time - start_time))/CLOCKS_PER_SEC << " sec" << endl;
}

double System::runGradientDescent(double stepLength, double initialAlphaValue)
{
    /* Returns the optimal alpha value for the system. 
     *
     * 
     */
    double initial_value        = initialAlphaValue;
    int max_number_of_steps     = 10000;
    double tol                  = 0.001;
    int numberOfMetropolisSteps = 1000;

    double current_alpha_value  = initial_value;

    for(int i = 0; i < max_number_of_steps; i++)
    {
        m_particles                 = m_initialState->getParticles();
        m_sampler                   = new Sampler(this);
        m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
        m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    

        for (int i=0; i < numberOfMetropolisSteps; i++) {
            bool acceptedStep = metropolisStep();

            /* Here you should sample the energy (and maybe other things using
            * the m_sampler instance of the Sampler class. Make sure, though,
            * to only begin sampling after you have let the system equilibrate
            * for a while. You may handle this using the fraction of steps which
            * are equilibration steps; m_equilibrationFraction.
            */
            m_sampler->sample(acceptedStep);
        }

        m_sampler->computeAverages();

        

        double gradient = 2*(m_sampler->getEnergyTimesAlphaDerivative() - (m_sampler->getEnergy()*m_sampler->getAlphaDerivative()));
        cout << "Grad times step-length :   "<< stepLength*gradient << endl;
        cout << "Current wf value       :   " << current_alpha_value << endl;

        m_waveFunction->setAlpha(current_alpha_value);

        current_alpha_value = current_alpha_value - stepLength*gradient;

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
    outfile.open("InteractingOneBodyDensity_" + to_string(m_numberOfParticles) + ".txt");

  
    clock_t start_time = clock();


    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();

        for(int i = 0; i < m_numberOfParticles; i++)
        {
            double norm = 0;
            std::vector<double> particle_position = m_particles[i]->getPosition();

            for(int dim = 0; dim < m_numberOfDimensions; dim++)
            {
                norm += particle_position[i]*particle_position[i];
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


