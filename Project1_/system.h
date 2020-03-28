#pragma once
#include <vector>

class System {
public:
    bool metropolisStep             ();
    bool importanceSamplingStep     (double timestep);
    void runMetropolisSteps         (int numberOfMetropolisSteps);
    void runImportanceSamplingSteps (int numberOfImportanceSamplingSteps, double timestep);
    void setNumberOfParticles       (int numberOfParticles);
    void setNumberOfDimensions      (int numberOfDimensions);
    void setStepLength              (double stepLength);
    void setA                       (double a);
    void setBeta                    (double beta);
    void setGamma                   (double gamma);
    void setDiffusionConstant       (double D);
    void setEquilibrationFraction   (double equilibrationFraction);
    void setHamiltonian             (class Hamiltonian* hamiltonian); 
    void setWaveFunction            (class WaveFunction* waveFunction);
    void setInitialState            (class InitialState* initialState);
    double runGradientDescent       (double stepLength, double initialAlphaValue, int numberOfGDSteps, int numMetropolisSteps, double tolerance);
    void runComputeOneBodyDensity   (int numberOfMetropolisSteps);
    class WaveFunction*             getWaveFunction()   { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()    { return m_hamiltonian; }
    class Sampler*                  getSampler()        { return m_sampler; }
    std::vector<class Particle*>    getParticles()      { return m_particles; }
    int getNumberOfParticles()          { return m_numberOfParticles; }
    int getNumberOfDimensions()         { return m_numberOfDimensions; }
    double getA()                       { return m_a; }
    double getBeta()                    { return m_beta; }
    double getGamma()                   { return m_gamma; }
    int getNumberOfMetropolisSteps()    { return m_numberOfMetropolisSteps; }
    double getEquilibrationFraction()   { return m_equilibrationFraction; }
    double getStepLength()              { return m_stepLength; }
    double getCurrentWavefunctionValue()     { return m_currentWavefunctionValue; }
    void setCurrentWavefunctionValue(double wfValue)     { m_currentWavefunctionValue = wfValue; }
    std::vector<double> getCurrentPosition() { return m_currentPosition; }
    void setCurrentPosition(std::vector<double> pos)                { m_currentPosition = pos; }

private:
    int                             m_numberOfParticles = 0;
    int                             m_numberOfDimensions = 0;
    int                             m_numberOfMetropolisSteps = 0;
    double                          m_a = 0;
    double                          m_beta = 1; 
    double                          m_gamma = 1;
    double                          m_equilibrationFraction = 0.0;
    double                          m_stepLength = 0.1;
    double                          m_D = 1; 
    double                          m_currentWavefunctionValue = 0;
    std::vector<double>             m_currentPosition = {0,0,0};
    class WaveFunction*             m_waveFunction = nullptr;
    class Hamiltonian*              m_hamiltonian = nullptr;
    class InitialState*             m_initialState = nullptr;
    class Sampler*                  m_sampler = nullptr;
    std::vector<class Particle*>    m_particles = std::vector<class Particle*>();
};

