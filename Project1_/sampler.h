#pragma once

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void computeAverages();
    double getEnergy()                         { return m_energy; }
    double getLocalEnergy()                    { return m_localEnergy; }
    double getEnergyTimesAlphaDerivative()     { return m_energyTimesAlphaDerivative; }
    double getAlphaDerivative()                { return m_alphaDerivative; }

private:
    int     m_numberOfMetropolisSteps = 0;
    int     m_stepNumber = 0;
    double  m_energy = 0;
    double  m_alphaDerivative = 0;
    double  m_energyTimesAlphaDerivative = 0;
    double  m_localEnergy = 0;
    double  m_energySquared = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergySquared = 0;
    double  m_cumulativeAlphaDerivative = 0;
    double  m_cumulativeLocalEnergyTimesAlphaDerivative = 0;
    class System* m_system = nullptr;
};
