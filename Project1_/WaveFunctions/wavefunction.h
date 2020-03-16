#pragma once
#include <vector>

class WaveFunction {
public:
    WaveFunction(class System* system, double alpha);
    int     getNumberOfParameters() { return m_numberOfParameters; }
    std::vector<double> getParameters() { return m_parameters; }
    void setAlpha(double alpha)         { m_parameters[0] = alpha; }
    virtual double evaluate(std::vector<class Particle*> particles) = 0;
    virtual double computeDoubleDerivative(std::vector<class Particle*> particles) = 0;
    virtual void computeDerivative(double *derivative, std::vector<class Particle*> particles, int particle_number) = 0;
    virtual void computeDriftForce(double *drift_force, double *gradient, int particle_number) = 0;
    virtual double computeAlphaDerivative(std::vector<class Particle*> particles) = 0;

protected:
    int     m_numberOfParameters = 0;
    std::vector<double> m_parameters = std::vector<double>();
    class System* m_system = nullptr;
};

