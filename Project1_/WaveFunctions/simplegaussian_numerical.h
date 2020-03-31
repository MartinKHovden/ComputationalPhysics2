#pragma once
#include "wavefunction.h"

class SimpleGaussianNumerical : public WaveFunction {
public:
    SimpleGaussianNumerical(class System* system, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    void computeDerivative(double *derivative, std::vector<class Particle*> particles, int particle_number);
    void computeDriftForce(double *drift_force, double * gradient, int particle_number);
    double computeAlphaDerivative(std::vector<Particle*> particles);

}; 
