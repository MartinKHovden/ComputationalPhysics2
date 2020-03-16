#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

class EllipticalGaussian : public WaveFunction
{
public:
    EllipticalGaussian(class System *system, double alpha);
    double evaluate(std::vector<class Particle*> particles);
    double evaluate_non_interacting_part(std::vector<class Particle*> particles);
    double evaluate_correlation_part(std::vector<class Particle*> particles);
    double calculate_distance(std::vector<double> r1_position, std::vector<double> r2_position);
    double computeDoubleDerivative(std::vector<class Particle*> particles);
    void computeDerivative(double *derivative, std::vector<class Particle*> particles, int particle_number);
    void computeDriftForce(double *drift_force, double * gradient, int particle_number);
    double computeAlphaDerivative(std::vector<class Particle*> particles);

};