#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillatorInteracting : public Hamiltonian {
public:
    HarmonicOscillatorInteracting(System* system, double omega);
    double computeLocalEnergy(std::vector<Particle*> particles);
    double computeKineticEnergy();
    double computeTotalPotentialInt(std::vector<Particle*> particles);
    double computePotentialInt(std::vector<Particle*> particles, int i , int j, double a);
    double computeDistance(std::vector<double> p1, std::vector<double> p2);




private:
    double m_omega = 0;
};

