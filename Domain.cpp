#include "Domain.h"

void Domain::applyBC(std::vector<Particle>& particles) {

    for (Particle& p : particles) {
        std::array<double, 3>& pos = p.getPos();
        std::array<double, 3>& vel = p.getVel();
        
        for (int i = 0; i < 3; ++i) {
            if (pos[i] <= 0.0) {
                pos[i] = -pos[i];
                vel[i] = std::abs(vel[i]);
            } 
            else if (pos[i] >= this->domainSize[i]) {
                pos[i] = 2.0 * this->domainSize[i] - pos[i];
                vel[i] = - std::abs(vel[i]);
            }
        }
    }
}

const std::array<double, 3>& Domain::getDims() const {
    return this->domainSize;
}