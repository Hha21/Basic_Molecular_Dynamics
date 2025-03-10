#include "Domain.h"

 /**
 * @brief Adds a force component to the particle force.
 * 
 * @param particles Reference to vector of all particles in 
 * simulation.
 * 
 * The function applies solid wall boundary conditions:
 * if (x_i < 0) :       x_i = -x_i
 *                      u_{i,x} = |u_{i,x}|
 * else if (x_i > L_x): x_i = 2L_x - x_i
 *                      u_{i,x} = -|u_{i,x}|
 */
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

/**
 * @brief Gets the domain dimensions.
 * @return A const reference to the domain size {Lz, Ly, Lz}.
 */
const std::array<double, 3>& Domain::getDims() const {
    return this->domainSize;
}