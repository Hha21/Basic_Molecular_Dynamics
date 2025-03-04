#include "Particle.h"
#include <array>

void Particle::printID() {
    for (int i = 0; i < 3; ++i) {
        std::cout << "Pos[" << i << "] = " << this->position[i] << std::endl;
        std::cout << "Frc[" << i << "] = " << this->force[i] << std::endl;
    }
}

void Particle::resetForce() {
    constexpr std::array<double, 3> f0 = {0.0, 0.0, 0.0};
    this->force = f0;
}