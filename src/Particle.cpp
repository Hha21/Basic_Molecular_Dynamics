#include "Particle.h"

/// Size definitions for @ref ParticleProp
constexpr double ParticleProp::mass[2];
constexpr double ParticleProp::epsilon[2][2];
constexpr double ParticleProp::sigma6[2][2];
constexpr double ParticleProp::sigma12[2][2];


/**
 * @name Getters
 * @brief Retrieve particle private members
 */


/**
 * @brief Gets the particle position (const).
 * @return A reference to the position array {x, y, z}.
 */
const std::array<double, 3>& Particle::getPos() const {
    return this->position;
}

/**
 * @brief Gets the particle velocity (const).
 * @return A reference to the velocity array {x, y, z}.
 */
const std::array<double, 3>& Particle::getVel() const {
    return this->velocity;
}

/**
 * @brief Gets the particle type (const).
 * @return An unsigned int (particle type).
 */
const unsigned int Particle::getType() const {
    return this->type;
}

/**
 * @brief Gets the particle mass (const).
 * @return double (particle mass).
 */
const double Particle::getMass() const {
    return this->mass;
}

/**
 * @brief Gets the particle ID (const).
 * @return An unsigned int (particle ID).
 */
const unsigned int Particle::getID() const {
    return this->ID;
}

/**
 * @brief Gets the particle position.
 * @return A reference to the position array {x, y, z}.
 */
std::array<double, 3>& Particle::getPos() {
    return this->position;
}

/**
 * @brief Gets the particle velocity.
 * @return A reference to the velocity array {u, v, w}.
 */
std::array<double, 3>& Particle::getVel() {
    return this->velocity;
}

/**
 * @name Setters
 * @brief Set particle private members
 */


 /**
 * @brief Adds a force component to the particle force.
 * 
 * @param k Index of the force component (0 = x, 1 = y, 2 = z).
 * @param value Magnitude of the force to be added.
 * 
 * This function updates the internal force array by adding 
 * the specified force to the k-th component.
 *
 * @see Solver::computeForces
 */
void Particle::addForceComp(int k, double value) {
    this->force[k] += value;
}

/**
 * @brief Sets the particle position.
 * @param newPos A reference to std::array<double, 3> = {x, y, z}.
 */
void Particle::setPos(const std::array<double, 3>& newPos) {
    this->position = newPos;
}

/**
 * @brief Sets the particle velocity.
 * @param newVel A reference to std::array<double, 3> = {u, v, w}.
 */
void Particle::setVel(const std::array<double, 3>& newVel) {
    this->velocity = newVel;
}

/**
 * @brief Resets particle force to {0.0,0.0,0.0}.
 */
void Particle::resetForce() {
    constexpr std::array<double, 3> f0 = {0.0, 0.0, 0.0};
    this->force = f0;
}

 /**
 * @brief Integrates Velocity with F/m.
 * 
 * @param dt Size of timestep (const double)
 */
void Particle::updateVel(const double dt) {
    for (int k = 0; k < 3; ++k) {
        this->velocity[k] += dt * this->force[k] * this->invmass;
    }
}

 /**
 * @brief Integrates Position with Velocity.
 * 
 * @param dt Size of timestep (const double)
 */
void Particle::updatePos(const double dt) {
    for (int k = 0; k < 3; ++k) {
        this->position[k] += dt * this->velocity[k];
    }
}

void Particle::printPos() {
    std::cout << "ID:" << this->ID << " Pos {";
    for (int i = 0; i < 3; ++i) {
        std::cout << " " << this->position[i] << " ,";
    }
    std::cout << "}" << std::endl;
}