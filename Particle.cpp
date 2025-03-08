#include "Particle.h"


// GETTERS
// CONST GETTERS
const std::array<double, 3>& Particle::getPos() const {
    return this->position;
}
const std::array<double, 3>& Particle::getVel() const {
    return this->velocity;
}
const unsigned int Particle::getType() const {
    return this->type;
}
const double Particle::getMass() const {
    return this->mass;
}
const unsigned int Particle::getID() const {
    return this->ID;
}

// NON-CONST GETTERS
std::array<double, 3>& Particle::getPos() {
    return this->position;
}

std::array<double, 3>& Particle::getVel() {
    return this->velocity;
}


// SETTERS

void Particle::setPos(const std::array<double, 3>& newPos) {
    this->position = newPos;
}

void Particle::setVel(const std::array<double, 3>& newVel) {
    this->velocity = newVel;
}

void Particle::resetForce() {
    constexpr std::array<double, 3> f0 = {0.0, 0.0, 0.0};
    this->force = f0;
}

void Particle::addForce(const std::array<double, 3>& force, bool inv) {
    if (!inv) {
        for (int k = 0; k < 3; ++k) {
            this->force[k] += force[k];
        }
    } else {
        for (int k = 0; k < 3; ++k) {
            this->force[k] -= force[k];
        }
    }
}

//EULER FOR VELOCITY
void Particle::updateVel(const double dt) {
    for (int k = 0; k < 3; ++k) {
        this->velocity[k] += dt * this->force[k] * this->invmass;
    }
}

void Particle::updatePos(const double dt) {
    for (int k = 0; k < 3; ++k) {
        this->position[k] += dt * this->velocity[k];
    }
}


void Particle::printID() {

    std::cout << "Pos = {";
    for (int i = 0; i < 3; ++i) {
        std::cout << this->position[i] << " ";
    }
    std::cout << "}" << std::endl;
    std::cout << "Vel = {";
    for (int i = 0; i < 3; ++i) {
        std::cout << this->velocity[i] << " ";
    }
    std::cout << "}" << std::endl;
}

void Particle::printForce() {
    std::cout << "Frc_" << this->ID << " = {";
    for (int i = 0; i < 3; ++i) {
        std::cout << this->force[i] << " ";
    }
    std::cout << "}" << std::endl;
}

void Particle::printPosVel() {
    std::cout << "Pos_" << this->ID << " = {";
    for (int i = 0; i < 3; ++i) {
        std::cout << this->position[i] << " ";
    }
    std::cout << "}" << std::endl;
    std::cout << "Vel_" << this->ID << " = {";
    for (int i = 0; i < 3; ++i) {
        std::cout << this->velocity[i] << " ";
    }
}