#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <array>
#include <iostream>

class Particle {

private:
    const unsigned int ID;
    const unsigned int type;
    const double mass;
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    std::array<double, 3> force;

public:
    
    // CONSTRUCTOR
    Particle(unsigned int ID_, 
        unsigned int type_, 
        double mass_, 
        std::array<double, 3>& pos, 
        std::array<double, 3>& vel, 
        std::array<double, 3>& frc) 
    : ID(ID_), type(type_), mass(mass_), 
    position(pos), velocity(vel), force(frc) {
        std::cout << "Instance " << ID_ << " of Particle created" << std::endl;
    };

    void printID();
};

#endif // __PARTICLE_H__