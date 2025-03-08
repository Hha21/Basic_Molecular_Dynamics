#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <array>
#include <iostream>

// STRUCT TO ASSIGN PROPERTIES
struct ParticleProp {
    static constexpr double mass[2] = {1.0, 10.0};
    static constexpr double epsilon[2][2] = {
        {3.0, 5.0},                             //epsilon[0][0], epsilon[0][1]
        {15.0, 60.0}                            //epsilon[1][0], epsilon[1][1]
    };
    //STORE SIGMA^6 & SIGMA^12 NOT SIGMA
    static constexpr double sigma6[2][2] = {
        {1.0, 64.0},
        {64.0, 729.0}
    };
    static constexpr double sigma12[2][2] = {
        {1.0, 4096.0},
        {4096.0, 531441.0}
    };
};

class Particle {

private:
    const unsigned int ID;
    const unsigned int type;
    const double mass;
    const double invmass;                       // PRE-COMPUTE;

    std::array<double, 3> position;
    std::array<double, 3> velocity;
    std::array<double, 3> force;

public:
    
    // CONSTRUCTOR
    Particle(unsigned int ID_, 
        unsigned int type_, 
        std::array<double, 3> pos = {0.0, 0.0, 0.0}, 
        std::array<double, 3> vel = {0.0, 0.0, 0.0}, 
        std::array<double, 3> frc = {0.0, 0.0, 0.0}) 
    :   ID(ID_), 
        type(type_), 
        mass(ParticleProp::mass[type_]), 
        invmass(1.0 / (ParticleProp::mass[type_])),
        position(pos), 
        velocity(vel), 
        force(frc) {
            std::cout << "Instance " << ID_ << ", Type " << type << " of Particle created" << std::endl;
    };

    void printID();
    void printForce();
    void printPosVel();

    // GETTERS
    const std::array<double, 3>& getPos() const;
    const std::array<double, 3>& getVel() const;
    const unsigned int getType() const;
    const double getMass() const;
    const unsigned int getID() const;
    std::array<double, 3>& getPos();
    std::array<double, 3>& getVel();
    
    // SETTERS
    void setPos(const std::array<double, 3>& newPos);
    void setVel(const std::array<double, 3>& newVel);
    void resetForce();

    void addForce(const std::array<double, 3>& force, bool inv);
    void updateVel(const double dt);
    void updatePos(const double dt);
};

#endif // __PARTICLE_H__