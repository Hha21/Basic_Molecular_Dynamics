#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <array>


/**
 * @struct ParticleProp
 * @brief Contains defined constants for two types of particles.
 */

struct ParticleProp {
    /// Mass[Type0, Type1]
    static constexpr double mass[2] = {1.0, 10.0};
    /// epsilon[T0][T0], epsilon[T0][T1], epsilon[T1][T0], epsilon[T1][T1]
    static constexpr double epsilon[2][2] = {
        {3.0, 5.0},                             
        {15.0, 60.0}                 
    };
    ///Store σ^6: sigma6[T0][T0], sigma6[T0][T1], sigma6[T1][T0], sigma6[T1][T1]
    static constexpr double sigma6[2][2] = {
        {1.0, 64.0},                       
        {64.0, 729.0}
    };
    ///Store σ^12: sigma12[T0][T0], sigma12[T0][T1], sigma12[T1][T0], sigma12[T1][T1]
    static constexpr double sigma12[2][2] = {
        {1.0, 4096.0},
        {4096.0, 531441.0}
    };
};

/**
 * @class Particle
 * @brief Contains variables and methods for a single particle.
 * 
 * Particles interact with Lennard-Jones Forces...
 */

class Particle {

private:
    const unsigned int ID;                      ///< Particle ID (0,...,N-1)
    const unsigned int type;                    ///< Particle Type (0,1)
    const double mass;                          ///< See @ref ParticleProp for masses
    const double invmass;                       ///< Pre-computed

    std::array<double, 3> position;             ///< Particle Position {x,y,z}
    std::array<double, 3> velocity;             ///< Particle Velocity {u,v,w}
    std::array<double, 3> force;                ///< Particle Force {Fx,Fy,Fz}

public:
    
    /**
     * @brief Constructor for Particle Class
     * @param ID_       Particle ID - unsigned int
     * @param type_     Particle Type - unsigned int
     * @param pos       Init. Particle Position (opt.) - std::array<double, 3> = {0.0, 0.0, 0.0}
     * @param vel       Init. Particle Velocity (opt.) - std::array<double, 3> = {0.0, 0.0, 0.0}
     * @param frc       Init. Particle Force (opt.) - std::array<double, 3> = {0.0, 0.0, 0.0}
     */
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
        force(frc) {};


    /// GETTERS
    const std::array<double, 3>& getPos() const;                        ///< Get Position in a const instance 
    const std::array<double, 3>& getVel() const;                        ///< Get Velocity in a const instance
    const unsigned int getType() const;                                 ///< Get Type in a const instance
    const double getMass() const;                                       ///< Get Mass in a const instance
    const unsigned int getID() const;                                   ///< get ID in a const instance
    std::array<double, 3>& getPos();                                    ///< get Position in non-const instance
    std::array<double, 3>& getVel();                                    ///< get Velocity in non-const instance
    
    /// SETTERS
    void addForceComp(int k, double value);                             ///< Set force[k] += value
    void setPos(const std::array<double, 3>& newPos);                   ///< Set Position to newPos
    void setVel(const std::array<double, 3>& newVel);                   ///< Set Velocity to newVel
    void resetForce();                                                  ///< Set Force to [0,0,0]

    void updateVel(const double dt);                                    ///< Integrate Velocity with F/m
    void updatePos(const double dt);                                    ///< Integrate Position with V
};

#endif // __PARTICLE_H__