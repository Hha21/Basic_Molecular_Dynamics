#ifndef __DOMAIN_H__
#define __DOMAIN_H__

#include <cstdlib>          //std::abs
#include <vector>
#include <array>
#include "Particle.h"

/**
 * @class Domain
 * @brief Represents the simulation domain and applies boundary conditions.
 *
 * The domain class manages the simulation space and ensures that particles
 * respect periodic boundary conditions when they move outside the domain limits.
 */
class Domain {

    private:
        
        const std::array<double, 3> domainSize;         ///< Domain Dimensions {Lx, Ly, Lz}

    public:

        /**
         * @brief Constructor for Domain Class.
         * @param Lx (double) Length of the domain in the x-direction.
         * @param Ly (double) Length of the domain in the y-direction.
         * @param Lz (double) Length of the domain in the z-direction.
         */
        Domain(double Lx_, double Ly_, double Lz_) :
            domainSize{Lx_, Ly_, Lz_} {}

        /**
         * @brief Applies boundary conditions (solid wall) to particles.
         * @param particles Reference to vector of all particles.
         */
        void applyBC(std::vector<Particle>& particles);
        
        ///GETTERS
        const std::array<double, 3>& getDims() const;   ///< Gets Domain Dimensions (const reference)
        
};







#endif //__DOMAIN_H__