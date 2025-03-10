#ifndef __DOMAIN_H__
#define __DOMAIN_H__

#include <vector>
#include <array>
#include "Particle.h"

class Domain {

    private:
        //{Lx, Ly, Lz} as 0 lower limit for all.
        const std::array<double, 3> domainSize;

    public:
        Domain(double Lx_, double Ly_, double Lz_) :
            domainSize{Lx_, Ly_, Lz_} {}

        void applyBC(std::vector<Particle>& particles);
        
        //GETTERS
        const std::array<double, 3>& getDims() const;
        
};







#endif //__DOMAIN_H__