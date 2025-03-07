#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <vector>
#include <iostream>

#include "ICScenario.h"
#include "Particle.h"
#include "Domain.h"
#include "Logger.h"

class Solver {

    private:
        std::vector<Particle> particles;    //Store of Particles
        Domain domain;                      //Domain for BC
        ICScenario scenario;                //Scenario for Init
        Logger logger;                      //Logger for file write

        const double dt;
        const double T;
        double time;

        double temp;
        const double percType1;
        const unsigned int N;

        double KE;                          //Kinetic Energy

    public:
        Solver(double Lx_, double Ly_, double Lz_, 
                double dt_, double T_, double temp_, 
                double percType1_, unsigned int N_, ICScenario scenario_);

        void initParticles();
        
        // DYNAMICS
        void computeForces();
        void computeKE();
        void step();
        void run();


        // HELPER
        void printForces();
};



#endif