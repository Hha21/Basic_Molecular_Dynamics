#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <vector>
#include <iostream>
#include <random>

#include "ICScenario.h"
#include "Particle.h"
#include "Domain.h"
#include "Logger.h"


/**
 * @class Solver
 * @brief Runs the Simulation of particle interaction.
 * 
 * This class initializes particles, applies forces, and integrates their motion over time.
 * It also manages boundary conditions and logging.
 */
class Solver {

    private:
        std::vector<Particle> particles;                        ///< Vector Store of Particles
        Domain domain;                                          ///< Simulation Domain for BC
        ICScenario scenario;                                    ///< Initial Condition Scenario
        Logger logger;                                          ///< Logger for output data

        std::mt19937 randGen;                                   ///<Random Number Generator
        std::uniform_real_distribution<double> distribution;    ///< Uniform Distribution [0,1]

        const double dt;                                        ///< Simulation Time-Step
        const double T;                                         ///< Total Simulation Time
        double time;                                            ///< Current Simulation Time

        double temp;                                            ///< Initial Temperature
        const double percType1;                                 ///< % of Type 1 Particles
        const unsigned int N;                                   ///< Number of particles simulated

        double KE;                                              ///< Kinetic Energy

        // GPU MEMORY POINTERS

        double* posX;
        double* posY;
        double* posZ;
        double* FORCE_BUFFER;
        unsigned int* types;


        void initPointers();
        void allocGPUMemory();
        void computeForcesCUDA();
        void cpHostToDevice();
        void cpDeviceToHost();
        void freeGPUMem();

    public:

        /**
         * @brief Constructor for Solver Class.
         * @param Lx_ Domain size in x-direction
         * @param Ly_ Domain size in y-direction
         * @param Lz_ Domain size in z-direction
         * @param dt_ Time step
         * @param T_ Total simulation time
         * @param temp_ Initial temperature
         * @param percType1_ Percentage of type 1 particles
         * @param N_ Number of particles
         * @param scenario_ Initial condition scenario
         */
        Solver(double Lx_, double Ly_, double Lz_, 
                double dt_, double T_, double temp_, 
                double percType1_, unsigned int N_, ICScenario scenario_);

        ~Solver();
        
        /**
        * @brief Initializes particles based on the selected initial condition scenario.
        */
        void initParticles();

        void resetParticles();
        
        /**
        * @brief Scales velocities to set desired intial temperature.
        */
        void setTemp();
        
        /**
        * @brief Computes Kinetic Energy within whole system.
        */
        void computeKE();

        /**
        * @brief Take one simulation step.
        */
        void step();

        /**
        * @brief Run simulation.
        */
        void run();


        double randDouble() {return this->distribution(this->randGen);};    ///< Random double ~[0,1]

        /**
         * @brief Generates a random position within the domain.
         * @return std::vector<double, 3> candidate random position {x,y,z}.
         */
        std::array<double, 3> getRandPos();

        /**
         * @brief Generates a random velocity vector with u_{i} ~ [-0.5,0.5].
         * @return std::vector<double, 3> random velocity {u,v,w}.
         */
        std::array<double, 3> getRandVel();

        /**
        * @brief Checks if an initial position is valid (ensures distance >= 0.5 between new particles).
        * @param pos new position to check.
        * @return true if valid, false otherwise.
        */
        bool isValidPos(const std::array<double, 3>& pos);

        const std::array<double, 3>& getFinalPosK(unsigned int ID);
        const double getKE();
};

// __global__ void LJKernel(double* posX, double* posY, double* posZ, 
//                         unsigned int* types, double* forceX, 
//                         double* forceY, double* forceZ, unsigned int N);

#endif