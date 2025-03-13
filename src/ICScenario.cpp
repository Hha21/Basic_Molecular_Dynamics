/**
 * @file ICScenario.cpp
 * @brief Implements initialization of predefined test scenarios.
 */

#include "ICScenario.h"
#include "Solver.h"
#include "SolverSerial.h" 

/**
 * @brief Initializes a predefined test scenario.
 * 
 * Sets up the simulation parameters for a chosen scenario, modifying
 * domain dimensions, time step, temperature, and particle count accordingly.
 * Calls @ref Solver to run simulation
 * 
 * @param scenario The chosen initial condition scenario.
 * @param Lx Domain length in the x-direction.
 * @param Ly Domain length in the y-direction.
 * @param Lz Domain length in the z-direction.
 * @param dt Time step size.
 * @param T Total simulation time.
 * @param temp Initial system temperature.
 * @param percType1 Percentage of type-1 particles.
 * @param N Number of particles.
 */
void initScenario(ICScenario scenario, const double& Lx, const double& Ly, const double& Lz, const double& dt, double& T, double& temp, double& percType1, unsigned int& N) {
    switch(scenario) {

        case ICScenario::ONE: {
            // IF PARAMS PROVIDED USE, ELSE USE DEFAULT FOR TEST CASE;
            T == 0.0 ? T = 1.0 : T = T;

            // ENFORCE MANDATORY TEST CASE PARAMS:
            percType1 = 100.0;
            N = 1;
            
            Solver solveTest1(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
            break;
        }
        
        case ICScenario::ONE_VEL: {
            // IF PARAMS PROVIDED USE, ELSE USE DEFAULT FOR TEST CASE;
            T == 0.0 ? T = 20.0 : T = T;

            // ENFORCE MANDATORY TEST CASE PARAMS:
            percType1 = 100.0;
            N = 1;

            Solver solveTest2(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
            break;
        }

        case ICScenario::TWO: {
            // IF PARAMS PROVIDED USE, ELSE USE DEFAULT FOR TEST CASE;
            T == 0.0 ? T = 50.0 : T = T;

            // ENFORCE MANDATORY TEST CASE PARAMS:
            percType1 = 100.0;
            N = 2;

            Solver solveTest3(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
            break;
        }

        case ICScenario::TWO_PASS1: {
            // IF PARAMS PROVIDED USE, ELSE USE DEFAULT FOR TEST CASE;
            T == 0.0 ? T = 50.0 : T = T;

            // ENFORCE MANDATORY TEST CASE PARAMS:
            percType1 = 100.0;
            N = 2;

            Solver solveTest4(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
            break;
        }

        case ICScenario::TWO_PASS2: {
            // IF PARAMS PROVIDED USE, ELSE USE DEFAULT FOR TEST CASE;
            T == 0.0 ? T = 50.0 : T = T;

            // ENFORCE MANDATORY TEST CASE PARAMS:
            percType1 = 100.0;
            N = 2;

            Solver solveTest5(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
            break;
        }

        case ICScenario::TWO_PASS3: {
            // IF PARAMS PROVIDED USE, ELSE USE DEFAULT FOR TEST CASE;
            T == 0.0 ? T = 50.0 : T = T;

            // ENFORCE MANDATORY TEST CASE PARAMS:
            percType1 = 0.0;
            N = 2;

            Solver solveTest6(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
            break;
        }

        case ICScenario::RANDOM: {

            Solver solveRand(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
            break;
        }

        default:
            std::cerr << "Error: No valid --ic-* option was provided." << std::endl;
            exit(1);

    }
}