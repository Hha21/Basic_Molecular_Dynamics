/**
 * @file ICScenario.h
 * @brief Parses initial condition scenarios for the simulation.
 */

#ifndef IC_SCENARIO_H
#define IC_SCENARIO_H

#include <iostream>
#include <array>

/**
 * @enum ICScenario
 * @brief Enumeration for different initial condition scenarios.
 * 
 * The enum is used to specify the type of initial condition setup.
 */
enum class ICScenario {
    NONE, 
    ONE,
    ONE_VEL,
    TWO,
    TWO_PASS1,
    TWO_PASS2,
    TWO_PASS3,
    RANDOM
};

/**
 * @brief Initializes a predefined test scenario.
 * 
 * This function sets up the initial particle configuration based on a given test case.
 * 
 * @param scenario The chosen ICScenario.
 * @param Lx Domain length in the x-direction.
 * @param Ly Domain length in the y-direction.
 * @param Lz Domain length in the z-direction.
 * @param dt Time step size.
 * @param T Total simulation time.
 * @param temp Initial system temperature.
 * @param percType1 Percentage of type-1 particles.
 * @param N Number of particles.
 */
void initScenario(ICScenario scenario, const double& Lx, const double& Ly, const double& Lz, const double& dt, double& T, double& temp, double& percType1, unsigned int& N);

//void initScenarioRandom(ICScenario scenario, double temp);

#endif // IC_SCENARIO_H