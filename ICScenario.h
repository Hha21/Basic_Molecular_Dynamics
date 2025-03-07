#ifndef IC_SCENARIO_H
#define IC_SCENARIO_H

#include <iostream>
#include <array>

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

//INITIALISE SCENARIOS FOR TEST CASES
void initScenarioTest(ICScenario scenario, const double& Lx, const double& Ly, const double& Lz, const double& dt, double& T, double& temp, double& percType1, unsigned int& N);

//void initScenarioRandom(ICScenario scenario, double temp);

#endif // IC_SCENARIO_H