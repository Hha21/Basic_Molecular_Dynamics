#include "ICScenario.h"
#include "Particle.h"
#include <iostream>
#include <array>

void initScenarioTest(ICScenario scenario, const double& Lx, const double& Ly, const double& Lz, const double& dt, double& T, double& temp) {
    switch(scenario) {

        case ICScenario::ONE: {
            std::cout << "Initiliasing ONE Stationary Particle..." << std::endl;

            // IF T PROVIDED LEAVE, ELSE DEFAULT IS T = 1.0;
            T == 0.0 ? T = 1.0 : T = T;
            
            std::array<double, 3> pos = {10.0, 10.0, 10.0};
            std::array<double, 3> vel = {0.0, 0.0, 0.0};
            std::array<double, 3> frc = {1.0, 1.0, 1.0};

            Particle p1(0, 0, 1.0, pos, vel, frc);
            p1.printID();
            

            break;
        }
        case ICScenario::ONE_VEL: {
            std::cout << "Initializing one moving particle..." << std::endl;
            break;
        }
        case ICScenario::TWO: {
            std::cout << "Initializing two bouncing particles..." << std::endl;
            break;
        }
        case ICScenario::TWO_PASS1: {
            std::cout << "Initializing two passing particles..." << std::endl;
            break;
        }
        default:
            std::cerr << "Error: No valid --ic-* option was provided." << std::endl;
            exit(1);

    }
}

// void initScenarioRandom(ICScenario scenario, double temp = 0.0) {
//     switch(scenario) {
//         case ICScenario::RANDOM:
//             std::cout << "Initializing random particle distribution..." << std::endl;
//             break;
//     }
// }