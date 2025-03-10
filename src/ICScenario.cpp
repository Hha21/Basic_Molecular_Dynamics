#include "ICScenario.h"
#include "Solver.h"

//Lx, Ly, Lz, dt, T, temp

void initScenario(ICScenario scenario, const double& Lx, const double& Ly, const double& Lz, const double& dt, double& T, double& temp, double& percType1, unsigned int& N) {
    switch(scenario) {

        case ICScenario::ONE: {
            // IF PARAMS PROVIDED USE, ELSE USE DEFAULT FOR TEST CASE;
            T == 0.0 ? T = 1.0 : T = T;

            // ENFORCE MANDATORY TEST CASE PARAMS:
            percType1 = 100.0;
            N = 1;
            
            Solver solveTest1(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
            // solveTest1.initParticles();

            break;
        }
        
        case ICScenario::ONE_VEL: {
            // IF PARAMS PROVIDED USE, ELSE USE DEFAULT FOR TEST CASE;
            T == 0.0 ? T = 20.0 : T = T;

            // ENFORCE MANDATORY TEST CASE PARAMS:
            percType1 = 100.0;
            N = 1;

            Solver solveTest2(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
            // solveTest2.initParticles();

            break;
        }

        case ICScenario::TWO: {
            // IF PARAMS PROVIDED USE, ELSE USE DEFAULT FOR TEST CASE;
            T == 0.0 ? T = 50.0 : T = T;

            // ENFORCE MANDATORY TEST CASE PARAMS:
            percType1 = 100.0;
            N = 2;

            Solver solveTest3(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
            // solveTest3.initParticles();


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
            percType1 = 50.0;
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

// void initScenarioRandom(ICScenario scenario, double temp = 0.0) {
//     switch(scenario) {
//         case ICScenario::RANDOM:
//             std::cout << "Initializing random particle distribution..." << std::endl;
//             break;
//     }
// }