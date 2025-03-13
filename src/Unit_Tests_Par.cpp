#include "ICScenario.h"
#include "Solver.h"
#include "Particle.h"

#include <cmath>
#include <iostream>
#include <omp.h>

#define BOOST_TEST_MODULE MOLECULAR_DYNAMIC_TEST
#include <boost/test/included/unit_test.hpp>

const double tolerance = 0.01;


BOOST_AUTO_TEST_CASE(test_molecular_dynamics) {
    std::cout << "RUNNING UNIT TESTS [PARALLEL]..." << std::endl;
    omp_set_num_threads(8);
    for (int i = 1; i < 7; ++i) {
        switch(i) {
            case(1): {
                
                const double Lx = 20.0, Ly = 20.0, Lz = 20.0;
                double percType1 = 100.0;
                const double dt = 0.001;
                const double T = 1.0;
                const double temp = 0.0;
                const unsigned int N = 1;

                ICScenario scenario = ICScenario::ONE;

                Solver solveTest1(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
                std::array<double, 3> finalPos = solveTest1.getFinalPosK(0);
                const double KE = solveTest1.getKE();
                
                for (int i = 0; i < 3; ++i) {
                    BOOST_CHECK_CLOSE(finalPos[i], 10.0, tolerance);
                }

                BOOST_CHECK(KE == 0.0);

                break;
            }

            case(2): {
                
                const double Lx = 20.0, Ly = 20.0, Lz = 20.0;
                double percType1 = 100.0;
                const double dt = 0.001;
                const double T = 20.0;
                const double temp = 0.0;
                const unsigned int N = 1;

                ICScenario scenario = ICScenario::ONE_VEL;

                Solver solveTest2(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
                std::array<double, 3> finalPos = solveTest2.getFinalPosK(0);
                const double KE = solveTest2.getKE();
    
                BOOST_CHECK_CLOSE(finalPos[0], 10.0, tolerance);
                BOOST_CHECK_CLOSE(finalPos[1], 10.0, tolerance);
                BOOST_CHECK_CLOSE(finalPos[2], 10.0, tolerance);

                BOOST_CHECK(KE == 0.5 * (5.0 * 5.0 + 2.0 * 2.0 + 1.0 * 1.0));

                break;
            }

            case(3): {

                const double Lx = 20.0, Ly = 20.0, Lz = 20.0;
                const double dt = 0.001;
                const double T = 50.0;
                const double temp = 0.0;
                double percType1 = 100.0;
                const unsigned int N = 2;

                ICScenario scenario = ICScenario::TWO;

                Solver solveTest3(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
                std::array<double, 3> finalPos = solveTest3.getFinalPosK(0);
                std::array<double, 3> finalPos2 = solveTest3.getFinalPosK(1);
                double KE = solveTest3.getKE();

                //CHECK POS Particle_1. NO FORCING IN Y OR Z DIRECTIONS
                BOOST_CHECK_CLOSE(finalPos[0], 8.507, tolerance);
                BOOST_CHECK_EQUAL(finalPos[1], 10.0);
                BOOST_CHECK_EQUAL(finalPos[2], 10.0);

                //CHECK POS Particle_2, SHOULD BE MIRROR OF 1;
                BOOST_CHECK_CLOSE(finalPos2[0], 10.0 + (10.0 - 8.507), tolerance);
                BOOST_CHECK_EQUAL(finalPos2[1], 10.0);
                BOOST_CHECK_EQUAL(finalPos2[2], 10.0);

                //KE approx 0
                BOOST_CHECK_SMALL(KE, 0.001);

                break;
            }

            case(4): {
                const double Lx = 20.0, Ly = 20.0, Lz = 20.0;
                const double dt = 0.001;
                const double T = 50.0;
                const double temp = 0.0;
                double percType1 = 100.0;
                const unsigned int N = 2;

                ICScenario scenario = ICScenario::TWO_PASS1;

                Solver solveTest4(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
                std::array<double, 3> finalPos = solveTest4.getFinalPosK(0);
                std::array<double, 3> finalPos2 = solveTest4.getFinalPosK(1);
                double KE = solveTest4.getKE();

                //CHECK POS Particle_1. NO FORCING IN Z DIRECTION
                BOOST_CHECK_CLOSE(finalPos[0], 7.2485, tolerance);
                BOOST_CHECK_CLOSE(finalPos[1], 5.7155, tolerance);
                BOOST_CHECK_EQUAL(finalPos[2], 10.0);

                //CHECK POS Particle_2. MIRRORS Particle_1
                BOOST_CHECK_CLOSE(finalPos2[0], 10.0 + (10.0 - finalPos[0]), tolerance);
                BOOST_CHECK_CLOSE(finalPos2[1], 10.0 + (10.0 - finalPos[1]), tolerance);
                BOOST_CHECK_EQUAL(finalPos[2], 10.0);

                //CHECK KE
                BOOST_CHECK_CLOSE(KE, 0.24794, tolerance);

                break;
            }

            case(5): {
                const double Lx = 20.0, Ly = 20.0, Lz = 20.0;
                const double dt = 0.001;
                const double T = 50.0;
                const double temp = 0.0;
                double percType1 = 100.0;
                const unsigned int N = 2;

                ICScenario scenario = ICScenario::TWO_PASS2;

                Solver solveTest5(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
                std::array<double, 3> finalPos = solveTest5.getFinalPosK(0);
                std::array<double, 3> finalPos2 = solveTest5.getFinalPosK(1);
                double KE = solveTest5.getKE();

                //CHECK POS Particle_1. NO FORCING IN Z DIRECTION
                BOOST_CHECK_CLOSE(finalPos[0], 12.754, tolerance);
                BOOST_CHECK_CLOSE(finalPos[1], 18.173, tolerance);
                BOOST_CHECK_EQUAL(finalPos[2], 10.0);

                //CHECK POS Particle_2. MIRRORS Particle_1
                BOOST_CHECK_CLOSE(finalPos2[0], 20.0 - finalPos[0], tolerance);
                BOOST_CHECK_CLOSE(finalPos2[1], 20.0 - finalPos[1], tolerance);
                BOOST_CHECK_EQUAL(finalPos[2], 10.0);

                //CHECK KE
                BOOST_CHECK_CLOSE(KE, 0.24694, tolerance);

                break;
            }

            case(6): {
                const double Lx = 20.0, Ly = 20.0, Lz = 20.0;
                const double dt = 0.001;
                const double T = 50.0;
                const double temp = 0.0;
                double percType1 = 100.0;
                const unsigned int N = 2;

                ICScenario scenario = ICScenario::TWO_PASS2;

                Solver solveTest6(Lx, Ly, Lz, dt, T, temp, percType1, N, scenario);
                std::array<double, 3> finalPos = solveTest6.getFinalPosK(0);
                std::array<double, 3> finalPos2 = solveTest6.getFinalPosK(1);
                double KE = solveTest6.getKE();

                //CHECK POS Particle_1. NO FORCING IN Z DIRECTION
                BOOST_CHECK_CLOSE(finalPos[0], 12.754, tolerance);
                BOOST_CHECK_CLOSE(finalPos[1], 18.173, tolerance);
                BOOST_CHECK_EQUAL(finalPos[2], 10.0);

                //CHECK POS Particle_2. CHECK CIRCULAR ORBIT
                BOOST_CHECK_CLOSE(std::abs(10.0 - finalPos2[0]), std::abs(10.0 - finalPos[0]), tolerance);
                BOOST_CHECK_CLOSE(std::abs(10.0 - finalPos2[1]), std::abs(10.0 - finalPos[1]), tolerance);
                BOOST_CHECK_EQUAL(finalPos[2], 10.0);

                //CHECK KE
                BOOST_CHECK_CLOSE(KE, 0.24694, tolerance);
            }
        }
    }
    std::cout << "UNIT TESTS COMPLETED" << std::endl;
}
