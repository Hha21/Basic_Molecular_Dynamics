#include <boost/program_options.hpp>
#include <iostream>
#include <array>

#include "ICScenario.h"
// #include "particle.h"

namespace po = boost::program_options;

int main(int argc, char** argv) {

    po::options_description opts("Allowed Options: ");
    opts.add_options()
        ("help",    "Print available options.")
        ("Lx", po::value<double>()->default_value(20.0), 
            "x length (Angstroms)")
        ("Ly", po::value<double>()->default_value(20.0), 
            "y length (Angstroms)")
        ("Lz", po::value<double>()->default_value(20.0), 
            "z length (Angstroms)")
        ("dt", po::value<double>()->default_value(0.001), 
            "Time-step")
        ("T", po::value<double>(), 
            "Final time")

        //TEST-CASES
        ("ic-one",        "Initial condition: one stationary particle")
        ("ic-one-vel",    "Initial condition: one moving particle")
        ("ic-two",        "Initial condition: two bouncing particles")
        ("ic-two-pass1",  "Initial condition: two passing particles")
        ("ic-two-pass2",  "Initial condition: two passing particles close")
        ("ic-two-pass3",  "Initial condition: two passing particles close, heavy")

        //RANDOM
        ("ic-random",     "Initial condition: N random particles")

        //OTHER PARAMETERS
        ("percent-type1", po::value<double>()->default_value(10.0), 
            "Percentage of type 1 particles with random IC")
        ("N", po::value<int>(), 
            "Number of particles to spawn with random IC")
        ("temp", po::value<double>(), 
            "Temperature (degrees Kelvin)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << opts << std::endl;
        return 0;
    }

    //INIT VARIABLES:
    //  --ESSENTIAL ALL
    const double Lx = vm["Lx"].as<double>();
    const double Ly = vm["Ly"].as<double>();
    const double Lz = vm["Lz"].as<double>();
    const double dt = vm["dt"].as<double>();
    //  --OPTIONAL ALL
    double T = vm.count("T") ? vm["T"].as<double>() : 0.0;
    double temp = vm.count("temp") ? vm["temp"].as<double>() : 0.0;

    //DETERMINE SCENARIO:
    ICScenario scenario = ICScenario::NONE;

    //PRE-DEFINED TEST-CASES:
    if (vm.count("ic-one")) {
        scenario = ICScenario::ONE;
        initScenarioTest(scenario, Lx, Ly, Lz, dt, T, temp);
    } 
    if (vm.count("ic-one-vel")) scenario = ICScenario::ONE_VEL;
    if (vm.count("ic-two")) scenario = ICScenario::TWO;
    if (vm.count("ic-two-pass1")) scenario = ICScenario::TWO_PASS1;
    if (vm.count("ic-two-pass2")) scenario = ICScenario::TWO_PASS2;
    if (vm.count("ic-two-pass3")) scenario = ICScenario::TWO_PASS3;

    //RANDOM TEST:
    if (vm.count("ic-random")) {
        if (!vm.count("percent-type1") || !vm.count("N") || !vm.count("T")) {
            std::cerr << "Error: --ic-random called but insufficient arguments provided" << std::endl;
            return 1;
        } else {
            double percType1 = vm["percent-type1"].as<double>();
            const int N = vm["N"].as<int>();
            const double T = vm["T"].as<double>();
            scenario = ICScenario::RANDOM;
        }
    }
    //ELSE THROW ERROR:
    if (scenario == ICScenario::NONE) {
        std::cerr << "Error: No --ic-* option was provided. Use --help for available options.\n";
        return 1;
    }

    return 0;
}