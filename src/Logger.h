#ifndef __LOGGER_H__
#define __LOGGER_H__

#include <fstream>
#include <vector>
#include "Particle.h"

/**
 * @class Logger
 * @brief Handles logging of particle data and kinetic energy to output files.
 *
 * Class adds functionality for writing simulation data to .txt files. If a test
 * case is defined particle.txt and kinetic_energy.txt are written to. Else in 
 * random just kinetic_energy.txt is written to.
 */
class Logger {

    private:
        std::ofstream particleFile;
        std::ofstream kineticEnergyFile;

    public:
        /**
        * @brief Constructs a Logger object with specified output file names.
        * @param particleFile Name of the file to log particle data.
        * @param energyFile Name of the file to log kinetic energy data.
        */
        Logger(const std::string& keFilename, const std::string& particleFilename = "");

        /**
        * @brief Desconstructor to close file streams.
        */
        ~Logger();

        /**
         * @brief Writes particle data to output file.
         * @param time Current simulation time.
         * @param particles Reference to the vector of particles.
         */
        void logParticleData(double time, const std::vector<Particle>& particles);

        /**
         * @brief Writes kinetic energy data to output file.
         * @param time Current simulation time.
         * @param kineticEnergy Computed kinetic energy of the system.
         */
        void logKineticEnergy(double time, double kineticEnergy);
};

#endif //__LOGGER_H__