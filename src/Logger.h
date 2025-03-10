#ifndef __LOGGER_H__
#define __LOGGER_H__

#include <fstream>
#include <vector>
#include "Particle.h"

class Logger {

    private:
        std::ofstream particleFile;
        std::ofstream kineticEnergyFile;

    public:
        Logger(const std::string& keFilename, const std::string& particleFilename = "");
        ~Logger();

        void logParticleData(double time, const std::vector<Particle>& particles);
        void logKineticEnergy(double time, double kineticEnergy);
        void closeFiles();
};

#endif //__LOGGER_H__