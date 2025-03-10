#include "Logger.h"

/**
 * @brief Constructs a Logger object, opening file streams for writing.
 * @param particleFile Name of the file to store particle data.
 * @param energyFile Name of the file to store kinetic energy data.
 */
Logger::Logger(const std::string& keFilename, const std::string& particleFilename) {
    kineticEnergyFile.open(keFilename);
    kineticEnergyFile << "TIME KE\n";

    if (!particleFilename.empty()) {
        particleFile.open(particleFilename);
        particleFile << "TIME ID X Y Z U V W\n";
    }
}

/**
 * @brief Logs particle data including time, position, and velocity.
 * @param time Current time in the simulation.
 * @param particles Reference to vector containing all particles.
 */
void Logger::logParticleData(double time, const std::vector<Particle>& particles) {
    for (const Particle& p : particles) {
        particleFile << time << " "
                    << p.getID() << " "
                    << p.getPos()[0] << " " << p.getPos()[1] << " " << p.getPos()[2] << " "
                    << p.getVel()[0] << " " << p.getVel()[1] << " " << p.getVel()[2] << "\n";
    }
}

/**
 * @brief Logs kinetic energy at a given time step.
 * @param time Current time in the simulation.
 * @param kineticEnergy Computed kinetic energy of the system.
 */
void Logger::logKineticEnergy(double time, double kineticEnergy) {
    
    kineticEnergyFile << time << " " << kineticEnergy << "\n";  
}

/**
 * @brief Closes the file streams when the logger is destroyed.
 */
Logger::~Logger() {
    this->particleFile.close();
    this->kineticEnergyFile.close();
}