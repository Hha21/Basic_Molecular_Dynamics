#include "Logger.h"

Logger::Logger(const std::string& keFilename, const std::string& particleFilename) {
    kineticEnergyFile.open(keFilename);
    kineticEnergyFile << "TIME KE\n";

    if (!particleFilename.empty()) {
        particleFile.open(particleFilename);
        particleFile << "TIME ID X Y Z U V W\n";
    }
}

void Logger::logParticleData(double time, const std::vector<Particle>& particles) {
    for (const Particle& p : particles) {
        particleFile << time << " "
                    << p.getID() << " "
                    << p.getPos()[0] << " " << p.getPos()[1] << " " << p.getPos()[2] << " "
                    << p.getVel()[0] << " " << p.getVel()[1] << " " << p.getVel()[2] << "\n";
    }
}

void Logger::logKineticEnergy(double time, double kineticEnergy) {
    
    kineticEnergyFile << time << " " << kineticEnergy << "\n";  
}

void Logger::closeFiles() {
    this->particleFile.close();
    this->kineticEnergyFile.close();
}

Logger::~Logger() {
    Logger::closeFiles();
}