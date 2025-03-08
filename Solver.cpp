#include "Solver.h"

//HELPER FUNCTIONS
static inline double pow8(double r2) {
    return r2 * r2 * r2 * r2;
}

static inline double pow14(double r2) {
    return pow8(r2) * r2 * r2 * r2;
}

std::array<double, 3> Solver::getRandPos() {
    //GET {Lx, Ly, Lz} ONCE
    static const std::array<double, 3> dims = this->domain.getDims();

    return  {dims[0] * (Solver::randDouble()),
            dims[1] * (Solver::randDouble()),
            dims[2] * (Solver::randDouble())};
}

std::array<double, 3> Solver::getRandVel() {
    return {-0.5 + (Solver::randDouble()),
            -0.5 + (Solver::randDouble()),
            -0.5 + (Solver::randDouble())};
}

bool Solver::isValidPos(const std::array<double, 3>& pos) {
    // IF NO PARTICLES THEN ADD
    if (this->particles.empty()) return true;

    static const double minDist2 = 0.5 * 0.5;

    for (const Particle& particle : this->particles) {
        const std::array<double, 3> pPos = particle.getPos();
        double dist2 =  (pos[0] - pPos[0]) * (pos[0] - pPos[0]) +
                    (pos[1] - pPos[1]) * (pos[1] - pPos[1]) +
                    (pos[2] - pPos[2]) * (pos[2] - pPos[2]);
        if (dist2 < minDist2) return false;
    }
    return true;
}

// SOLVER METHODS
Solver::Solver(double Lx_, double Ly_, double Lz_, 
                double dt_, double T_, double temp_, 
                double percType1_, unsigned int N_, ICScenario scenario_) 
        :   domain(Lx_, Ly_, Lz_),
            dt(dt_), T(T_), temp(temp_),
            percType1(percType1_), N(N_),
            time(0.0), scenario(scenario_),
            randGen(std::random_device{}()), 
            distribution(0.0, 1.0),
            logger("kinetic_energy.txt",                                    // IF RAND -> NO particles.txt
                scenario == ICScenario::RANDOM ? "" : "particles.txt") {
            
            Solver::initParticles();    
            std::cout << "SOLVER INITIALISED!" << std::endl;
            Solver::computeForces();
            Solver::run();

}

void Solver::initParticles() {
    // IF PARTICLES
    particles.clear();

    //INIT PARTICLES DEPENDING UPON TEST CASE
    switch(this->scenario) {    
        case ICScenario::ONE: {

            std::array<double, 3> p0Pos = {10.0, 10.0, 10.0};
            std::array<double, 3> p0Vel = {0.0, 0.0, 0.0};
            
            particles.emplace_back(0,0,p0Pos,p0Vel);

            break;
        }

        case ICScenario::ONE_VEL: {

            std::array<double, 3> p0Pos = {10.0, 10.0, 10.0};
            std::array<double, 3> p0Vel = {5.0, 2.0, 1.0};

            particles.emplace_back(0,0,p0Pos,p0Vel);

            break;
        }

        case ICScenario::TWO: {

            std::array<double, 3> p0Pos = {8.5, 10.0, 10.0};
            std::array<double, 3> p0Vel = {0.0, 0.0, 0.0};

            std::array<double, 3> p1Pos = {11.5, 10.0, 10.0};
            // VEL SAME AS p0;

            particles.emplace_back(0, 0, p0Pos, p0Vel);
            particles.emplace_back(1, 0, p1Pos, p0Vel);

            break;
        }

        case ICScenario::TWO_PASS1: {

            std::array<double, 3> p0Pos = {8.5, 11.5, 10.0};
            std::array<double, 3> p0Vel = {0.5, 0.0, 0.0};

            std::array<double, 3> p1Pos = {11.5, 8.5, 10.0};
            std::array<double, 3> p1Vel = {-0.5, 0.0, 0.0};

            particles.emplace_back(0, 0, p0Pos, p0Vel);
            particles.emplace_back(1, 0, p1Pos, p1Vel);

            break;
        }

        case ICScenario::TWO_PASS2: {

            std::array<double, 3> p0Pos = {8.5, 11.3, 10.0};
            std::array<double, 3> p0Vel = {0.5, 0.0, 0.0};

            std::array<double, 3> p1Pos = {11.5, 8.7, 10.0};
            std::array<double, 3> p1Vel = {-0.5, 0.0, 0.0};

            particles.emplace_back(0, 0, p0Pos, p0Vel);
            particles.emplace_back(1, 0, p1Pos, p1Vel);

            break;
        }

        case ICScenario::TWO_PASS3: {

            std::array<double, 3> p0Pos = {8.5, 11.3, 10.0};
            std::array<double, 3> p0Vel = {0.5, 0.0, 0.0};

            std::array<double, 3> p1Pos = {11.5, 8.7, 10.0};
            std::array<double, 3> p1Vel = {-0.5, 0.0, 0.0};

            particles.emplace_back(0, 1, p0Pos, p0Vel);
            particles.emplace_back(1, 1, p1Pos, p1Vel);

            break;
        }

        case ICScenario::RANDOM: {
            std::cout << "HERE RANDOM!" << std::endl;

            unsigned int ID_ = 0;

            while (this->particles.size() < this->N) {
                std::array<double, 3> newPos = Solver::getRandPos();
                if (Solver::isValidPos(newPos)) {
                    std::array<double, 3> newVel = Solver::getRandVel();

                    unsigned int type = (rand() / (double)RAND_MAX) < percType1 / 100.0 ? 1 : 0;
                    particles.emplace_back(ID_++, type, newPos, newVel);
                }
            }

            break;
        }

        default:
            std::cerr << "ERR: INVALID SCENARIO!" << std::endl;
            exit(1);
    }
}

void Solver::computeForces() {
    // RESET FORCES EACH TIME CALLED
    for (Particle& p : this->particles) {
        p.resetForce();
    }

    // LOOP FOR UNIQUE PAIRS (i, j) WITH i < j
    for (unsigned int i = 0; i < this->N; ++i) {
        for (unsigned int j = i + 1; j < this->N; ++j) {

            Particle& p1 = this->particles[i];
            Particle& p2 = this->particles[j];

            //PARTICLE TYPES
            unsigned int typ1 = p1.getType();
            unsigned int typ2 = p2.getType();

            std::array<double, 3> r;
            for (int k = 0; k < 3; ++k) {
                r[k] = p2.getPos()[k] - p1.getPos()[k];
            }

            // SQUARED DISTANCE (1 / r^2)
            double inv_r2 = 1.0 / (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        
            // LENNARD-JONES FORCE
            double eps_ij = ParticleProp::epsilon[typ1][typ2];
            double sigma6_ij = ParticleProp::sigma6[typ1][typ2];
            double sigma12_ij = ParticleProp::sigma12[typ1][typ2];

            double force_mag = -24.0 * eps_ij * ((sigma12_ij * pow14(inv_r2)) - (sigma6_ij * pow8(inv_r2)));

            //FORCE VECTOR:
            std::array<double, 3> force;
            for (int k = 0; k < 3; ++k) {
                force[k] = force_mag * r[k];
            }

            //APPLY TO BOTH:
            p1.addForce(force, false);
            p2.addForce(force, true);
        }
    }
}

void Solver::computeKE() {
    //INIT KE
    double KE_ = 0.0;
    for (const Particle& p : this->particles) {
        double speed = 0.0;

        for (int k = 0; k < 3; ++k) {
            speed += p.getVel()[k] * p.getVel()[k];
        }

        KE_ += 0.5 * p.getMass() * speed; 
    } 
    this->KE = KE_;  
}

void Solver::step() {
    // 1 - APPLY BC
    this->domain.applyBC(this->particles);

    // 2 - COMPUTE FORCES (RESET FORCES INSIDE)
    Solver::computeForces();

    // 3 - EULER METHOD FOR UPDATE
    for (Particle& p : particles) {
        p.updateVel(this->dt);
        p.updatePos(this->dt);
    }

    // 4 - UPDATE t
    this->time += this->dt;
}

void Solver::run() {
    std::cout << "RUNNING SIMULATION... " << std::endl;

    double lastOutputTime = -0.1; 

    while (this->time < this->T) { 
        if (time >= lastOutputTime + 0.1) {
            Solver::computeKE();
            logger.logParticleData(this->time, this->particles);
            logger.logKineticEnergy(this->time, this->KE);
            lastOutputTime = this->time;
        }

        Solver::step();
    }
    Solver::computeKE();
    logger.logParticleData(this->time, this->particles);
    logger.logKineticEnergy(this->time, this->KE);
    lastOutputTime = this->time;

    std::cout << "SIMULATION COMPLETE" << std::endl;
}

//TEMP

void Solver::printForces() {
    for (unsigned int i = 0; i < this->N; ++i) {
        particles[i].printForce();
    }
}

void Solver::printPosVels() {
    for (unsigned int i = 0; i < this->N; ++i) {
        particles[i].printPosVel();
    }
}