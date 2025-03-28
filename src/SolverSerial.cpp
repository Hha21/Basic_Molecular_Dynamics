#include "SolverSerial.h"

///Boltzmann Constant k_b
static constexpr double BOLTZMANN = 0.8314459920816467;


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

/**
 * @brief Constructs the solver, intiliases problem and runs upon being called.
 */
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
            std::cout << "SERIAL SOLVER INITIALISED" << std::endl;
            Solver::run();
}

Solver::~Solver(){

}

/**
 * @brief Initializes particles according to the specified initial condition scenario.
 */
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

            unsigned int type1Num = static_cast<unsigned int>(std::ceil((percType1 / 100.0) * this->N));
            unsigned int currType1 = 0;
            unsigned int ID_ = 0;

            while (this->particles.size() < this->N) {
                //GENERATE RANDOM POS THEN CHECK
                std::array<double, 3> newPos = Solver::getRandPos();
                if (Solver::isValidPos(newPos)) {
                    std::array<double, 3> newVel = Solver::getRandVel();
                    unsigned int type = 0;

                    if (currType1 < type1Num) {
                        type = 1;
                        currType1++;
                    }
                    
                    particles.emplace_back(ID_++, type, newPos, newVel);
                }
            }
            break;
        }

        default:
            std::cerr << "ERR: INVALID SCENARIO!" << std::endl;
            exit(1);
    }

    Solver::setTemp();
}


//DYNAMICS 

void Solver::setTemp() {
    if (this->temp <= 0.0) return;

    //GET CURRENT KE
    Solver::computeKE();

    double temp0 = (2.0 * this->KE) / (3.0 * BOLTZMANN * this->N);

    double lambda = std::sqrt(this->temp / temp0);

    for (Particle& p : this->particles) {
        std::array<double, 3> v0 = p.getVel();
        for (int k = 0; k < 3; ++k) v0[k] *= lambda;
        p.setVel(v0);
    }
}


/**
 * @brief Computes the Lennard-Jones force between all particle pairs.
 */
void Solver::computeForces() {
    // RESET FORCES EACH TIME CALLED
    for (Particle& p : this->particles) {
        p.resetForce();
    }

    // LOOP FOR UNIQUE PAIRS (i, j) WITH i < j
    for (unsigned int i = 0; i < this->N; ++i) {

        unsigned int typ1 = this->particles[i].getType();
        const std::array<double, 3>& pos1 = this->particles[i].getPos();
        for (unsigned int j = i + 1; j < this->N; ++j) {

            unsigned int typ2 = this->particles[i].getType();

            const std::array<double, 3>& pos2 = this->particles[j].getPos();
            
            double rx = pos1[0] - pos2[0];
            double ry = pos1[1] - pos2[1];
            double rz = pos1[2] - pos2[2];
            double r2 = rx * rx + ry * ry + rz * rz;

            // SQUARED DISTANCE (1 / r^2)
            double inv_r2 = 1.0 / r2;
            double inv_r4 = inv_r2 * inv_r2;
            double inv_r8 = inv_r4 * inv_r4;
            double inv_r14 = inv_r8 * inv_r4 * inv_r2; 
        
            // LENNARD-JONES FORCE
            double eps_ij = ParticleProp::epsilon[typ1][typ2];
            double sigma6_ij = ParticleProp::sigma6[typ1][typ2];
            double sigma12_ij = ParticleProp::sigma12[typ1][typ2];

            double force_mag = -24.0 * eps_ij * ((2.0 * sigma12_ij * inv_r14) - (sigma6_ij * inv_r8));

            //FORCE VECTOR:
            double fx = force_mag * rx;
            double fy = force_mag * ry;
            double fz = force_mag * rz;
            this->particles[i].addForceComp(0, -fx);
            this->particles[i].addForceComp(1, -fy);
            this->particles[i].addForceComp(2, -fz);
            this->particles[j].addForceComp(0, fx);
            this->particles[j].addForceComp(1, fy);
            this->particles[j].addForceComp(2, fz);
        }
    }
}

void Solver::computeKE() {
    //INIT KE
    double KE_ = 0.0;
    for (const Particle& p : this->particles) {
        double speed = 0.0;
        const std::array<double, 3> pVel = p.getVel();
        for (int k = 0; k < 3; ++k) {
            speed += pVel[k] * pVel[k];
        }

        KE_ += 0.5 * p.getMass() * speed; 
    } 
    this->KE = KE_;  
}

void Solver::step() {

    // 1 - SET TEMP
    Solver::setTemp();

    // 2 - APPLY BC
    this->domain.applyBC(this->particles);

    // 3 - COMPUTE FORCES (RESET FORCES INSIDE)
    Solver::computeForces();

    // 4 - EULER METHOD FOR UPDATE
    for (Particle& p : particles) {
        p.updateVel(this->dt);
        p.updatePos(this->dt);
    }

    // 5 - INCREMENT t
    this->time += this->dt;
}

void Solver::run() {
    std::cout << "RUNNING SIMULATION... " << std::endl;

    double lastOutputTime = -0.1; 

    while (this->time < this->T) { 
        if (this->time >=  lastOutputTime + 0.1) {
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

const std::array<double, 3>& Solver::getFinalPosK(unsigned int ID) {
    //const std::array<double, 3> pos = this->particles[ID].getPos();
    return this->particles[ID].getPos();
}

const double Solver::getKE(){
    return this->KE;
}