#include "SolverSerialCUDA.h"

///Boltzmann Constant k_b
static constexpr double BOLTZMANN = 0.8314459920816467;

// GPU STRUCT FOR PROPERTIES (CANNOT ACCESS STATIC MEMORY)
struct GPU_ParticleProp {
    double mass[2];
    double epsilon[2][2];
    double sigma6[2][2];
    double sigma12[2][2];
};

//DECLARE GPU CONST
__constant__ GPU_ParticleProp GPU_PARTICLE_PROPERTIES;

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
            Solver::allocGPUMemory();
            Solver::initPointers();
            
            std::cout << "CUDA SOLVER INITIALISED" << std::endl;
            Solver::run();
}

void Solver::allocGPUMemory() {

    // TEMP TO COPY PARTICLE PROPERTY VALUES
    GPU_ParticleProp TEMP_PARTICLE_PROPERTIES;

    for (int i = 0; i < 2; i++) {
        TEMP_PARTICLE_PROPERTIES.mass[i] = ParticleProp::mass[i];
        for (int j = 0; j < 2; j++) {
            TEMP_PARTICLE_PROPERTIES.epsilon[i][j] = ParticleProp::epsilon[i][j];
            TEMP_PARTICLE_PROPERTIES.sigma6[i][j] = ParticleProp::sigma6[i][j];
            TEMP_PARTICLE_PROPERTIES.sigma12[i][j] = ParticleProp::sigma12[i][j];
        }
    }

    // COPY PROPERTIES
    cudaMemcpyToSymbol(GPU_PARTICLE_PROPERTIES, &TEMP_PARTICLE_PROPERTIES, sizeof(GPU_ParticleProp));

    // SHARE THESE BETWEEN CPU AND GPU
    cudaMallocManaged(&posX, this->N * sizeof(double));
    cudaMallocManaged(&posY, this->N * sizeof(double));
    cudaMallocManaged(&posZ, this->N * sizeof(double));
    cudaMallocManaged(&FORCE_BUFFER, this->N * 3 * sizeof(double));
    cudaMallocManaged(&types, this->N * sizeof(unsigned int));
}

__device__
void getGPUProperties(const int typ1, const int typ2, double* PROPERTIES) {

    // NOT USED 

    int typSum = typ1 + typ2;
    switch(typSum) {
        case(0): {
            PROPERTIES[0] = 3.0;
            PROPERTIES[1] = 1.0;
            PROPERTIES[2] = 1.0; 
            break;
        }
        case(1): {
            PROPERTIES[0] = 15.0;
            PROPERTIES[1] = 64.0;
            PROPERTIES[2] = 4096.0; 
            break;
        }
        case(2): {
            PROPERTIES[0] = 60.0;
            PROPERTIES[1] = 729.0;
            PROPERTIES[2] = 531441.0; 
            break;
        }
    }
}

__global__ 
void LJKernel(double* posX, double* posY, double* posZ, unsigned int* types,
                            double* FORCE_BUFFER, unsigned int N) {
    
    // COMPUTE FORCES ON i FROM ALL OTHER PARTICLES j

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N)  return;

    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    //double* PROPERTIES = new double[3];

    for (int j = 0; j < N; j++) {

        if (j == i) continue;

        unsigned int typ1 = types[i];
        unsigned int typ2 = types[j];

        double rx = posX[j] - posX[i];
        double ry = posY[j] - posY[i];
        double rz = posZ[j] - posZ[i];
        double r2 = rx * rx + ry * ry + rz * rz;

        double inv_r2 = 1.0 / r2;
        double inv_r4 = inv_r2 * inv_r2;
        double inv_r8 = inv_r4 * inv_r4;
        double inv_r14 = inv_r8 * inv_r4 * inv_r2;

        // getGPUProperties(typ1, typ2, PROPERTIES);

        // double eps_ij = PROPERTIES[0];
        // double sigma6_ij = PROPERTIES[1];
        // double sigma12_ij = PROPERTIES[2];

        double eps_ij = GPU_PARTICLE_PROPERTIES.epsilon[typ1][typ2];
        double sigma6_ij = GPU_PARTICLE_PROPERTIES.sigma6[typ1][typ2];
        double sigma12_ij = GPU_PARTICLE_PROPERTIES.sigma12[typ1][typ2];

        double force_mag = -24.0 * eps_ij * ((2.0 * sigma12_ij * inv_r14) - (sigma6_ij * inv_r8));

        fx += force_mag * rx;
        fy += force_mag * ry;
        fz += force_mag * rz;

        //printf("Value rx, fx, i, j = %f %f %d %d\n", rx, fx, i, j);
    }

    //delete[] PROPERTIES;

    // NOW CPU CAN ACCESS FORCE_BUFFER FOR PARTICLE i UPDATE
    FORCE_BUFFER[i * 3 + 0] = fx;
    FORCE_BUFFER[i * 3 + 1] = fy;
    FORCE_BUFFER[i * 3 + 2] = fz;
}

void Solver::computeForcesCUDA() {
    int THREADS_PER_BLOCK = 256;
    int NUM_INTERACTIONS = (this->N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    int NUM_BLOCKS = (NUM_INTERACTIONS + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

    LJKernel<<<NUM_BLOCKS, THREADS_PER_BLOCK>>>(this->posX, this->posY, this->posZ, this->types, this->FORCE_BUFFER, this->N);

    cudaDeviceSynchronize();
}

void Solver::cpHostToDevice() {

    //DEVICE NEEDS THE NEW POSITIONS FROM EACH PARTICLE

    for (const Particle& p: this->particles) {
        const std::array<double, 3>& pos = p.getPos();
        const unsigned int ID_ = p.getID();
        this->posX[ID_] = pos[0];
        this->posY[ID_] = pos[1];
        this->posZ[ID_] = pos[2];
    }
}

void Solver::cpDeviceToHost() {

    // HOST NEEDS THE NEW FORCES COMPUTED BY DEVICE

    for (unsigned int i = 0; i < this->N; ++i) {
        Particle& p = this->particles[i];
        p.addForceComp(0, this->FORCE_BUFFER[i * 3 + 0]);
        p.addForceComp(1, this->FORCE_BUFFER[i * 3 + 1]);
        p.addForceComp(2, this->FORCE_BUFFER[i * 3 + 2]);
    }
}

void Solver::freeGPUMem() {

    // FREE MEMORY ALLOCATED TO CPU

    cudaDeviceSynchronize();

    cudaFree(posX);
    cudaFree(posY);
    cudaFree(posZ);
    cudaFree(FORCE_BUFFER);
    cudaFree(types);

}


Solver::~Solver(){
    //Solver::freeGPUMem();
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
            Solver::setTemp();
            break;
        }

        default:
            std::cerr << "ERR: INVALID SCENARIO!" << std::endl;
            exit(1);
    }
}


void Solver::resetParticles() {
    for (Particle& p : this->particles) {
        p.resetForce();
    }
}

void Solver::initPointers() {

    // ALLOCATE INITIAL VALUES FOR POS & VALUES FOR TYPES
    for (const Particle& p: this->particles) {
        const std::array<double, 3>& pos = p.getPos();
        const unsigned int ID_ = p.getID();

        this->posX[ID_] = pos[0];
        this->posY[ID_] = pos[1];
        this->posZ[ID_] = pos[2];
        this->types[ID_] = p.getType();
    }
}

//DYNAMICS 

void Solver::setTemp() {
    if (this->temp <= 0.0) return;

    std::cout << "TEMP APPLIED!" << std::endl;

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
    Solver::computeForcesCUDA();
    Solver::cpDeviceToHost();
    // 3 - EULER METHOD FOR UPDATE
    for (Particle& p : particles) {
        p.updateVel(this->dt);
        p.updatePos(this->dt);
        p.resetForce();
    }
    Solver::cpHostToDevice();

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
    cudaDeviceSynchronize();
}

const std::array<double, 3>& Solver::getFinalPosK(unsigned int ID) {
    //const std::array<double, 3> pos = this->particles[ID].getPos();
    return this->particles[ID].getPos();
}

const double Solver::getKE(){
    return this->KE;
}