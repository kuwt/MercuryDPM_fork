//made by irana, 2018

#include "Mercury3D.h"
#include <iostream>
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <cstdlib> //needed to be cygwin compatible (system not found)
#include "Boundaries/PeriodicBoundary.h"
#include "Boundaries/SubcriticalMaserBoundaryTEST.h"
#include "Walls/IntersectionOfWalls.h"
#include "Chute.h"


//Copy of MpiMaserChuteTest, with rough bottom.
class ParallelMaserWithRoughBottom : public Mercury3D
{
public:
    ///Standard maser without rough bottom: use setupInitialConditions
    void setupInitialConditions() override
    {
        //Add a particle
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        
        
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 8; j++)
            {
                for (int k = 0; k < 10; k++)
                {
                    p0.setPosition(Vec3D(1.0 + i * 2.0, 1.0 + j * 2.0, 1.0 + k * 2.0));
                    p0.setRadius(random.getRandomNumber(0.9, 1.0));
                    p0.setVelocity(Vec3D(3.0, random.getRandomNumber(-0.1, 0.1), 0.0));
                    particleHandler.copyAndAddObject(p0);
                }
            }
        }
        
        //Add maser boundary
        boundaryHandler.clear();
        SubcriticalMaserBoundaryTEST m0;
        m0.set(Vec3D(1.0, 0.0, 0.0), 0.0, 10.0);
        m0.setActivationTime(10);
        boundaryHandler.copyAndAddObject(m0);
        
        //Add periodic boundary at side
        PeriodicBoundary b0;
        b0.set(Vec3D(0.0, 1.0, 0.0), 0.0, 20.0);
        boundaryHandler.copyAndAddObject(b0);
        
        //Add bottom wall
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0, 0, -1), Vec3D(0.0, 0, 0));
        wallHandler.copyAndAddObject(w0);
    }
    
    void actionsOnRestart() override
    {
        setParticlesWriteVTK(true);
        dataFile.setFileType(FileType::ONE_FILE);
    }
    
    void printTime() const override
    {
        MPIContainer& communicator = MPIContainer::Instance();
        std::cout << "domain no. " << communicator.getProcessorID()
                  << ", t=" << std::setprecision(3) << std::left << std::setw(6) << getTime()
                  << ", tmax=" << std::setprecision(3) << std::left << std::setw(6) << getTimeMax()
                  << ", N = " << particleHandler.getNumberOfRealObjectsLocal()
                  << std::endl;
        std::cout.flush();
    }
};

int main(int argc, char* argv[])
{
    //Create the problem
    ParallelMaserWithRoughBottom problem;
    problem.setName("ParallelMaserWithRoughBottom");
    
    //Create a species
    LinearViscoelasticFrictionSpecies* species = problem.speciesHandler.copyAndAddObject(
            LinearViscoelasticFrictionSpecies());
    species->setDensity(1.0);
    
    //Set domain accordingly (domain boundaries are not walls!)
    problem.setXMin(-10.0);
    problem.setXMax(150.0);
    problem.setYMin(-10.0);
    problem.setYMax(100.0);
    problem.setZMin(-10.0);
    problem.setZMax(100.0);
    
    //specify particle properties
    species->setDensity(6.0 / constants::pi);
    
    //specify body forces
    problem.setGravity(Vec3D(1.0 * 0.37460659341, 0.0, -1.0 * 0.92718385456));
    
    //Set the number of domains for parallel decomposition
    //todo make this more than 1 for true parallel code
    problem.setNumberOfDomains({2, 1, 1});
    
    //specify contact properties
    //normal forces
    species->setStiffness(1e5);
    species->setDissipation(0.63);
    //tangential (sliding) forces
    species->setSlidingFrictionCoefficient(0.5);
    species->setSlidingStiffness(1.2e4);
    species->setSlidingDissipation(0.16);
    //tangential (rolling) torques
    species->setRollingFrictionCoefficient(0.2);
    species->setRollingStiffness(1.2e4);
    species->setRollingDissipation(6.3e-2);
    //normal (torsion/spin) torques
    species->setTorsionFrictionCoefficient(0.1);
    species->setTorsionStiffness(1.2e4);
    species->setSlidingDissipation(6.3e-2);
    
    //set other simulation parameters
    Mdouble MinParticleMass = species->getDensity() * 4.0 / 3.0 * constants::pi;
    //Mdouble tc = species->getCollisionTime(MinParticleMass);
    Mdouble tc = 0.02;
    logger(INFO, "tc: %", tc);
    problem.setTimeStep(tc / 50.0);
    problem.setTimeMax(50);
    problem.setSaveCount(250); //used to be 500
    
    //Set output to paraview
    problem.setParticlesWriteVTK(true);
    
    problem.solve(argc, argv);
    
    return 0;
}

