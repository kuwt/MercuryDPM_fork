// Granular column collapse test (from 'MercurySimpleDemos/InsertionBoundarySelfTest.cpp')

#include <iostream>
#include "Mercury3D.h"
#include "Species/HertzianViscoelasticMindlinSpecies.h"
#include "Walls/InfiniteWall.h"
//#include "Particles/BaseParticle.h"
//#include "Species/LinearViscoelasticFrictionSpecies.h"

class GranularCollapse
        : public Mercury3D // [granularCollapse : class] from 'MercurySimpleDemos/InsertionBoundarySelfTest.cpp'
{
public:
    GranularCollapse() //from 'MercurySimpleDemos/HGridUnitTest.cpp'
    {
        species = new HertzianViscoelasticMindlinSpecies;
        speciesHandler.addObject(species);
    }
    
    void setupInitialConditions()
    {
        // [GranularCollapse : initial]
        setSystemDimensions(3);
        setGravity(Vec3D(0., 0., -9.8)); //body forces
        setXMin(-.025);
        setYMin(.0);
        setZMin(.0);
        setXMax(.025);
        setYMax(.05);
        setZMax(.35);
        // [GranularCollapse : initial]
        
        // [GranularCollapse : particles]
        SuperQuadric* ParticleToCopy;
        ParticleToCopy = new SuperQuadric;
        ParticleToCopy->setSpecies(species);
        // [GranularCollapse : particles]
        
        // [GranularCollapse : insertion]
        ParticleToCopy->setAxesAndExponents(radMax * radScale, radMax / std::sqrt(radScale),
                                            radMax / std::sqrt(radScale), 1.0, 1.0);
        ParticleToCopy->setInertia();
        ParticleToCopy->setVelocity(Vec3D(0., 0., 0.));
        unsigned nx = (getXMax() - getXMin()) / (2 * radMax * radScale);
        unsigned ny = (getYMax() - getYMin()) / (2 * radMax * radScale);
        unsigned nz = 1 + (nMax / nx / ny);
        Mdouble dx = (getXMax() - getXMin()) / (nx + 1);
        Mdouble dy = (getYMax() - getYMin()) / (ny + 1);
        Mdouble dz = 2 * radMax * radScale;
        for (unsigned k = 0; k < nz; k++)
            for (unsigned j = 0; j < ny; j++)
                for (unsigned i = 0; i < nx; i++)
                {
                    ParticleToCopy->setPosition(Vec3D(getXMin() + dx * i, //(i+random.getRandomNumber(-.5,.5)),
                                                      getYMin() + dy * j, //(j+random.getRandomNumber(-.5,.5)),
                                                      getZMin() + dz * k));
                    ParticleToCopy->setOrientationViaNormal(Vec3D(random.getRandomNumber(-1.0, 1.0),
                                                                  random.getRandomNumber(-1.0, 1.0),
                                                                  random.getRandomNumber(-1.0, 1.0)));
                    if (particleHandler.getNumberOfObjects() < nMax)
                        particleHandler.copyAndAddObject(ParticleToCopy);
                }
        // std::cout << "N = " << particleHandler.getNumberOfObjects() << std::endl;        
        // [GranularCollapse : insertion]
        
        // [GranularCollapse : walls]
        InfiniteWall wall;
        wall.setSpecies(species);
        wall.set(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(wall); //bottom (-1)
        wall.set(Vec3D(1, 0, 0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(wall); //front (-2)
        wall.set(Vec3D(-1, 0, 0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(wall); //back (-3)
        wall.set(Vec3D(0, -1, 0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(wall); //leftside (-4)
        wall.set(Vec3D(0, 1, 0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(wall); //rightside (-5)
        // [GranularCollapse : walls]
    }
    
    void actionsAfterTimeStep()
    {
        if (problemtype != 0 && problemtype != 3)
        {
            if (((getTime() < stopInsert && getTime() + getTimeStep() > stopInsert) ||
                 particleHandler.getNumberOfObjects() >= nMax) && !once)
            {
                std::cout << "Ending particle insertion" << std::endl;
                boundaryHandler.removeLastObject();
                once = true;
            }
        }
        if (getTime() < openGate && getTime() + getTimeStep() > openGate)
        {
            std::cout << "Shifting right wall rightward" << std::endl;
            dynamic_cast<InfiniteWall*>(wallHandler.getLastObject())->set(Vec3D(0, 1, 0), Vec3D(0, 2 + getYMax(), 0));
            //dynamic_cast<InfiniteWall*>(wallHandler.getLastObject())->setVelocity(Vec3D(0,0.1,0));
        }
    }
    
    // [GranularCollapse : datamembers]
    HertzianViscoelasticMindlinSpecies* species;
    int problemtype, nMax;
    double radMin, radMax, radScale, stopInsert, openGate;
    bool once = false;
    // [GranularCollapse : datamembers]
}; // [GranularCollapse : class]

int main(int argc, char* argv[]) // [GranularCollapse : main]
{
    GranularCollapse problem;
    problem.setName("granular_collapse_sqellipsoids");
    std::cout << "Modified from 'MercurySimpleDemos/InsertionBoundarySelfTest.cpp'" << std::endl;
    
    // [GranularCollapse : particle properties]
    problem.species->setDensity(950); // 2000;
    
    double input1 = 50; //1000; // 500; 1000; 15500;
    double input2 = 1.5; // magnifying scale (geq 1)
    if (argc > 1)
    {
        input1 = atof(argv[1]);
        input2 = atof(argv[2]);
    }
    problem.nMax = input1;
    problem.radScale = input2;
    std::cout << "Number of particles = " << problem.nMax << std::endl;
    std::cout << "Half-length scale = " << problem.radScale << std::endl;
    
    problem.radMax = 5.e-3; //1.e-3; 4.5e-3; 4.e-3
    problem.radMin = problem.radMax / problem.radScale / std::sqrt(problem.radScale);
    // [GranularCollapse : particle properties]
    
    // [GranularCollapse : contact properties]
    problem.species->setElasticModulusAndRestitutionCoefficient(1.e7, 0.7); //1.e5
//    //normal forces
//    problem.species->setStiffness(2e4); //1e5; 2e4; 8e3; 1e4
//    problem.species->setDissipation(0.5); //0.5; 0.25; 0.3; dissipation <= sqrt(2*stiffness*MinParticleMass) to avoid overdamping
//    //tangential forces : sliding
//    problem.species->setSlidingFrictionCoefficient(0.5);
//    problem.species->setSlidingStiffness(problem.species->getStiffness()*0.2);
//    problem.species->setSlidingDissipation(problem.species->getDissipation()*0.2);
//    //normal torques : spin
//    problem.species->setTorsionFrictionCoefficient(problem.species->getSlidingFrictionCoefficient()*0.05);
//    problem.species->setTorsionStiffness(problem.species->getStiffness()*0.1);
//    problem.species->setTorsionDissipation(problem.species->getDissipation()*0.05);
//    //tangential torques : rolling
//    problem.species->setRollingFrictionCoefficient(problem.species->getSlidingFrictionCoefficient()*0.05);
//    problem.species->setRollingStiffness(problem.species->getStiffness()*0.1);
//    problem.species->setRollingDissipation(problem.species->getDissipation()*0.05);
    // [GranularCollapse : contact properties]
    
    // [GranularCollapse : test normal forces]
    //set cout parameters
    std::cout.precision(3);
    std::cout << std::scientific;
    //(from 'MercurySimpleDemos/HourGlass3DDemo.cpp')
    //Calculates collision time for two copies of a particle of given dissipation_, k, effective mass
    double MinParticleMass = problem.species->getDensity() * 4. / 3. * constants::pi * mathsFunc::cubic(problem.radMin);
    std::cout << "MinParticleMass = " << MinParticleMass << std::endl;
    //Calculates collision time for two copies of a particle of given dissipation_, k, effective mass
    double tc = problem.species->getCollisionTime(2 * problem.radMin, problem.species->getDensity(), 1.0);
    std::cout << "tc = " << tc << std::endl;
//    //Calculates restitution coefficient for two copies of given dissipation_, k, effective mass
//    double r = problem.species->getRestitutionCoefficient(MinParticleMass);
//    std::cout << "r = " << r << std::endl;
    //Calculates the maximum relative velocity allowed for a normal collision of two particles of radius r and particle mass m (for higher velocities particles could pass through each other)
    //std::cout << "vmax=" << helpers::getMaximumVelocity(species->getStiffness(), HGgetSpecies(0)->getDissipation(), HG.MinParticleRadius, MinParticleMass) << std::endl;
    //unset cout parameters
    std::cout << std::defaultfloat;
    // [GranularCollapse : test normal forces]
    
    // [GranularCollapse : time parameters]
    problem.setTimeStep(tc / 50); // 1e-4; // (tc / 50);
    problem.setTimeMax(2.5); //5; //15;
    problem.openGate = .7 * problem.getTimeMax();
    problem.setSaveCount(500);
    // [GranularCollapse : time parameters]
    
    // [GranularCollapse : contact detection]
    problem.setHGridMaxLevels(2);
    // [GranularCollapse : contact detection]
    
    problem.setSuperquadricParticlesWriteVTK(true);
    problem.setWallsWriteVTK(FileType::ONE_FILE);
    problem.solve(); //argc,argv
    return 0;
} // [GranularCollapse : main]