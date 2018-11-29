#include "CSK.h"
#include<sstream>
#include<chrono>
#include<iostream>
/**
 * Simulates the mixing stage
 */

int main()
{
    //start measuring time
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    startClock = std::chrono::system_clock::now();

    Mdouble radius = 35.0e-3; //each particle has 5g
    Mdouble fillVolume = 0.6;
    Mdouble outflowRate = 0.0; //a zero value also signifies that no outflow pipe is created (i.e. the outflow is closed)
    Mdouble rpm = 14.4;

    Mdouble density = 1010.0;
    Mdouble collisionTime = 0.005; //default 0.01
    Mdouble restitution = 0.6;
    Mdouble friction = 0.5;

    //define the particle properties
    LinearViscoelasticSlidingFrictionSpecies particleSpecies;
    Mdouble mass = 4./3.*density*pi*cubic(0.4*radius);

    particleSpecies.setDensity(density);
    particleSpecies.setCollisionTimeAndRestitutionCoefficient(collisionTime, restitution, mass);
    particleSpecies.setSlidingFrictionCoefficient(friction);
    particleSpecies.setSlidingStiffness(2.0/7.0*particleSpecies.getStiffness());
    particleSpecies.setSlidingDissipation(2.0/7.0*particleSpecies.getDissipation());
    //define the properties of the particle-wall contacts
    auto particleWallSpecies = particleSpecies; //set species equal to the particle species

    //set the optional parameters
    bool delayedStartOfRotation=false;
    bool polydispersity = true;
    bool batchFill = true;
    Mdouble bagsPerRevolution = 200;
    Mdouble revolutions = 0.3;

    //setup simulation
    CSK csk(&particleSpecies, &particleWallSpecies, fillVolume, outflowRate, rpm, radius,
            delayedStartOfRotation, polydispersity, batchFill, bagsPerRevolution);
    std::stringstream ss;
    csk.setName("55initialFillSingleFillReducedRangeOnlyFilling_TestCase");
    logger(INFO,"Output file name: %",csk.getName());
    csk.setTimeStep(0.1 * collisionTime);
    csk.setSaveCount(20.0*collisionTime/csk.getTimeStep()); //save every 500 collisions
    csk.setTimeMax(revolutions*60.0/rpm);

    //uncomment to view in paraview
    //csk.setParticlesWriteVTK(true);
    //csk.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    csk.setNumberOfDomains({2,2,1});

    csk.solve();
    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << std::endl;
    logger(INFO, "Elapsed time for solving the PDE: % s", elapsed_seconds.count());
    return 0;
}
