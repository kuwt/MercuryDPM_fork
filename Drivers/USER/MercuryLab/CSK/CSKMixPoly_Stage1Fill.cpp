#include "CSK.h"
#include<sstream>
#include<chrono>
#include<iostream>
/**
 * Simulates the initial filling stage
 */

int main()
{
    
    //start measuring time
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    startClock = std::chrono::system_clock::now();
 
    Mdouble radius = 27.5e-3; //each particle has 5g
    Mdouble fillVolume = 0.6;
    Mdouble outflowRate = 0.0; //a zero value also signifies that no outflow pipe is created (i.e. the outflow is closed)
    Mdouble rpm = 2.0;

    unsigned nDoms = 6; //the number of domains used for parallelisation

    Mdouble density = 1010.0;
    Mdouble collisionTime = 0.005 / 2.5; //default 0.01
    Mdouble restitution = 0.6;
    Mdouble friction = 0.5;

    //define the particle properties
    LinearViscoelasticSlidingFrictionSpecies particleSpecies;
    Mdouble mass = 4./3.*density*pi*cubic(0.2*radius);
    
    particleSpecies.setDensity(density);
    particleSpecies.setCollisionTimeAndRestitutionCoefficient(collisionTime, restitution, mass);
    particleSpecies.setSlidingFrictionCoefficient(friction);
    particleSpecies.setSlidingStiffness(2.0/7.0*particleSpecies.getStiffness());
    particleSpecies.setSlidingDissipation(2.0/7.0*particleSpecies.getDissipation());
    //define the properties of the particle-wall contacts
    auto particleWallSpecies = particleSpecies; //set species equal to the particle species

    //setup simulation
    CSK csk(&particleSpecies, &particleWallSpecies, fillVolume, outflowRate, rpm, radius,false,true,true);
    //Set the number of domains for parallel decomposition
    csk.setNumberOfDomains({nDoms,1,1});
    std::stringstream ss;
	std::string baseName = "new_para_phase3_fill";
	ss << baseName << "-r" << radius << "-rpm" << rpm << "-nDomains" << nDoms;
	csk.setName(ss.str());
    csk.setTimeStep(0.02 * collisionTime); 
    csk.setSaveCount(10.0*collisionTime/csk.getTimeStep()); //normally 500* -- save every 500 collisiions
    csk.setTimeMax(6.0*600.0/rpm);
    //for quick testing
    //csk.setTimeMax(1.0);






    //csk.setParticlesWriteVTK(true);
    //csk.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    //csk.wallTest(500e-3);
    //csk.guilloTest(6e-3);
    csk.solve();
    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << std::endl;
    logger(INFO, "Elapsed time for solving the PDE: % s", elapsed_seconds.count());
    return 0; 
}
