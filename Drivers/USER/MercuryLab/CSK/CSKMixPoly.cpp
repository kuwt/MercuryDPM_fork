#include "CSK.h"
#include<sstream>
#include<chrono>
/**
 * Simulates the mixing stage
 */
int main()
{
    //start measuring time
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    startClock = std::chrono::system_clock::now();

    //set all necessary parameters
    Mdouble radius = 40e-3;   
	Mdouble fillVolume = 0.6;
    Mdouble outflowRate = 0.0; //a zero value also signifies that no outflow pipe is created (i.e. the outflow is closed)
    Mdouble rpm = 36;

    Mdouble density = 1010;
    Mdouble collisionTime = 0.005/1.5;
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
    CSK csk(&particleSpecies, &particleWallSpecies, fillVolume, outflowRate, rpm, radius,true,true);
    std::stringstream ss;
	std::string baseName = "phase3_mix";
	ss << baseName << "-r" << radius << "-rpm" << rpm;
	csk.setName(ss.str());
    csk.setTimeStep(0.02 * collisionTime); 
    csk.setSaveCount(200.*collisionTime/csk.getTimeStep()); //save every 500 collisions
    csk.setTimeMax(600.0);

    //define the properties of the second particle type
    auto particle2Species = particleSpecies; //set species equal to the first particle species
    auto particle2Particle1Species = particleSpecies; //set species equal to the first particle species
    auto particle2WallSpecies = particleSpecies; //set species equal to the first particle species
    //convert all particles with negative x-positions to the new species.
    auto conversionCondition = [] (const BaseParticle* const p) {return p->getPosition().X<0;};
    //add new particle species and convert particles
    csk.addSecondParticleSpecies(&particle2Species,&particle2Particle1Species,&particle2WallSpecies,conversionCondition);

    csk.setParticlesWriteVTK(false);
    //csk.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    //csk.wallTest(500e-3);
    //csk.guilloTest(6e-3);
    csk.solve();
    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << std::endl;
    logger(INFO, "Elapsed time for solving the PDE: % s", elapsed_seconds.count());
}
