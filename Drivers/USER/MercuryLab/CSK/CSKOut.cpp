#include "CSK.h"

/**
 * Simulates the outflow stage
 */
int main()
{
    Mdouble radius = 30e-3;
    Mdouble fillVolume = 0.1;
    Mdouble outflowRate = 0.001;
    Mdouble rpm = 6;

    Mdouble density = 300;
    Mdouble collisionTime = 0.01;
    Mdouble restitution = 0.5;
    Mdouble friction = 0.5;

    //define the particle properties
    LinearViscoelasticSlidingFrictionSpecies particleSpecies;
    Mdouble mass = density*4./3.*pi*cubic(radius);
    particleSpecies.setDensity(density);
    particleSpecies.setCollisionTimeAndRestitutionCoefficient(collisionTime, restitution, mass);
    particleSpecies.setSlidingFrictionCoefficient(friction);
    particleSpecies.setSlidingStiffness(2.0/7.0*particleSpecies.getStiffness());
    particleSpecies.setSlidingDissipation(2.0/7.0*particleSpecies.getDissipation());
    //define the properties of the particle-wall contacts
    auto particleWallSpecies = particleSpecies; //set species equal to the particle species

    //setup simulation
    CSK csk(&particleSpecies, &particleWallSpecies, fillVolume, outflowRate, rpm, radius);
    csk.setName("CSKOut");
    csk.setTimeStep(0.1 * collisionTime); //10 times too high collision rate
    csk.setSaveCount(5.*collisionTime/csk.getTimeStep()); //save every 20 collisions
    csk.setTimeMax(12.0);

    csk.setParticlesWriteVTK(true);
    csk.setWallsWriteVTK(FileType::NO_FILE);
    //csk.wallTest(5e-3);
    csk.solve();
}
