#include "GCG.h"
#include "Species/LinearViscoelasticSlidingFrictionReversibleAdhesiveSpecies.h"

/**
 * Simulates the mixing
 */
int main()
{
    //set all necessary parameters
    Mdouble radius = 5e-3; // m
    Mdouble poly = 0.25;
	Mdouble maxFillRate = 100e-6; // m^3/s
    Mdouble rpm = 120;

    Mdouble density = 1500; //kg/m^3
    Mdouble collisionTime = 0.003; // s
    Mdouble restitution = 0.2;
    Mdouble friction = 0.5;

    //define the particle properties
    const Mdouble mass = 4./3.*density*pi*cubic(radius);
    const Mdouble minMass = 4./3.*density*pi*cubic((1.-1.5*poly)*radius);

    LinearViscoelasticSlidingFrictionReversibleAdhesiveSpecies particleSpecies;
    particleSpecies.setDensity(density);
    particleSpecies.setCollisionTimeAndRestitutionCoefficient(collisionTime, restitution, minMass);
    particleSpecies.setSlidingFrictionCoefficient(friction);
    particleSpecies.setSlidingStiffness(2.0/7.0*particleSpecies.getStiffness());
    particleSpecies.setSlidingDissipation(2.0/7.0*particleSpecies.getDissipation());
    particleSpecies.setAdhesionStiffness(particleSpecies.getStiffness());
    particleSpecies.setAdhesionForceMax(1*mass);
    auto wallSpecies = particleSpecies;

    //setup simulation
    GCG gcg(maxFillRate, rpm, radius, poly, &particleSpecies, &wallSpecies);
    gcg.removeOldFiles();
    gcg.setTimeStep(0.02 * collisionTime);
    gcg.setSaveCount(0.01/gcg.getTimeStep()); //save every 0.02s
    gcg.setTimeMax(120);

    gcg.setParticlesWriteVTK(true);
    gcg.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    gcg.solve();
}
