#include <Particles/LiquidFilmParticle.h>
#include "Mercury3D.h"
#include "Species/LinearViscoelasticFrictionLiquidMigrationWilletSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Walls/BasicIntersectionOfWalls.h"
#include "Boundaries/DropletBoundary.h"
using constants::pi;
using mathsFunc::cubic;

/**
 * A new test for the insertiion of liquid droplets.
 * Note, you need to create a species for the droplets that has liquidVolumeMax-0
 */
int main()
{
    //general properties
    Mercury3D dpm;
    dpm.setName("DropletBoundarySelfTest");
    double collisionTime = 0.005;
    dpm.setTimeStep(collisionTime/50.);
    dpm.setTimeMax(1.0);
    dpm.setSaveCount(200);
    dpm.setGravity({0,0,-1});
    dpm.setMax({1,1,0.3});
    dpm.setMin({0,0,0});

    //contact law and density
    LinearViscoelasticFrictionLiquidMigrationWilletSpecies species;
    species.setDensity(6./pi);
    species.setCollisionTimeAndRestitutionCoefficient(0.5e-2, 0.5, 1);
    species.setLiquidBridgeVolumeMax(1e-3);
    auto s = dpm.speciesHandler.copyAndAddObject(species);
    species.setLiquidBridgeVolumeMax(0);
    dpm.speciesHandler.copyAndAddObject(species);

    //walls
    dpm.wallHandler.copyAndAddObject(InfiniteWall({0,0,-1},dpm.getMin(),s));

    //two rows of particles
    double n = 3;
    double dx = (dpm.getXMax()-dpm.getXMin())/n;
    LiquidFilmParticle particle;
    particle.setSpecies(s);
    particle.setRadius(0.3*dx);
    particle.fixParticle();
    for (double x = dpm.getMin().X; x<=dpm.getMax().X+dx/2; x+=dx) {
        particle.setPosition({x,dpm.getYMax(),dpm.getZMin()});
        dpm.particleHandler.copyAndAddObject(particle);
    }
    for (double x = dpm.getMin().X+dx/2; x<=dpm.getMax().X; x+=dx) {
        particle.setPosition({x,dpm.getYMax(),dpm.getZMin()+0.5*dx});
        dpm.particleHandler.copyAndAddObject(particle);
    }

    //droplets
    double flowRate = 1e-3/60.; //m^3/s
    double dropletRadius = 4e-3;
    double dropletVolume = std::pow(dropletRadius,3)*constants::pi/6.0;
    DropletBoundary d;
    d.setGenerateDroplets([&dpm,flowRate,dropletRadius,dropletVolume] (DropletBoundary& d) {
        static double accumulatedDropletVolume = 0;
        while (accumulatedDropletVolume < flowRate*dpm.getTime()) {
            Vec3D position {dpm.random.getRandomNumber(dpm.getXMin(), dpm.getXMax()), dpm.getYMax(), dpm.getZMax()};
            Vec3D velocity {0,0,0};
            d.droplets_.emplace_back(DropletBoundary::Droplet{position, velocity, dropletRadius});
            accumulatedDropletVolume +=dropletVolume;
        }
    });
    dpm.boundaryHandler.copyAndAddObject(d);

    dpm.setParticlesWriteVTK(true);
//    dpm.setWallsWriteVTK(FileType::ONE_FILE);
    dpm.setInteractionsWriteVTK(true);
    dpm.boundaryHandler.setWriteVTK(true);
    dpm.solve();
    return 0;
}
