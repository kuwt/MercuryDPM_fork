#include "Mercury3D.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Walls/BasicIntersectionOfWalls.h"
using constants::pi;
using mathsFunc::cubic;

/**
 *
 */
int main()
{
    //general properties
    Mercury3D test;
    test.setName("BasicIntersectionOfWallsSelfTest");
    test.setTimeStep(1e-4);
    test.setTimeMax(1.0);
    test.setSaveCount(100);
    test.setGravity({0,0,-1});
    test.setMax({.5,.5,.5});
    test.setMin({-.5,-.5,-.5});

    //contact law and density
    LinearViscoelasticSpecies species;
    species.setDensity(6./pi);
    species.setCollisionTimeAndRestitutionCoefficient(0.5e-2, 0.5, 1);
    auto ps = test.speciesHandler.copyAndAddObject(species);
    auto s = test.speciesHandler.copyAndAddObject(species);
    ps->setStiffness(0);

    //walls
    InfiniteWall wallDown;
    wallDown.setSpecies(s);
    wallDown.set({0,0,-1},{0,0,0});

    InfiniteWall wallBack;
    wallBack.setSpecies(s);
    wallBack.set({0,-1,0},{0,0,0});

    InfiniteWall wallLeft;
    wallLeft.setSpecies(s);
    wallLeft.set({-1,0,0},{0,0,0});

    //test.wallHandler.copyAndAddObject(wallLeft);
    //test.wallHandler.copyAndAddObject(wallDown);
    BasicIntersectionOfWalls wall;
    wall.setSpecies(s);
    wall.add(wallDown);
    wall.add(wallBack);
    wall.add(wallLeft);
    test.wallHandler.copyAndAddObject(wall);

    //particles
    BaseParticle particle;
    particle.setSpecies(ps);

    const bool simpletest = true;
    if (simpletest) {
        particle.setRadius(0.2);
        particle.setPosition({0.1,0.1,0.1});
        test.particleHandler.copyAndAddObject(particle);
        particle.setPosition({-0.4,0.1,0.1});
        test.particleHandler.copyAndAddObject(particle);
        particle.setPosition({0.1,0.1,-0.4});
        test.particleHandler.copyAndAddObject(particle);
        particle.setPosition({0.5,0.1,0.5});
        test.particleHandler.copyAndAddObject(particle);
        test.setTimeStep(1e-10);
        test.setTimeMax(test.getTimeStep());
    } else {
        const unsigned int n = 20;
        const Mdouble dx = (test.getXMax()-test.getXMin())/double(n);
        const Mdouble dz = (test.getZMax()-test.getZMin())/double(n);
        particle.setRadius(0.49999*std::min(dx,dz));
        for (Mdouble x=test.getXMin()+0.5*dx; x<test.getXMax(); x+=dx) {
            for (Mdouble z=test.getZMin()+0.5*dx; z<test.getZMax(); z+=dz) {
                particle.setPosition({x,0,z});
                test.particleHandler.copyAndAddObject(particle);
            }
        }
        test.setTimeStep(1e-10);
        test.setTimeMax(test.getTimeStep());
        test.setXBallsAdditionalArguments(" ");
    }

    //test.write(std::cout);
    test.solve();
    return 0;
}
