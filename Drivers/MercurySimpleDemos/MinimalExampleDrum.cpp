

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <CG/Functions/Lucy.h>
#include <CG/TimeAveragedCG.h>
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include <Walls/SimpleDrumSuperquadrics.h>
#include "Chute.h"
#include "Walls/InfiniteWall.h"
#include "Boundaries/PeriodicBoundary.h"
#include "MercuryTime.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"


class MinimalExampleDrum : public Mercury3D
{
public:
    MinimalExampleDrum()
    {
        LinearViscoelasticFrictionSpecies species;
        species.setCollisionTimeAndRestitutionCoefficient(0.005, 0.1, 1);
        species.setDensity(constants::pi / 6);
        
        species.setSlidingDissipation(species.getDissipation() * 2. / 7.);
        species.setSlidingStiffness(species.getStiffness() * 2. / 7.);
        species.setSlidingFrictionCoefficient(0.5);
        
        species.setRollingStiffness(species.getStiffness() * 2.0 / 7.0);
        species.setRollingFrictionCoefficient(0.5);
        species.setRollingDissipation(species.getDissipation() * 2. / 7.);
        speciesHandler.copyAndAddObject(species);
        
        setTimeStep(1e-4);
        setTimeMax(10);
        setMin(-10, -0, -10);
        setMax(10, 2, 10);
        setGravity({0, 0, -10});
        setSaveCount(500);
    }
    
    void actionsAfterTimeStep() override
    {
        static bool isActivated = false;
        if (!isActivated && getTime() > 5)
        {
            for (BaseWall* w : wallHandler)
            {
                w->setAngularVelocity({0, -0.1, 0});
            }
            isActivated = true;
        }
        
    }
    
    void setupInitialConditions() override
    {
        
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(1.0);
        
        /*p0.setPosition(Vec3D(0.0, 0.0, -8.9));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);*/
        for (unsigned int i = 0; i < 5 ; ++i)
        {
            unsigned failCounter = 0;
            do
            {
                Mdouble r = random.getRandomNumber(-0.9, 0.9) * 10;
                Mdouble theta = random.getRandomNumber(0, constants::pi * 2.0);
                Mdouble y = random.getRandomNumber(getYMin(), getYMax());
                Vec3D position;
                position.X = r * std::cos(theta);
                position.Z = y;
                position.Y = r * std::sin(theta);
        
                p0.setPosition(position);
                p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        
                failCounter++;
        
                if (failCounter == 10000)
                {
                    break;
                }
        
            } while (!checkParticleForInteraction(p0));
    
            particleHandler.copyAndAddObject(p0);
        }
        
        SimpleDrumSuperquadrics w0;
        w0.setRadius(10);
        w0.setOrientationViaNormal(Vec3D(0, 1, 0));
        w0.setPosition(Vec3D(0, 0, 0));
        w0.setAngularVelocity(Vec3D(0, 0, 0));
        w0.setSpecies(speciesHandler.getObject(0));
    
        AxisymmetricIntersectionOfWalls w;
        w.addObject(Vec3D(1,0,0), Vec3D(10,0,0));
        w.setOrientationViaNormal(Vec3D(0, -1, 0));
        w.setPosition(Vec3D(0, 0, 0));
        w.setAngularVelocity(Vec3D(0, 0, 0));
        w.setSpecies(speciesHandler.getObject(0));
        wallHandler.copyAndAddObject(w);
        
        PeriodicBoundary b;
        b.set({0,1,0}, getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b);
        
    }
};


int main()
{
    MinimalExampleDrum problem;
    problem.setName("MinimalExampleDrum");
    problem.setParticlesWriteVTK(true);
    //problem.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    problem.solve();
    return 0;
}