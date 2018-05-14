//
// Created by irana on 5/24/17.
//

#include <Mercury3D.h>
#include <Walls/AxisymmetricIntersectionOfWalls.h>
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <Boundaries/PeriodicBoundary.h>
#include <Walls/CylindricalWall.h>
#include "ScrewRectangularSectionClean.h"
#include "../Helicoid05.h"

class ScrewRectangularTester : public Mercury3D
{
public:
    ScrewRectangularTester()
    {
        setMin({-10, -10, 0});
        setMax({10, 10, 10});
        setName("ScrewRectangularSelfTest");
    }
    
    void setupInitialConditions()
    {
        makeSpecies();
        makeScrew();
        makeBoundaries();
        insertParticles();
        setTimeMax(.01);
        setTimeStep(1e-4);
        setWallsWriteVTK(FileType::ONE_FILE);
        setParticlesWriteVTK(true);
        setSaveCount(1000);
        setGravity({0, 0, 0});
    }

private:
    // makes the screw
    void makeScrew()
    {
        // make the helicoid
        ScrewRectangularSectionClean helicoid;
        Vec3D startPosition(0, 0, 0);
        Mdouble screwLength = 5;
        Mdouble bladeRadius = 5;
        Mdouble numberOfTurns = 2;
        Mdouble screwSpeed = 10*constants::pi;
        Mdouble bladeThickness = 1;
        Mdouble shaftRadius = 1;
        helicoid.set(startPosition, screwLength, bladeRadius, shaftRadius, numberOfTurns, screwSpeed, bladeThickness, true);
        helicoid.setSpecies(speciesHandler.getObject(0));
        wallHandler.copyAndAddObject(helicoid);
        
        // make outer casing
        AxisymmetricIntersectionOfWalls casing;
        casing.setSpecies(speciesHandler.getObject(0));
        casing.setPosition({0, 0, 0});
        casing.setOrientation(Vec3D(0., 0., 1.));
        casing.addObject(Vec3D(1., 0., 0.), Vec3D(5, 0., 0.));
        casing.setAngularVelocity(Vec3D(0., 0., 0.));
        wallHandler.copyAndAddObject(casing);
    }
    
    void makeSpecies()
    {
        //dummy species
        LinearViscoelasticSlidingFrictionSpecies species;
        species.setDensity(6 / constants::pi);
        species.setCollisionTimeAndRestitutionCoefficient(0.005, 0.01, 1);
        
        species.setSlidingDissipation(species.getDissipation() * 2.0 / 7.0);
        species.setSlidingStiffness(species.getStiffness() * 2.0 / 7.0);
        species.setSlidingFrictionCoefficient(0.75);
        speciesHandler.copyAndAddObject(species);
    }
    
    void insertParticles()
    {
        Mdouble dispersity = 0.01;
        Mdouble particleRadius = 0.25;
        
        
        particleHandler.clear();
        BaseParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        
        for (unsigned int i = 0; i < 100000; i++)
        {
            p0.setRadius(particleRadius * random.getRandomNumber(1.0 - dispersity, 1.0 + dispersity));
            Mdouble r = random.getRandomNumber(1 + p0.getRadius(), 5 - p0.getRadius());
            Mdouble phi = random.getRandomNumber(.0, 2.0 * constants::pi);
            p0.setPosition(Vec3D(r * cos(phi), r * sin(phi), random.getRandomNumber(.0, 5)));
            if(checkParticleForInteraction(p0))
                particleHandler.copyAndAddObject(p0);
        }
    }
    
    void makeBoundaries()
    {
        PeriodicBoundary boundary;
        boundary.set({0, 0, 1}, 0, 5);
        boundaryHandler.copyAndAddObject(boundary);
    }
};

int main()
{
    ScrewRectangularTester problem;
    problem.solve();
    
}