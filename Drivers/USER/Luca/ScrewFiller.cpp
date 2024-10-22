//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://www.MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//! [CST:headers]
#include "Mercury3D.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Walls/InfiniteWallWithHole.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/Screw.h"
#include "Helicoid.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
//! [CST:headers]

//! [CST:class]
class ScrewFiller : public Mercury3D {
    
private:
    const double rMax = .35/2.;
    const double rMin = .19/2.;
    const double length = 1.5;
    const double thickness = 0.06;
    
    const int nP = 3000;

    void setupInitialConditions() override {
        // gravity, particle radius
        setGravity(Vec3D(0.0, 0.0, 0.0));
        particleRadius = 0.015;
        
        // set problem geometry
        setXMax(rMax);
        setYMax(rMax);
        setZMax(length);
        setXMin(-rMax);
        setYMin(-rMax);

        // set problem species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        species->setDensity(2000);
        
        double stiffness = 1.0e4;
        species->setStiffnessAndRestitutionCoefficient(stiffness, 0.8, pow(particleRadius, 3) * constants::pi * 4.0 / 3.0 *species->getDensity());
        species->setSlidingStiffness(2.0 / 7.0 * stiffness);
        species->setSlidingFrictionCoefficient(0.1);
        species->setRollingStiffness(2.0 / 5.0 * stiffness);
        species->setRollingFrictionCoefficient(0.1);
        species->setTorsionStiffness(2.0 / 5.0 * stiffness);
        species->setTorsionFrictionCoefficient(0.1);

        std::cout << "\nCollision time: " << species->getCollisionTime(pow(particleRadius, 3) * constants::pi * 4.0 / 3.0 *species->getDensity()) << std::endl;
        
        // set problem walls
        // walls with hole along the z direction
//        leftWall = wallHandler.copyAndAddObject(InfiniteWallWithHole());
//        leftWall->set(Vec3D(0, 0, 1), getZMin(), rMax);
//        rightWall = wallHandler.copyAndAddObject(InfiniteWallWithHole());
//        rightWall->set(Vec3D(0, 0, 1), getZMax(), rMax);
        
        // infinite walls at the z edges of the screw
//        zMinWall = wallHandler.copyAndAddObject(InfiniteWall());
//        zMinWall->set(Vec3D(0, 0, -1), getZMin());
//        zMaxWall = wallHandler.copyAndAddObject(InfiniteWall());
//        zMaxWall->set(Vec3D(0, 0, 1), getZMax());
        
        // periodic boundary in the screw's axis direction
        b0 = boundaryHandler.copyAndAddObject(PeriodicBoundary());
        b0->set(Vec3D(0,0,1), getZMin(), getZMax());
        
        // outer cylinder (the screw case) [I guess I have to use the axisymmetric wall]
        screwCase = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        screwCase->setPosition(Vec3D(0,0,0));
        screwCase->setOrientation(Vec3D(0, 0, 1));
        screwCase->addObject(Vec3D(1,0,0),Vec3D(rMax,0,0));
        
        // inner cylinder (the screw shaft) [I guess I have to use the axisymmetric wall]
        screwShaft = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        screwShaft->setPosition(Vec3D(0,0,0));
        screwShaft->setOrientation(Vec3D(0, 0, 1));
        screwShaft->addObject(Vec3D(-1,0,0),Vec3D(rMin,0,0));

        // creation of the screw and setting of its properties
        screw = wallHandler.copyAndAddObject(Helicoid());
        // set(Start position, Length, Radius, Number of turns, Rotation speed, Thickness)
        screw->set(Vec3D(0, 0, 0), 0.5*length, rMax, 2.0, constants::pi, 0.);
        
        // particle creation
        // particles are created randomly in a cylindrical shell of radius r in [rMin; rMax]
        // particle radius is uniformly randomly picked depending on the dispersity
        // dispersity of 1 means radius +- 1*radius
        Mdouble distance;
        Mdouble r, phi;
        Mdouble dispersity = 0.1;
        
        particleHandler.clear();
        SphericalParticle p0;
        p0.setSpecies(species);
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        
        for (int i = 0; i < nP; i++)
        {
            p0.setRadius(particleRadius*random.getRandomNumber(1.0-dispersity,1.0+dispersity));
            r = random.getRandomNumber(rMin + p0.getRadius(), rMax - p0.getRadius());
            phi = random.getRandomNumber(.0, 2.0*constants::pi);
            p0.setPosition(Vec3D(r*cos(phi),r*sin(phi),random.getRandomNumber(.0, length)));
            particleHandler.copyAndAddObject(p0);
        }
        
        /*
         ToDo:
         - should check for ovelap during the insertion (not really useful if blade AND particles are inflated)
         - inflate the particle radius
         - inflate the thickness also
         - correct the magic number in the radius change
         
         int particlesPlaced = 0;
         int maxAttempts = 10;
         Vec3D normal;
        */

    }

    void actionsOnRestart() override
    {
        // set problem species
        species = dynamic_cast<LinearViscoelasticFrictionSpecies*>(speciesHandler.getObject(0));
        b0 = dynamic_cast<PeriodicBoundary*>(boundaryHandler.getObject(0));
        screwCase = dynamic_cast<AxisymmetricIntersectionOfWalls*>(wallHandler.getObject(0));
        screwShaft = dynamic_cast<AxisymmetricIntersectionOfWalls*>(wallHandler.getObject(1));
        screw = dynamic_cast<Helicoid*>(wallHandler.getObject(2));
    }

    
    //! [CST:beforetime]
    void actionsBeforeTimeStep() override
    {
        // continuously changes the radius of the blade from rMin to rMax in a time timeMax - 1.5
        // according to: (t - t0)/(tMax - t0) * (rMax - rMin) + rMin
        if (getTime() > 0.0 && getTime() < 2.0)
        {
            screw->setThickness((getTime() - 0.0)/(2.0 - 0.0)*thickness);
        }
        
        // activates gravity after a given time
        
        if (getTime() > 2.0 && getTime() < 4.0)
        {
            setGravity(Vec3D(0.0, -9.81*(getTime() - 2.0)/(4.0 - 2.0), 0.0));
        }
        
        // rotates the helicoid
        if (getTime() > 4.0) screw->move_time(getTimeStep());
    }
    //! [CST:beforetime]

    BaseWall* readUserDefinedWall(const std::string& type) const override
    {
        return new Helicoid();
    }



    //! [CST:datamembers]
public:
    double particleRadius;
    LinearViscoelasticFrictionSpecies* species;
    Helicoid* screw;
//    InfiniteWallWithHole* leftWall;
//    InfiniteWallWithHole* rightWall;
//    InfiniteWall *zMinWall, *zMaxWall;
    AxisymmetricIntersectionOfWalls* screwCase;
    AxisymmetricIntersectionOfWalls* screwShaft;
    PeriodicBoundary* b0;
    //! [CST:datamembers]
};
//! [CST:class]

//! [CST:main]
int main(int argc, char *argv[]) {
    
    // create CoilFiller object
    ScrewFiller problem;
    
    // set some basic problem properties
    problem.setName("ScrewFiller");
    problem.setSystemDimensions(3);
    problem.setTimeStep(0.004/50.0/10.0);
    problem.setTimeMax(10.0);
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(1000, problem.getTimeMax(), problem.getTimeStep()));
    
    // actually solving the problem
    problem.solve(argc,argv);
    
}
//! [CST:main]
// the end
