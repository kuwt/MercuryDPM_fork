//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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

#include "Mercury3D.h"
#include "Particles/BaseParticle.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Walls/InfiniteWall.h"
#include "Walls/InfiniteWallWithHole.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "Walls/Screw.h"
#include "Shaft.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"

/*
 ToDo:
 - should check for ovelap during the insertion (not really useful if blade AND particles are inflated)
 - inflate the particle radius
 - inflate the thickness also
 - correct the magic number in the radius change
 
 - correct the particle-screw_side interaction
 - adjust material properties
 - clean everything and make it consistent
 - fix the collision rule
 
 */

class ShaftTester : public Mercury3D {
    
private:
    const double radius = 0.2;
    const double length = 1.0;
    
    const int nP = 500;
    
    void setupInitialConditions() {
        // gravity, particle radius
        setGravity(Vec3D(0.0,0.0,-9.81));
        particleRadius = 0.025;
        
        // set problem geometry
        setXMax(2.0*radius);
        setYMax(2.0*radius);
        setZMax(radius);
        setXMin(-2.0*radius);
        setYMin(-2.0*radius);
        setZMin(.0*radius);
        
        // set problem species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        species->setDensity(2000);
        
        double stiffness = 1.0e4;
        species->setStiffnessAndRestitutionCoefficient(stiffness, 0.8, pow(particleRadius, 3) * constants::pi * 4.0 / 3.0 *species->getDensity());
        species->setSlidingStiffness(2.0 / 7.0 * stiffness);
        species->setSlidingFrictionCoefficient(0.5);
        species->setRollingStiffness(2.0 / 5.0 * stiffness);
        species->setRollingFrictionCoefficient(0.5);
        species->setTorsionStiffness(2.0 / 5.0 * stiffness);
        species->setTorsionFrictionCoefficient(0.5);
        
        std::cout << "\nCollision time: " << species->getCollisionTime(pow(particleRadius, 3) * constants::pi * 4.0 / 3.0 *species->getDensity()) << std::endl;
        
        // periodic boundary in the screw's axis direction
        
        xMinWall = wallHandler.copyAndAddObject(InfiniteWall());
        xMinWall->set(Vec3D(-1,0,0),Vec3D(getXMin(),0,0));
        
        xMaxWall = wallHandler.copyAndAddObject(InfiniteWall());
        xMaxWall->set(Vec3D(1,0,0),Vec3D(getXMax(),0,0));
        
        yMinWall = wallHandler.copyAndAddObject(InfiniteWall());
        yMinWall->set(Vec3D(0,-1,0),Vec3D(0,getYMin(),0));
        
        yMaxWall = wallHandler.copyAndAddObject(InfiniteWall());
        yMaxWall->set(Vec3D(0,1,0),Vec3D(0,getYMax(),0));
        
        zMinWall = wallHandler.copyAndAddObject(InfiniteWall());
        zMinWall->set(Vec3D(0,0,-1),Vec3D(0,0,getZMin()));
        
        zMaxWall = wallHandler.copyAndAddObject(InfiniteWall());
        zMaxWall->set(Vec3D(0,0,1),Vec3D(0,0,getZMax()));
        
//        // outer cylinder (the screw case) [I guess I have to use the axisymmetric wall]
//        screwCase = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
//        screwCase->setSpecies(species);
//        screwCase->setPosition(Vec3D(0,0,0));
//        screwCase->setOrientation(Vec3D(0, 0, 1));
//        screwCase->addObject(Vec3D(1,0,0),Vec3D(2.0*radius,0,0));
        
        wallHandler.clear();
        
        screwCase = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        screwCase -> setSpecies(species);
        screwCase -> setPosition(Vec3D(0,0,0));
        screwCase -> setOrientation(Vec3D(0.0,1.0,0.0));
        screwCase -> addObject(Vec3D(1,0,0), Vec3D(2.0*radius,0.0,0.0));
//        screwCase -> setAngularVelocity(Vec3D(0.0,constants::pi,0.0));
        
        // inner cylinder
        screwShaft = wallHandler.copyAndAddObject(AxisymmetricIntersectionOfWalls());
        screwShaft->setSpecies(species);
        screwShaft->setPosition(Vec3D(0,0,0));
        screwShaft->setOrientation(Vec3D(0,1,0));
        screwShaft->addObject(Vec3D(-1,0,0),Vec3D(-radius,0,0));
        
        // particle creation
        // particles are created randomly in a cylindrical shell of radius r in [rMin; rMax]
        // particle radius is uniformly randomly picked depending on the dispersity
        // dispersity of 1 means radius +- 1*radius
        Mdouble distance;
        Mdouble r, phi;
        Mdouble dispersity = 0.1;
        
        particleHandler.clear();
        BaseParticle p0;
        p0.setSpecies(species);
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        
        for (int i = 0; i < nP; i++)
        {
            p0.setRadius(particleRadius*random.getRandomNumber(1.0-dispersity,1.0+dispersity));
            r = random.getRandomNumber(radius + p0.getRadius(), 2.0*radius - p0.getRadius());
            phi = random.getRandomNumber(.0, 2.0*constants::pi);
//            p0.setPosition(Vec3D(0.0,2.0*rMax-p0.getRadius(),length/4.0+0.5*particleRadius));
            p0.setPosition(Vec3D(r*cos(phi),random.getRandomNumber(.0, length),r*sin(phi)));
//            p0.setPosition(Vec3D(random.getRandomNumber(getXMin(),getXMax()),random.getRandomNumber(rMax,getYMax()),random.getRandomNumber(.0, length)));
            particleHandler.copyAndAddObject(p0);
        }
        
        /*
         ToDo:
         - should check for ovelap during the insertion (not really useful if blade AND particles are inflated)
         - inflate the particle radius
         - inflate the thickness also
         - correct the magic number in the radius change
         */
    }
    
    //! [CST:beforetime]
    void actionsBeforeTimeStep()
    {
//        wallHandler.getObject(0)->setOrientation(Vec3D(0.0,1.0,0.0));
//        wallHandler.getObject(1)->setOrientation(Vec3D(0.0,1.0,0.0));
//        wallHandler.getObject(2)->setOrientation(Vec3D(0.0,1.0,0.0));
//        
//        // rotates the helicoid
//        if (getTime() > 2.0)
//        {
//            wallHandler.getObject(0)->setAngularVelocity(Vec3D(0.0,0.5*constants::pi,0.0));
//            wallHandler.getObject(1)->setAngularVelocity(Vec3D(0.0,0.5*constants::pi,0.0));
//            wallHandler.getObject(2)->setAngularVelocity(Vec3D(0.0,0.5*constants::pi,0.0));
//            
////            screwCase->setAngularVelocity(Vec3D(.0,0.5*constants::pi,.0));
//            
//            wallHandler.getObject(0)->setOrientation(Vec3D(0.0,1.0,0.0));
//            wallHandler.getObject(1)->setOrientation(Vec3D(0.0,1.0,0.0));
//            wallHandler.getObject(2)->setOrientation(Vec3D(0.0,1.0,0.0));
//        }
    }

public:
    double particleRadius;
    LinearViscoelasticFrictionSpecies* species;
    Shaft* shaft;
    InfiniteWall *xMinWall, *xMaxWall;
    InfiniteWall *yMinWall, *yMaxWall;
    InfiniteWall *zMinWall, *zMaxWall;
    AxisymmetricIntersectionOfWalls* screwShaft;
    AxisymmetricIntersectionOfWalls* screwCase;
};

int main(int argc UNUSED, char *argv[] UNUSED) {
    
    // create CoilFiller object
    ShaftTester problem;
    
    // set some basic problem properties
    problem.setName("ShaftTester");
    problem.setSystemDimensions(3);
    problem.setTimeStep(0.008/50.0);
    problem.setTimeMax(10.0);
    problem.setSaveCount(helpers::getSaveCountFromNumberOfSavesAndTimeMaxAndTimeStep(1000, problem.getTimeMax(), problem.getTimeStep()));
    
    // actually solving the problem
    problem.solve();
    
}

