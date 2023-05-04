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

#include <Mercury3D.h>
#include <Species/LinearViscoelasticSpecies.h>
#include <Walls/IntersectionOfWalls.h>
using constants::pi;

/*
* This is our problem description. Everything is set up here.
* We inherit from Mercury3D, since this gives full flexibility.
* For more predefined problems (for instance, chutes), please see the
* documentation.
*/
class Drum : public Mercury3D
{
public:
    /* We define our own 'setupInitialConditions' function here,
     * which defines all the specifics about our simulation here.
     */
    void setupInitialConditions() override
    {
        //particle radius
        Mdouble radius = 0.1;

        //The first step: set any properties which are always true for your system.
        // (for instance, if gravity remains constant, set it here)setName("Drum");
        setGravity(Vec3D(0,0,-9.8));
        setTimeMax(30);
        setTimeStep(0.002);
        //visualised length
        setXMax(100);
        //visualised height
        setZMax(50);
        //visualised width
        setYMax(5);
        setXMin(0);
        setZMin(0);
        setYMin(-5);

        //Now, decide what Species you need for your system.
        LinearViscoelasticSpecies species;
        species.setDensity(2000);
        Mdouble mass = species.getDensity() * 4.0 / 3.0 * pi * radius * radius * radius;
        species.setCollisionTimeAndRestitutionCoefficient(50*getTimeStep(),0.1,mass);
        LinearViscoelasticSpecies* s = speciesHandler.copyAndAddObject(species);

        //Add your walls below, and don't forget to set the species!

        //place a wall on the ground
        InfiniteWall ground;
        ground.setSpecies(speciesHandler.getObject(0));
        Vec3D position = Vec3D(0,0,0);
        Vec3D normal = Vec3D(0,0,-1);
        ground.set(normal,position);
        wallHandler.copyAndAddObject(ground);

        //place a wall at the front
        InfiniteWall front;
        front.setSpecies(speciesHandler.getObject(0));
        position = Vec3D(0,0,0);
        normal = Vec3D(-1,0,0);
        front.set(normal,position);
        wallHandler.copyAndAddObject(front);

        //place a wall at the back
        InfiniteWall back;
        back.setSpecies(speciesHandler.getObject(0));
        position = Vec3D(getXMax(),0,0);
        normal = Vec3D(1,0,0);
        back.set(normal,position);
        wallHandler.copyAndAddObject(back);

        //place an obstacle wall at 50% of distance, 50% height
        IntersectionOfWalls obstacle;
        obstacle.setSpecies(speciesHandler.getObject(0));
        //add the left side of the wall
        position = Vec3D(0.49*getXMax(),0,0.5*getZMax());
        normal = Vec3D(1,0,0);
        obstacle.addObject(normal,position);
        //add the upper side of the wall
        normal = Vec3D(0,0,-1);
        obstacle.addObject(normal,position);
        //add the right side of the wall
        position = Vec3D(0.51*getXMax(),0,0.5*getZMax());
        normal = Vec3D(-1,0,0);
        obstacle.addObject(normal,position);
        wallHandler.copyAndAddObject(obstacle);

        //place a target of 6 radii at 100% distance, 10% height
        IntersectionOfWalls target;
        target.setSpecies(speciesHandler.getObject(0));
        //add the left side of the wall
        position = Vec3D(0.9*getXMax(),0,0.1*getZMax());
        normal = Vec3D(1,0,0);
        target.addObject(normal,position);
        //add the upper side of the wall
        normal = Vec3D(0,0,-1);
        target.addObject(normal,position);
        //add the right side of the wall
        position = Vec3D(getXMax()-6*radius,0,0.1*getZMax());
        normal = Vec3D(-1,0,0);
        target.addObject(normal,position);
        wallHandler.copyAndAddObject(target);

        //Place the cannon ball at 0% distance, 0% height
        SphericalParticle cannonBall;
        cannonBall.setSpecies(speciesHandler.getObject(0));
        cannonBall.setRadius(radius);
        cannonBall.setPosition(Vec3D(0,0,radius));
        Mdouble inclination = 86.7 * pi / 180;
        Mdouble speed = 99.0;
        cannonBall.setVelocity(speed*Vec3D(cos(inclination),0,sin(inclination)));
        particleHandler.copyAndAddObject(cannonBall);

    }

    // uncomment to add drag force
//    void computeExternalForces(BaseParticle* CI) override
//    {
//        DPMBase::computeExternalForces(CI);
//        //Stokes drag: D = Cd * .5 * rho * V^2 * A
//        Mdouble dragCoefficient = 0.5*0.5*1*CI->getVelocity().getLength()*pi*square(CI->getRadius());
//        CI->addForce(- dragCoefficient * CI->getVelocity());
//    }
};

int main(int argc, char **argv)
{
    Drum problem;
    problem.setName("Cannon");
    problem.setSaveCount(10);
    problem.wallHandler.setWriteVTK(FileType::ONE_FILE);
    problem.setParticlesWriteVTK(true);
    problem.setXBallsAdditionalArguments(" -rmult 10");

    //solve the system, the single particle will now bounce on the plate
    problem.solve(argc, argv);
    return 0;
}
