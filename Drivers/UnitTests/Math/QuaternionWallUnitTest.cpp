//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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

#include "Species/LinearViscoelasticFrictionBondedSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include "DPMBase.h"

using constants::pi;

/**
 * This driver tests the implementation of the rotating walls
 *
 * https://youtu.be/gObNthb3Szo
 */
class QuaternionWallUnitTest : public DPMBase
{
    //bonds all particle interactions, making them stick to the wall
    void actionsBeforeTimeStep() override
    {
        if (getTime()==0)
        {
            for (auto& i: interactionHandler)
            {
                dynamic_cast<BondedInteraction*>(i)->bond();
            }
        }
    }

//    void printTime() const override
//    {
//        std::cout
//         << "t " << std::setprecision(3) << std::left << std::setw(4) << getTime()
//         << ",\t q " << std::setprecision(3) << std::left << wallHandler.getObject(0)->getOrientation()
//         << ",\t n " << std::setprecision(3) << std::left << wallHandler.getObject(0)->getOrientation().getAxis()
//         << std::endl;
//    }

};

int main() {

    //setup problem
    QuaternionWallUnitTest q;
    q.setName("QuaternionWallUnitTest");
    q.setMin(-6,-6,-6);
    q.setMax(6,6,6);
    q.setDimension(3);
    q.setGravity(Vec3D(0,0,0));

    //set time parameters
    q.setTimeStep(0.01);
    q.setTimeMax(20*pi);
    q.setSaveCount(100);

    //set species
    auto s = q.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionBondedSpecies());
    s->setDensity(6/pi);
    s->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(1,0.5,0.5,1);
    s->setBondForceMax(.5);
    s->setSlidingFrictionCoefficient(0.5);
    s->setRollingStiffness(.2*s->getStiffness());
    s->setRollingDissipation(.2*s->getDissipation());
    s->setRollingFrictionCoefficient(0.2);

    Mdouble angularVelocity = 0.025;

    logger(INFO, "Create cylindrical wall rotating around its axis");
    {
        Vec3D position = Vec3D(1.5, 0, 1.5); //position in top right corner

        //Create wall
        AxisymmetricIntersectionOfWalls w;
        w.setSpecies(s);
        w.setPosition(position);
        w.setAxis(Vec3D(0, 1, 0));
        w.addObject(Vec3D(-1, 0, 0), Vec3D(.5, 0, 0));
        w.addObject(Vec3D(0, 0, 1), Vec3D(0, 0, -1));
        w.addObject(Vec3D(0, 0, -1), Vec3D(0, 0, 1));
        w.setAngularVelocity(angularVelocity*Vec3D(0, 1, 0));
        q.wallHandler.copyAndAddObject(w);

        //Create particle
        SphericalParticle p(s);
        p.setRadius(0.50001);
        p.setVelocity(Vec3D(0, 0, 0));
        p.setPosition(position + Vec3D(0, 0, 1));
        q.particleHandler.copyAndAddObject(p);
        p.setPosition(position + Vec3D(0, 0, -1));
        q.particleHandler.copyAndAddObject(p);
        p.setPosition(position + Vec3D(1, 0, 0));
        q.particleHandler.copyAndAddObject(p);
        p.setPosition(position + Vec3D(-1, 0, 0));
        q.particleHandler.copyAndAddObject(p);
    }

    logger(INFO, "Create cubic wall rotating around its axis");
    {
        Vec3D position = Vec3D(-1.5, 0, -1.5); //bottom left corner

        //Create wall
        IntersectionOfWalls w;
        w.setSpecies(s);
        w.setPosition(position);
        w.addObject(Vec3D( 1, 0, 0), Vec3D(-.5, 0, 0));
        w.addObject(Vec3D(-1, 0, 0), Vec3D( .5, 0, 0));
        w.addObject(Vec3D( 0, 1, 0), Vec3D( 0,-.5, 0));
        w.addObject(Vec3D( 0,-1, 0), Vec3D( 0, .5, 0));
        w.addObject(Vec3D( 0, 0, 1), Vec3D( 0, 0,-.5));
        w.addObject(Vec3D( 0, 0,-1), Vec3D( 0, 0, .5));
        w.setAngularVelocity(angularVelocity*Vec3D(0, 1, 0));
        q.wallHandler.copyAndAddObject(w);

        //Create particle
        SphericalParticle p(s);
        p.setRadius(0.50001);
        p.setVelocity(Vec3D(0, 0, 0));
        p.setPosition(position + Vec3D(0, 0, 1));
        q.particleHandler.copyAndAddObject(p);
        p.setPosition(position + Vec3D(0, 0, -1));
        q.particleHandler.copyAndAddObject(p);
        p.setPosition(position + Vec3D(1, 0, 0));
        q.particleHandler.copyAndAddObject(p);
        p.setPosition(position + Vec3D(-1, 0, 0));
        q.particleHandler.copyAndAddObject(p);
    }

    logger(INFO, "Create square outer wall rotating around the other walls");
    {
        //Create Wall
        InfiniteWall w;
        w.setSpecies(s);
        w.setPosition(Vec3D(0,0,5));
        w.setAngularVelocity(angularVelocity*Vec3D(0, 1, 0));
        w.setNormal(Vec3D(0,0,1));
        w.setPrescribedPosition(
         [angularVelocity] (Mdouble time)
         { return 5.0*Vec3D(sin(angularVelocity*time),0.0,cos(angularVelocity*time)); });
        q.wallHandler.copyAndAddObject(w);

        ///\todo why do we need to set the position?
        w.setPosition(Vec3D(5,0,0));
        w.setNormal(Vec3D(1,0,0));
        w.setPrescribedPosition(
         [angularVelocity] (Mdouble time)
         { return 5.0*Vec3D(sin(angularVelocity*time+0.5*pi),0.0,cos(angularVelocity*time+0.5*pi)); });
        q.wallHandler.copyAndAddObject(w);

        w.setPosition(Vec3D(0,0,-5));
        w.setNormal(Vec3D(0,0,-1));
        w.setPrescribedPosition(
         [angularVelocity] (Mdouble time)
         { return 5.0*Vec3D(sin(angularVelocity*time+pi),0.0,cos(angularVelocity*time+pi)); });
        q.wallHandler.copyAndAddObject(w);

        w.setPosition(Vec3D(-5,0,0));
        w.setNormal(Vec3D(-1,0,0));
        w.setPrescribedPosition(
         [angularVelocity] (Mdouble time)
         { return 5.0*Vec3D(sin(angularVelocity*time+1.5*pi),0.0,cos(angularVelocity*time+1.5*pi)); });
        q.wallHandler.copyAndAddObject(w);

        //Create particle
        SphericalParticle p(s);
        p.setRadius(0.50001);
        p.setVelocity(Vec3D(0, 0, 0));
        p.setPosition(Vec3D(   0, 0,-4.5));
        q.particleHandler.copyAndAddObject(p);
        p.setPosition(Vec3D(-4.5, 0,   0));
        q.particleHandler.copyAndAddObject(p);
        p.setPosition(Vec3D(   0, 0, 4.5));
        q.particleHandler.copyAndAddObject(p);
        p.setPosition(Vec3D( 4.5, 0,   0));
        q.particleHandler.copyAndAddObject(p);
    }


    //q.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    //q.setParticlesWriteVTK(true);

    //solve
    q.solve();

    return 0;
}
