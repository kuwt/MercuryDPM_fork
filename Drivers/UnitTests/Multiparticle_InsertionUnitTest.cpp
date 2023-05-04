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

#include "Mercury3D.h"
#include "Boundaries/FixedClusterInsertionBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
//#include "Species/LinearPlasticViscoelasticFrictionSpecies.h"
#include "Species/LinearViscoelasticSlidingFrictionBondedSpecies.h"
#include "Walls/InfiniteWall.h"


/*!
This script contains the example to create multi-particles, which collide each other.
 */
class MultiParticlesInsertion : public Mercury3D
{
public:

    MultiParticlesInsertion(){

        setHGridMaxLevels(2);
        setMin(Vec3D(0, 0, 0));
        setMax(Vec3D(0.4,0.4,0.4));

        setName("MultiParticlesInsertion");
        setGravity(Vec3D(0, 0, -9.81));
        setXBallsAdditionalArguments("-solidf -v0");
        dataFile.setSaveCount(10);

        //Particle species
        particleSpecies = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionBondedSpecies());
        particleSpecies->setDensity(2000); //set the species density
        particleSpecies->setDissipation(0.01); //set the dissipation.
//        particleSpecies->setPlasticParameters(1e3,1e3*5,1e3,0.1);//set the spring stiffness.
        particleSpecies->setStiffness(1e3);
        particleSpecies->setBondForceMax(1.0e2);
    }

    LinearViscoelasticSlidingFrictionBondedSpecies* particleSpecies;

    void setupInitialConditions() override {

        Mdouble radiusParticle = 0.03;

        SphericalParticle p0,p1,p2,p3;
        p0.setRadius(radiusParticle);
        p1.setRadius(radiusParticle);
        p2.setRadius(radiusParticle);
        p3.setRadius(radiusParticle);

        p0.setSpecies(particleSpecies);
        p1.setSpecies(particleSpecies);
        p2.setSpecies(particleSpecies);
        p3.setSpecies(particleSpecies);

        p0.setPosition(Vec3D(0.5 * getXMax()*0.8, 0.5 * getYMax(), 0.5*getZMax()));
        p1.setPosition(Vec3D(0.5 * getXMax(), 0.5 * getYMax(), 0.5*getZMax()));

        p2.setPosition(Vec3D(0.5 * getXMax()*0.8, 0.5 * getYMax(), 0.1*getZMax()));
        p3.setPosition(Vec3D(0.5 * getXMax(), 0.5 * getYMax(), 0.1*getZMax()));

        p0.setVelocity(Vec3D(1.0, 0.0, 0.0));
//        p1.setVelocity(Vec3D(-1.0, 0.0, 0.0));

        particleHandler.copyAndAddObject(p0);
        particleHandler.copyAndAddObject(p1);
        particleHandler.copyAndAddObject(p2);
        particleHandler.copyAndAddObject(p3);

        //bond certain particles
        auto i = dynamic_cast<BondedInteraction*>(interactionHandler.getInteraction(particleHandler.getObject(0),particleHandler.getObject(1),0));
        i->bond();
        i = dynamic_cast<BondedInteraction*>(interactionHandler.getInteraction(particleHandler.getObject(2),particleHandler.getObject(3),0));
        i->bond();

        //set wall properties
        wallHandler.clear();
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0, 0, -1), Vec3D(0, 0, getZMin()));
        wallHandler.copyAndAddObject(w0);

    }

    void printTime() const override
    {
        logger(INFO,"t=%, tMax=%, N=%", getTime(),getTimeMax(), particleHandler.getSize());
    }

};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    logger(INFO,"Inserting multi-particles.");

    MultiParticlesInsertion insertion;

    insertion.setTimeStep(1e-4);
    insertion.setTimeMax(0.5);
    insertion.solve();

    helpers::check(4, insertion.particleHandler.getSize(), 0.1, "Number of particles check");
//    helpers::check(0, insertion.interactionHandler.getSize(), 0.1, "Number of interactions check");

}
