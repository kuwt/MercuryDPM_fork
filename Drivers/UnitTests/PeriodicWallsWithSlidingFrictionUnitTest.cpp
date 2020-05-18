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
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Mercury2D.h"
#include "Boundaries/PeriodicBoundary.h"
#include <Logger.h>

///This test is a UnitTest for:
///Periodic Particles in combination with HGrid
///Tangential Springs over periodic contacts and during transitions from periodic to normal and reverse
///Removal of Particles while tangential spring information has to be conserved

class PeriodicWallsWithSlidingFrictionUnitTest : public Mercury2D
{
    void setupInitialConditions() override {
        LinearViscoelasticFrictionSpecies species;
        const double stiffness = 1e5;
        species.setDensity(constants::pi / 6.0);
        species.setDensity(2500);
        species.setStiffness(stiffness);
        species.setSlidingStiffness(2.0 / 7.0 * stiffness);
        species.setSlidingFrictionCoefficient(0.5);
        species.setRollingStiffness(2.0 / 5.0 * stiffness);
        species.setRollingFrictionCoefficient(0.5);
        species.setTorsionStiffness(2.0 / 5.0 * stiffness);
        species.setTorsionFrictionCoefficient(0.5);
        speciesHandler.copyAndAddObject(species);

        //setParticleDimensions(2);

        setGravity(Vec3D(0.0, 0.0, 0.0));
        setTimeMax(8e-2);
        setTimeStep(2.0e-4);
        setFileType(FileType::NO_FILE);

        setXMin(0.0);
        setYMin(0.0);
        setXMax(1.0);
        setYMax(1.0);

        PeriodicBoundary b0;
        b0.set(Vec3D(1, 0, 0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(b0);
        b0.set(Vec3D(0, 1, 0), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(b0);

        SphericalParticle P0, P1, P2, P3, P4, P5, P6, P7;
        P0.setSpecies(speciesHandler.getObject(0));
        P1.setSpecies(speciesHandler.getObject(0));
        P2.setSpecies(speciesHandler.getObject(0));
        P3.setSpecies(speciesHandler.getObject(0));
        P4.setSpecies(speciesHandler.getObject(0));
        P5.setSpecies(speciesHandler.getObject(0));
        P6.setSpecies(speciesHandler.getObject(0));
        P7.setSpecies(speciesHandler.getObject(0));
        
        //Particle to be removed
        P0.setPosition(Vec3D(0.1, 0.1, 0.0));
        P1.setPosition(Vec3D(0.3, 0.3, 0.0));
        //Contacts starts normal becomes periodic
        P2.setPosition(Vec3D(0.6, 0.041, 0.0));
        P2.setAngularVelocity(Vec3D(0.3, 0.6, 0.9));
        P3.setPosition(Vec3D(0.4, 0.001, 0.0));
        //Normal case
        P4.setPosition(Vec3D(0.6, 0.82, 0.0));
        P4.setAngularVelocity(Vec3D(0.3, 0.6, 0.9));
        P5.setPosition(Vec3D(0.4, 0.78, 0.0));
        //Contact starts periodic becomes normal and periodic again
        P6.setPosition(Vec3D(0.02, 0.42, 0.0));
        P6.setAngularVelocity(Vec3D(0.3, 0.6, 0.9));
        P7.setPosition(Vec3D(0.82, 0.38, 0.0));

        P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        P1.setVelocity(Vec3D(0.0, 0.0, 0.0));
        P2.setVelocity(Vec3D(-1.0, 0.0, 0.0));
        P3.setVelocity(Vec3D(1.0, 0.0, 0.0));
        P4.setVelocity(Vec3D(-1.0, 0.0, 0.0));
        P5.setVelocity(Vec3D(1.0, 0.0, 0.0));
        P6.setVelocity(Vec3D(-1.0, 0.0, 0.0));
        P7.setVelocity(Vec3D(1.0, 0.0, 0.0));

        P0.setRadius(0.1);
        P1.setRadius(0.1);
        P2.setRadius(0.1);
        P3.setRadius(0.1);
        P4.setRadius(0.1);
        P5.setRadius(0.1);
        P6.setRadius(0.1);
        P7.setRadius(0.1);

        particleHandler.copyAndAddObject(P0);
        particleHandler.copyAndAddObject(P1);
        particleHandler.copyAndAddObject(P2);
        particleHandler.copyAndAddObject(P3);
        particleHandler.copyAndAddObject(P4);
        particleHandler.copyAndAddObject(P5);
        particleHandler.copyAndAddObject(P6);
        particleHandler.copyAndAddObject(P7);
    }

protected:

    void actionsBeforeTimeStep() override {
        static bool FirstParticleRemoved = false;
        static bool SecondParticlRemoved = false;
        if (getTime() > 3.5e-2 && !FirstParticleRemoved)
        {
            logger(VERBOSE, "Just before removing first particle");
            for (BaseInteraction* i : interactionHandler)
            {
                logger(VERBOSE, "Interaction between particle (id=% index=%) and particle (id=% index=%) delta=%)",
                       i->getP()->getId(), i->getP()->getIndex(), i->getI()->getId(),
                       i->getI()->getIndex(), i->getTimeStamp());
            }

            FirstParticleRemoved = true;
            particleHandler.removeObject(1);

            logger(VERBOSE, "Just after removing first particle");
            for (BaseInteraction* i : interactionHandler)
            {
                logger(VERBOSE, "Interaction between particle (id=% index=%) and particle (id=% index=%) delta=%)",
                       i->getP()->getId(), i->getP()->getIndex(), i->getI()->getId(),
                       i->getI()->getIndex(), i->getTimeStamp());
            }
            logger(VERBOSE, "What should have happened:");
            logger(VERBOSE, "First the particle (id=1,index=1) is removed");
            logger(VERBOSE, "This causes the particle (id=7,index=7) to move to (id=7,index=1)");
            logger(VERBOSE, "So the tangential spring witch was between particle (id=7,index=1) and (id=6,index=6) "
                    "should be move to the last particle and reversed");
        }

        if (getTime() > 4.5e-2 && !SecondParticlRemoved)
        {
            logger(VERBOSE, "Just before removing second particle");
            for (BaseInteraction* i : interactionHandler)
            {
                logger(VERBOSE, "Interaction between particle (id=% index=%) and particle (id=% index=%) delta=%)",
                       i->getP()->getId(), i->getP()->getIndex(), i->getI()->getId(),
                       i->getI()->getIndex(), i->getTimeStamp());
            }

            SecondParticlRemoved = true;
            particleHandler.removeObject(0);

            logger(VERBOSE, "Just after removing second particle");
            //TW: here this code seems to create a segmentation fault on my Mac
            for (BaseInteraction* i : interactionHandler)
            {
                logger(VERBOSE, "Interaction between particle (id=% index=%) and particle (id=% index=%) delta=%)",
                       i->getP()->getId(), i->getP()->getIndex(), i->getI()->getId(),
                       i->getI()->getIndex(), i->getTimeStamp());
            }
            logger(VERBOSE, "What should have happened:");
            logger(VERBOSE, "First the particle (id=0,index=0) is removed");
            logger(VERBOSE, "This causes the particle (id=6,index=6) to move to (id=0,index=0)");
            logger(VERBOSE, "So the tangential spring between particle (id=6,index=6) and (id=7,index=1) "
                    "should be move to the last particle and reversed");
        }
    }

    void printTime() const override {
    }
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    ///Start off by solving the default problem
    PeriodicWallsWithSlidingFrictionUnitTest periodicWallsWithSlidingFrictionUnitTest;
    periodicWallsWithSlidingFrictionUnitTest.setName("PeriodicWallsWithSlidingFrictionUnitTest");
    periodicWallsWithSlidingFrictionUnitTest.solve();

    BaseParticle* pNormal = periodicWallsWithSlidingFrictionUnitTest.particleHandler.getObject(4);
    BaseParticle* pCase1 = periodicWallsWithSlidingFrictionUnitTest.particleHandler.getObject(0);
    BaseParticle* pCase2 = periodicWallsWithSlidingFrictionUnitTest.particleHandler.getObject(2);
    Vec3D shift;

    shift = Vec3D(-0.58, -0.4, 0.0);
    if (!pCase1->getPosition().isEqualTo(pNormal->getPosition() + shift, 1e-10))
    {
        logger(FATAL, "E0 The particle is in the wrong position. It is %, however is should be %",
               pCase1->getPosition(), pNormal->getPosition() + shift);
    }
    if (!pCase1->getVelocity().isEqualTo(pNormal->getVelocity(), 1e-10))
    {
        logger(FATAL, "E1 The particle has the wrong velocity. It is %, however is should be %", pCase1->getVelocity(),
               pNormal->getVelocity());
    }

    shift = Vec3D(0.0, -0.779, 0.0);
    if (!pCase2->getPosition().isEqualTo(pNormal->getPosition() + shift, 1e-10))
    {
        logger(FATAL, "E2 The particle is in the wrong position. It is %, however is should be %",
               pCase2->getPosition(), pNormal->getPosition() + shift);
    }
    if (!pCase2->getVelocity().isEqualTo(pNormal->getVelocity(), 1e-10))
    {
        logger(FATAL, "E3 The particle has the wrong velocity. It is %, however is should be %", pCase2->getVelocity(),
               pNormal->getVelocity());
    }

    pNormal = periodicWallsWithSlidingFrictionUnitTest.particleHandler.getObject(5);
    pCase1 = periodicWallsWithSlidingFrictionUnitTest.particleHandler.getObject(1);
    pCase2 = periodicWallsWithSlidingFrictionUnitTest.particleHandler.getObject(3);

    shift = Vec3D(1.0 - 0.58, -0.4, 0.0);
    if (!pCase1->getPosition().isEqualTo(pNormal->getPosition() + shift, 1e-10))
    {
        logger(FATAL, "E4 The particle is in the wrong position. It is %, however is should be %",
               pCase1->getPosition(), pNormal->getPosition() + shift);
    }
    if (!pCase1->getVelocity().isEqualTo(pNormal->getVelocity(), 1e-10))
    {
        logger(FATAL, "E5 The particle has the wrong velocity. It is %, however is should be %", pCase1->getVelocity(),
               pNormal->getVelocity());
    }

    shift = Vec3D(0.0, 1.0 - 0.779, 0.0);
    if (!pCase2->getPosition().isEqualTo(pNormal->getPosition() + shift, 1e-10))
    {
        logger(FATAL, "E6 The particle is in the wrong position. It is %, however is should be %",
               pCase2->getPosition(), pNormal->getPosition() + shift);
    }
    if (!pCase2->getVelocity().isEqualTo(pNormal->getVelocity(), 1e-10))
    {
        logger(FATAL, "E7 The particle has the wrong velocity. It is %, however is should be %", pCase2->getVelocity(),
               pNormal->getVelocity());
    }
}
