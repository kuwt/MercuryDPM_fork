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

#include<iostream>
#include "Mercury3D.h"
#include <Species/LinearViscoelasticFrictionSpecies.h>
#include <Walls/InfiniteWall.h>
#include <Boundaries/PeriodicBoundary.h>

class MercuryLogo : public Mercury3D
{
public:
    void setupInitialConditions() override {
        logger(INFO, "Setting up a simulation around the Mercury logo");

        const double rho = constants::pi;
        const double restitutionCoefficient = 0.05;

        const double mass_1 = 4 / 3 * constants::pi * pow(0.5, 3.0) * rho;
        const double mass_2 = 4 / 3 * constants::pi * pow(0.35, 3.0) * rho;

        auto species0 = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        auto species1 = speciesHandler.copyAndAddObject(species0);
        auto species01 = speciesHandler.getMixedObject(species0, species1);

        species0->setDensity(rho);
        species0->setCollisionTimeAndRestitutionCoefficient(tc, restitutionCoefficient, mass_1);
        //  Set the tangential dissipation equal to the normal dissipation for small-small collisions
        species0->setSlidingDissipation(species0->getDissipation());
        species0->setSlidingStiffness(species0->getStiffness() * 2.0 / 7.0);
        species0->setSlidingFrictionCoefficient(0.5);
        //////
        //
        species1->setDensity(rho);
        species1->setCollisionTimeAndRestitutionCoefficient(tc, restitutionCoefficient, mass_2);
        //  Set the tangential dissipation equal to the normal dissipation for large-large collision
        species1->setSlidingDissipation(species1->getDissipation());
        species1->setSlidingStiffness(species1->getStiffness() * 2.0 / 7.0);
        species1->setSlidingFrictionCoefficient(0.5);
        //
        species01->setCollisionTimeAndRestitutionCoefficient(tc, restitutionCoefficient, mass_1, mass_2);
        //  Set the tangential dissipation equal to the normal dissipation for mixed collision
        species01->setSlidingDissipation(species01->getDissipation());
        species01->setSlidingFrictionCoefficient(0.5);
        species01->setSlidingStiffness(species01->getStiffness() * 2.0 / 7.0);
        //
        SphericalParticle p0;
        p0.setRadius(0.5);
        p0.fixParticle();
        p0.setSpecies(species0);
        ////////////////
        constructTextAsParticles(p0);
        ////////////////
        //Set the walls of the domain
        InfiniteWall w0;
        w0.setSpecies(species1);
        //
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0, getYMin(), 0.0));
        wallHandler.copyAndAddObject(w0);

        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0.0, getYMax(), 0.0));
        wallHandler.copyAndAddObject(w0);
        ///////////////////
        //Make the "free" particle to be inserted
        SphericalParticle p1;
        p1.setRadius(0.35);
        p1.setSpecies(species1);
        Vec3D pos;
        ////////////////////
        //Insert the moving particles
        for (unsigned int numberOfParticlesInserted = 0;
             numberOfParticlesInserted < numberOfParticlesToBeInserted;
             ++numberOfParticlesInserted)
        {
            pos.X = random.getRandomNumber(getXMin(), getXMax());
            pos.Y = random.getRandomNumber(getYMin(), getYMax());
            pos.Z = random.getRandomNumber(getZMin(), getZMax());
            p1.setPosition(pos);
            p1.setVelocity(Vec3D(random.getRandomNumber(0.7, 1.2),
                                 random.getRandomNumber(-0.2, 0.2),
                                 random.getRandomNumber(0.3, 1.4)));

            particleHandler.copyAndAddObject(p1);
        }
        ///////////////////
        PeriodicBoundary pb0;
        pb0.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(pb0);
        pb0.set(Vec3D(0.0, 0.0, 1.0), getZMin(), getZMax());
        boundaryHandler.copyAndAddObject(pb0);
        logger(INFO, "Setup finished");
    }

    ///Write the words "Mercury=fun" in the domain, where each symbol is constructed by some particles. The input particle
    /// p0 is the particle that is copied and inserted in the domain every time, it needs to be a fixed particle with a
    /// species and radius assigned to it.
    void constructTextAsParticles(BaseParticle& p0)
    {
        const double pDia = 2.0 * p0.getRadius();
        ///////////////////////
        // Letter M
        /////////////////////
        const double xPosLine1 = getXMin() + pDia;
        const double xPosLine2 = getXMin() + 5.0 * pDia;//x-location at the
        // 2nd edge of letter M, which is 3*partDia apart from 1st edge.
        for (unsigned int i = 0; i < 4; i++)
        {
            p0.setPosition(Vec3D(xPosLine1, 0.5 * getYMax(), getZMax() - (i + 1) * pDia));
            particleHandler.copyAndAddObject(p0);
            p0.setPosition(Vec3D(xPosLine2, 0.5 * getYMax(), getZMax() - (i + 1) * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        double finalZPos = getZMax() - 4.0 * pDia;
        // Being simple
        p0.setPosition(Vec3D(xPosLine1 + pDia, 0.5 * getYMax(), getZMax() - 1.5 * pDia));
        particleHandler.copyAndAddObject(p0);
        p0.setPosition(Vec3D(xPosLine1 + 2.0 * pDia, 0.5 * getYMax(), getZMax() - 2.0 * pDia));
        particleHandler.copyAndAddObject(p0);
        p0.setPosition(Vec3D(xPosLine2 - pDia, 0.5 * getYMax(), getZMax() - 1.5 * pDia));
        particleHandler.copyAndAddObject(p0);
        ////////////////////
        // Letter E
        ///////////////////
        // xPosLine2 is from letter 'M'
        // finalZPos is the z-location of the bottom most tip of from letter 'M'
        for (unsigned int i = 0; i < 5; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + 2.0 * pDia, 0.5 * getYMax(), finalZPos - (i + 1) * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        //
        for (unsigned int i = 0; i < 3; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + (i + 3) * pDia, 0.5 * getYMax(), finalZPos - pDia));
            particleHandler.copyAndAddObject(p0);
            p0.setPosition(Vec3D(xPosLine2 + (i + 3) * pDia, 0.5 * getYMax(), finalZPos - 5.0 * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        //
        for (unsigned int i = 0; i < 2; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + (i + 3) * pDia, 0.5 * getYMax(), finalZPos - 3.0 * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        finalZPos = finalZPos - 6.0 * pDia;// update final Z-position
        //////////////////
        // Letter R
        /////////////////
        // finalZPos is the z-location of the bottom most tip of from letter 'E'
        for (unsigned int i = 0; i < 5; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + 7.0 * pDia, 0.5 * getYMax(), finalZPos - i * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        //
        for (unsigned int i = 0; i < 3; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + (i + 8) * pDia, 0.5 * getYMax(), finalZPos));
            particleHandler.copyAndAddObject(p0);
            p0.setPosition(Vec3D(xPosLine2 + (i + 8) * pDia, 0.5 * getYMax(), finalZPos - 2.0 * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        // Again being simple and placing the dots manually.
        p0.setPosition(Vec3D(xPosLine2 + 10.0 * pDia, 0.5 * getYMax(), finalZPos - pDia));
        particleHandler.copyAndAddObject(p0);
        //
        p0.setPosition(Vec3D(xPosLine2 + 8.5 * pDia, 0.5 * getYMax(), finalZPos - 3.0 * pDia));
        particleHandler.copyAndAddObject(p0);
        //
        p0.setPosition(Vec3D(xPosLine2 + 9.5 * pDia, 0.5 * getYMax(), finalZPos - 4.0 * pDia));
        particleHandler.copyAndAddObject(p0);
        //
        finalZPos = finalZPos - 5.0 * pDia;
        ////////////////
        // Letter C
        ///////////////
        for (unsigned int i = 0; i < 5; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + 11.0 * pDia, 0.5 * getYMax(), finalZPos - i * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        //
        for (unsigned int i = 0; i < 3; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + (i + 12) * pDia, 0.5 * getYMax(), finalZPos));
            particleHandler.copyAndAddObject(p0);
            p0.setPosition(Vec3D(xPosLine2 + (i + 12) * pDia, 0.5 * getYMax(), finalZPos - 4.0 * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        finalZPos = finalZPos - 5.0 * pDia;
        ////////////////
        // Letter U
        //////////////
        for (unsigned int i = 0; i < 5; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + 6.0 * pDia, 0.5 * getYMax(), finalZPos - i * pDia));
            particleHandler.copyAndAddObject(p0);
            p0.setPosition(Vec3D(xPosLine2 + 9.0 * pDia, 0.5 * getYMax(), finalZPos - i * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        for (unsigned int i = 0; i < 2; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + (i + 7.0) * pDia, 0.5 * getYMax(), finalZPos - 4.0 * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        finalZPos = finalZPos - 5.0 * pDia;
        ///////////////
        // Letter R
        //////////////
        for (unsigned int i = 0; i < 5; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + pDia, 0.5 * getYMax(), finalZPos - i * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        //
        for (unsigned int i = 0; i < 3; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + (i + 2) * pDia, 0.5 * getYMax(), finalZPos));
            particleHandler.copyAndAddObject(p0);
            p0.setPosition(Vec3D(xPosLine2 + (i + 2) * pDia, 0.5 * getYMax(), finalZPos - 2.0 * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        // Again being simple and placing the dots manually.
        p0.setPosition(Vec3D(xPosLine2 + 4.0 * pDia, 0.5 * getYMax(), finalZPos - pDia));
        particleHandler.copyAndAddObject(p0);
        //
        p0.setPosition(Vec3D(xPosLine2 + 2.5 * pDia, 0.5 * getYMax(), finalZPos - 3.0 * pDia));
        particleHandler.copyAndAddObject(p0);
        //
        p0.setPosition(Vec3D(xPosLine2 + 3.5 * pDia, 0.5 * getYMax(), finalZPos - 4.0 * pDia));
        particleHandler.copyAndAddObject(p0);
        //
        finalZPos = finalZPos - 5.0 * pDia;
        ////////////
        // Letter Y
        ////////////
        for (unsigned int i = 0; i < 5; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 - 0.8 * (i + 1) * pDia, 0.5 * getYMax(), finalZPos - i * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        //
        for (unsigned int i = 0; i < 2; i++)
        {
            p0.setPosition(Vec3D(xPosLine1 + i * 0.8 * pDia, 0.5 * getYMax(), finalZPos - i * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        finalZPos = finalZPos - 5.0 * pDia;
        /////////////
        // symbol '='
        //////////////
        logger(DEBUG, "%", xPosLine2);
        for (unsigned int i = 0; i < 5; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + (i + 5) * pDia, 0.5 * getYMax(), finalZPos + 2.0 * pDia));
            particleHandler.copyAndAddObject(p0);
            p0.setPosition(Vec3D(xPosLine2 + (i + 5) * pDia, 0.5 * getYMax(), finalZPos));
            particleHandler.copyAndAddObject(p0);
        }
        finalZPos = finalZPos - 5.0 * pDia;
        ///////////////
        // letter 'F'
        ///////////////
        for (unsigned int i = 0; i < 5; i++)
        {
            p0.setPosition(Vec3D(xPosLine2, 0.5 * getYMax(), finalZPos - i * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        for (unsigned int i = 0; i < 3; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + (i + 1) * pDia, 0.5 * getYMax(), finalZPos));
            particleHandler.copyAndAddObject(p0);
        }
        for (unsigned int i = 0; i < 2; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + (i + 1) * pDia, 0.5 * getYMax(), finalZPos - 2.0 * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        ///////////////
        // letter 'U'
        ///////////////
        for (unsigned int i = 0; i < 5; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + 5.0 * pDia, 0.5 * getYMax(), finalZPos - i * pDia));
            particleHandler.copyAndAddObject(p0);
            p0.setPosition(Vec3D(xPosLine2 + 8.0 * pDia, 0.5 * getYMax(), finalZPos - i * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        for (unsigned int i = 0; i < 2; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + (i + 6.0) * pDia, 0.5 * getYMax(), finalZPos - 4.0 * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        //////////////
        // letter 'N'
        //////////////
        for (unsigned int i = 0; i < 5; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + 10.0 * pDia, 0.5 * getYMax(), finalZPos - i * pDia));
            particleHandler.copyAndAddObject(p0);
            p0.setPosition(Vec3D(xPosLine2 + 14.0 * pDia, 0.5 * getYMax(), finalZPos - i * pDia));
            particleHandler.copyAndAddObject(p0);
        }
        for (unsigned int i = 0; i < 5; i++)
        {
            p0.setPosition(Vec3D(xPosLine2 + (i + 19) * 0.58 * pDia, 0.5 * getYMax(), finalZPos - i * pDia));
            particleHandler.copyAndAddObject(p0);
        }
    }

public:
    unsigned int numberOfParticlesToBeInserted;
    const double tc = 0.005;
};

int main()
{
    MercuryLogo logo;
    logo.numberOfParticlesToBeInserted = 600;
    logo.setName("MercuryLogo");
    logo.setZMax(45.0);
    logo.setYMax(5.0);
    logo.setXMax(20.0);
    logo.setSystemDimensions(3);
    logo.setTimeStep(logo.tc * 0.02);
    logo.setTimeMax(1.0);
    logo.setGravity(Vec3D(0.0, 0.0, 0.0));

    logo.setSaveCount(100);
    logo.setFileType(FileType::ONE_FILE);
    logo.setXBallsAdditionalArguments("-solidf -v0");
    logo.solve();
    return 0;
}
