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

#include "Mercury2D.h"
#include "Particles/LiquidFilmParticle.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include "Species/Species.h"
#include "Species/LinearViscoelasticFrictionLiquidMigrationWilletSpecies.h"
#include <iomanip>

/// In this file two particles are symmetrically placed in a bi-axial box are allowed to jump around under gravity. It tests walls gravity and symmetry.

class LiquidMigrationSelfTest : public Mercury2D
{

    void setupInitialConditions() override {
        LiquidFilmParticle P0, P1, P2, P3, P4;
        P0.setSpecies(speciesHandler.getObject(0));
        P1.setSpecies(speciesHandler.getObject(0));
        P2.setSpecies(speciesHandler.getObject(0));
        P3.setSpecies(speciesHandler.getObject(0));
        P4.setSpecies(speciesHandler.getObject(0));
        
        P0.setLiquidVolume(2e-12);
        P1.setLiquidVolume(2e-12);
        P2.setLiquidVolume(2e-12);
        P3.setLiquidVolume(2e-12);
        P4.setLiquidVolume(2e-12);

        P0.setPosition(Vec3D(0.006, 0.005, 0.0));
        P1.setPosition(Vec3D(0.007, 0.005, 0.0));
        P2.setPosition(Vec3D(0.005, 0.005, 0.0));
        P3.setPosition(Vec3D(0.004, 0.005, 0.0));
        P4.setPosition(Vec3D(0.004, 0.006, 0.0));

        P0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        P1.setVelocity(Vec3D(0.0, 0.0, 0.0));
        P2.setVelocity(Vec3D(0.0, 0.0, 0.0));
        P3.setVelocity(Vec3D(0.03, 0.0, 0.0));
        P4.setVelocity(Vec3D(0.0, 0.0, 0.0));

        P0.setRadius(0.0005);
        P1.setRadius(0.0005);
        P2.setRadius(0.0005);
        P3.setRadius(0.0005);
        P4.setRadius(0.0005);
        particleHandler.copyAndAddObject(P0);
        particleHandler.copyAndAddObject(P1);
        particleHandler.copyAndAddObject(P2);
        particleHandler.copyAndAddObject(P3);
        particleHandler.copyAndAddObject(P4);

        setGravity(Vec3D(0, 0, 0));
        setMax({0.01,0.01,0.01});

        wallHandler.clear();
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(-1, 0, 0), Vec3D(getXMin(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(1, 0, 0), Vec3D(getXMax(), 0, 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, -1, 0), Vec3D(0, getYMin(), 0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0, 1, 0), Vec3D(0, getYMax(), 0));
        wallHandler.copyAndAddObject(w0);
    }

    void outputXBallsData(std::ostream& os) const override {
        os << particleHandler.getNumberOfObjects() * 2.0 + interactionHandler.getNumberOfObjects()
           << " " << getTime()
           << " " << getXMin()
           << " " << getYMin()
           << " " << getZMin()
           << " " << getXMax()
           << " " << getYMax()
           << " " << getZMax()
           << " " << std::endl;

        // This outputs the particle data
        for (unsigned int i = 0; i < particleHandler.getNumberOfObjects(); i++)
        {
            outputXBallsDataParticle(i, 14, os);
        }
        // This outputs the particle data
        for (BaseParticle * const p : particleHandler)
        {
            LiquidFilmParticle* j = dynamic_cast<LiquidFilmParticle*>(p);
            os << j->getPosition().X << " "
               << j->getPosition().Y << " "
               << j->getPosition().Z + 1 << " "
               << (j->getLiquidVolume() > 0) << " 0 0 "
               << sqrt(j->getLiquidVolume() / 2e-12) * 0.0005 / 5
               << " 0 0 0 0 0 0 0" << std::endl;
        }
        // This outputs the particle data
        for (BaseInteraction * const i : interactionHandler)
        {
            LiquidMigrationWilletInteraction* j = dynamic_cast<LiquidMigrationWilletInteraction*>(i);
            os << j->getContactPoint().X << " "
               << j->getContactPoint().Y << " "
               << j->getContactPoint().Z + 1 << " "
               << (j->getLiquidBridgeVolume() > 0) << " 0 0 "
               << sqrt(j->getLiquidBridgeVolume() / 2e-12) * 0.0005 / 5
               << " 0 0 0 0 0 0 0" << std::endl;
        }
        logger(DEBUG, "Have output the properties of the problem to disk ");
    }

    void printTime() const override {
        Mdouble volTotP = 0.0;
        unsigned int nLB = 0;
        for (BaseParticle* const p : particleHandler)
        {
            const LiquidFilmParticle* l = dynamic_cast<LiquidFilmParticle*>(p);
            if (l != nullptr)
            {
                volTotP += l->getLiquidVolume();
            }
            else
            {
                logger(ERROR, "The pointer to the LiquidFilmParticle is a null pointer");
            }
        }
        Mdouble volTotI = 0.0;
        for (BaseInteraction* const i : interactionHandler)
        {
            const LiquidMigrationWilletInteraction* l = dynamic_cast<LiquidMigrationWilletInteraction*>(i);
            if (l != nullptr)
            {
                volTotI += l->getLiquidBridgeVolume();
                if (l->getLiquidBridgeVolume() != 0.0)
                {
                    nLB++;
                }
            }
            else
            {
                logger(ERROR, "The pointer to the LiquidMigrationWilletInteraction is a null pointer");
            }

        }
        logger(INFO, "t = %", getTime());
        logger(INFO, "Vol = % (% in particle), #LB = %", volTotP + volTotI, volTotP, nLB);

        ///\todo how to set the precision in the logger? Below the former code.
        /*
        std::cout
        << "t=" << std::setprecision(6) << std::left << std::setw(6) << getTime()
        << ", Vol=" << std::setprecision(8) << std::left << std::setw(6) << volTotP + volTotI
        << "(" << std::setprecision(6) << std::left << std::setw(6) << volTotP
        << " in particle), #LB=" << std::setprecision(6) << std::left << std::setw(12) << nLB
        << std::endl;
        std::cout.flush();
         */
    }

    
};

int main(int argc UNUSED, char* argv[] UNUSED)
{
    LiquidMigrationSelfTest LiquidMigrationProblem;
    LiquidMigrationProblem.setName("LiquidMigrationSelfTest");
    LiquidMigrationProblem.setSystemDimensions(3);
    LiquidMigrationProblem.setXBallsAdditionalArguments("-v0 -solid -3dturn 1");
    auto species = new LinearViscoelasticFrictionLiquidMigrationWilletSpecies;
    LiquidMigrationProblem.speciesHandler.addObject(species);
    species->setDensity(2000);
    species->setStiffness(10000);
    species->setLiquidBridgeVolumeMax(5e-12);
    species->setDistributionCoefficient(0.8);
    species->setSurfaceTension(1.0);
    species->setContactAngle(20.0 * constants::pi / 180.0);

    LiquidMigrationProblem.setTimeMax(0.7);//2.5);
    LiquidMigrationProblem.setSaveCount(100);
    LiquidMigrationProblem.setTimeStep(2e-5);
    LiquidMigrationProblem.fStatFile.setFileType(FileType::NO_FILE);
    LiquidMigrationProblem.solve();
}
