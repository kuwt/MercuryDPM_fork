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

#include <Particles/LiquidFilmParticle.h>
#include "Mercury3D.h"
#include "Species/LinearViscoelasticFrictionLiquidMigrationWilletSpecies.h"
using constants::pi;
using mathsFunc::cubic;

/*!
 * A two particle pairs collide in the MPI transition zone
 */
class LiquidMigrationMPI2Test : public Mercury3D {
public:

    //sets up name, domain, parallel decomposition, time step, repulsive-adhesive species
    LiquidMigrationMPI2Test() {
        setName("LiquidMigrationMPI2Test");
        setMin(-Vec3D(2.1,1,1));
        setMax(Vec3D(2.1,1,1));
        setNumberOfDomains({NUMBER_OF_PROCESSORS,1,1});
        setTimeStep(1e-4);
        setTimeMax(3e-2);
        setSaveCount(10);
        eneFile.setFileType(FileType::NO_FILE);
        fStatFile.setFileType(FileType::NO_FILE);
        restartFile.writeFirstAndLastTimeStep();
        setXBallsAdditionalArguments("-solidf");

        //define species
        LinearViscoelasticFrictionLiquidMigrationWilletSpecies species;
        species.setDensity(6.0/pi);
        species.setStiffness(2e5);
        species.setDissipation(25);
        species.setContactAngle(0);
        double effectiveRadius = 0.25;
        species.setSurfaceTension(1./(2.0 * constants::pi * effectiveRadius ));
        species.setLiquidBridgeVolumeMax(1e-10);
        speciesHandler.copyAndAddObject(species);
    }

    //sets up two particle pairs at equilibirium overlap at edge of domain, moving towards each other
    void setupInitialConditions() override
    {
        auto s = speciesHandler.getLastObject();

        LiquidFilmParticle particle;
        particle.setSpecies(s);
        particle.setLiquidVolume(1.5e-10);
        particle.setRadius(0.5);
        Mdouble equilibriumOverlap = 5e-6;//species.getAdhesionForceMax()/species.getLoadingStiffness();
        Mdouble reducedRadius = particle.getRadius()-0.5*equilibriumOverlap;
        logger(INFO, "equilibriumOverlap % reducedRadius %", equilibriumOverlap, reducedRadius);

        particle.setVelocity(Vec3D(10,0,0));
        particle.setPosition(Vec3D(getXMin()+reducedRadius,0,0));
        particleHandler.copyAndAddObject(particle);
        particle.setPosition(Vec3D(getXMin()+3*reducedRadius,0,0));
        particleHandler.copyAndAddObject(particle);

        particle.setVelocity(Vec3D(-10,0,0));
        particle.setPosition(Vec3D(getXMax()-3*reducedRadius,0,0));
        particleHandler.copyAndAddObject(particle);
        particle.setPosition(Vec3D(getXMax()-reducedRadius,0,0));
        particleHandler.copyAndAddObject(particle);

//        write(std::cout,false);
    }

    void outputXBallsData(std::ostream& os) const override
    {
        os << particleHandler.getSize()*2.0 + interactionHandler.getSize()
           << " " << getTime()
           << " " << getXMin()
           << " " << getYMin()
           << " " << getZMin()
           << " " << getXMax()
           << " " << getYMax()
           << " " << getZMax()
           << " " << std::endl;

        // This outputs the particle data
        for (unsigned int i = 0; i < particleHandler.getSize(); i++)
        {
            outputXBallsDataParticle(i, 14, os);
        }
        // This outputs the particle data
        for (auto i : particleHandler)
        {
            LiquidFilmParticle* j = dynamic_cast<LiquidFilmParticle*>(i);
            os
                    << j->getPosition().X << " "
                    << j->getPosition().Y << " "
                    << j->getPosition().Z + 1 << " "
                    << 100*(j->getLiquidVolume() > 0) << " 0 0 "
                    << sqrt(j->getLiquidVolume() / 2e-11)*0.5 / 5
                    << " 0 0 0 0 0 0 0" << std::endl;
        }
        // This outputs the particle data
        for (auto i : interactionHandler)
        {
            LiquidMigrationWilletInteraction* j = dynamic_cast<LiquidMigrationWilletInteraction*>(i);
            os
                    << j->getContactPoint().X << " "
                    << j->getContactPoint().Y << " "
                    << j->getContactPoint().Z + 1 << " "
                    << 100*(j->getLiquidBridgeVolume() > 0) << " 0 0 "
                    << sqrt(j->getLiquidBridgeVolume() / 2e-11)*0.5 / 5
                    << " 0 0 0 0 0 0 0" << std::endl;
        }
    }

    //plot temporary values of max overlap and sliding spring as the contact moves across the domain
    void printTime() const override
    {
        double liquidBridgeVolume = interactionHandler.getLiquidBridgeVolume()*1e10;
        double liquidFilmVolume = particleHandler.getLiquidFilmVolume()*1e10;
        std::stringstream f;
        for (auto p : particleHandler)
            f << static_cast<LiquidFilmParticle*>(p)->getLiquidVolume()*1e10 << ' ';
        std::stringstream m;
        for (auto p : particleHandler)
            m << static_cast<LiquidFilmParticle*>(p)->isMPIParticle();
        std::stringstream b;
        for (auto  i : interactionHandler)
            b << dynamic_cast<LiquidMigrationWilletInteraction*>(i)->getLiquidBridgeVolume()*1e10 << ' ';
        logger(INFO,"t % F %(%) B % T %",getNumberOfTimeSteps(),f.str(),m.str(),b.str(),liquidFilmVolume+liquidBridgeVolume);
    }

};

int main()
{
    if (!std::is_base_of<MPILiquidFilmParticle,MPIParticle>()) {
        logger(WARN,"This test can only be ran if MPIParticles are derived from MPILiquidFilmParticle");
        return 0;
    }
    LiquidMigrationMPI2Test dpm;
    dpm.solve();
    return 0;
}
