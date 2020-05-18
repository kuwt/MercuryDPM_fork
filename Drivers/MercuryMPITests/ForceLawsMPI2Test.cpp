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

#include "Mercury3D.h"
#include "Species/LinearPlasticViscoelasticFrictionIrreversibleAdhesiveSpecies.h"
using constants::pi;
using mathsFunc::cubic;

class ForceLawsMPI2Test : public Mercury3D {
public:

    //sets up name, domain, parallel decomposition, time step, repulsive-adhesive species
    ForceLawsMPI2Test() {
        setName("ForceLawsMPI2Test");
        setMin(-Vec3D(5,1,1));
        setMax(Vec3D(5,1,1));
        setNumberOfDomains({NUMBER_OF_PROCESSORS,1,1});
        setTimeStep(1e-4);
        setTimeMax(1.0);
        setSaveCount((unsigned)(getTimeMax()/getTimeStep()/20.));
        //fStatFile.writeFirstAndLastTimeStep();
        eneFile.setFileType(FileType::NO_FILE);
        fStatFile.setFileType(FileType::NO_FILE);
        restartFile.writeFirstAndLastTimeStep();
        setXBallsAdditionalArguments("-solidf");

        LinearPlasticViscoelasticFrictionIrreversibleAdhesiveSpecies species;
        species.setDensity(6.0/pi);
        species.setPlasticParameters(2e5,10e5,2e5,.5);
        species.setDissipation(25);
        species.setSlidingFrictionCoefficient(0.5);
        species.setSlidingDissipation(2./7.*species.getDissipation());
        species.setSlidingStiffness(2./7.*species.getLoadingStiffness());
        species.setAdhesionForceMax(1);
        species.setAdhesionStiffness(species.getLoadingStiffness());
        speciesHandler.copyAndAddObject(species);



    }

    //sets up two particles with small overlap (such that forces are in equilibrium) at edge of domain, moving to the right
    void setupInitialConditions() override
    {
        auto s = speciesHandler.getLastObject();

        SphericalParticle particle;
        particle.setSpecies(s);
        particle.setRadius(0.5);
        particle.setVelocity(Vec3D(10,0,0));
        Mdouble equilibriumOverlap = 5e-6;//species.getAdhesionForceMax()/species.getLoadingStiffness();
        Mdouble reducedRadius = particle.getRadius()-0.5*equilibriumOverlap;
        logger(INFO, "equilibriumOverlap % reducedRadius %", equilibriumOverlap, reducedRadius);
        particle.setPosition(Vec3D(getXMin()+reducedRadius,0,0));
        particleHandler.copyAndAddObject(particle);

        particle.setPosition(Vec3D(getXMin()-reducedRadius,0,0));
        particleHandler.copyAndAddObject(particle);
        write(std::cout,false);
    }

    //set initial values of max overlap and sliding spring
    void actionsBeforeTimeStep() override
    {
        if (getNumberOfTimeSteps()==0) {
            //logger.assert_always(interactionHandler.getSize() != 0, "interaction not found");
            if (interactionHandler.getSize() != 0) {
                auto i = dynamic_cast<LinearPlasticViscoelasticInteraction *>(interactionHandler.getLastObject());
                logger.assert_always(i, "Interaction type needs to be LinearPlasticViscoelastic");
                i->setMaxOverlap(1e-5);
                logger(INFO,"changed max overlap to %",i->getMaxOverlap());

                auto j = dynamic_cast<FrictionInteraction *>(interactionHandler.getLastObject());
                logger.assert_always(j, "Interaction type needs to be Friction");
                j->setSlidingSpring({1e-5,0,0});
                logger(INFO,"changed sliding spring to %",j->getSlidingSpring().X);

            }
        }
   }

    //plot temporary values of max overlap and sliding spring as the contact moves across the domain
    void printTime() const override
    {
        if (interactionHandler.getSize()!=0) {
            auto i = dynamic_cast<const LinearPlasticViscoelasticInteraction *>(interactionHandler.getLastObject());
            logger.assert_always(i, "Interaction type needs to be LinearPlasticViscoelastic");
            auto j = dynamic_cast<const FrictionInteraction *>(interactionHandler.getLastObject());
            logger.assert_always(j, "Interaction type needs to be Friction");

            logger(INFO, "t % #inter % domain % middle % maxOverlap % overlap % sliding %", getTime(), interactionHandler.getSize(), PROCESSOR_ID, domainHandler.getCurrentDomain()->getMiddle().X, i->getMaxOverlap(), i->getOverlap(),j->getSlidingSpring().X);
        } else {
            logger(INFO, "t % domain % middle %", getTime(), PROCESSOR_ID, domainHandler.getCurrentDomain()->getMiddle().X);
        }
    }

    //check final values of max overlap and sliding spring after the contact moved from processor 0 to processor 1
    void actionsAfterSolve() override
    {
        if (PROCESSOR_ID!=1) return;
        logger.assert(interactionHandler.getSize()>0, "Contact was lost");
        auto i = dynamic_cast<const LinearPlasticViscoelasticInteraction *>(interactionHandler.getLastObject());
        logger.assert_always(i, "Interaction type needs to be LinearPlasticViscoelastic");
        auto j = dynamic_cast<const FrictionInteraction *>(interactionHandler.getLastObject());
        logger.assert_always(j, "Interaction type needs to be Friction");

        helpers::check(i->getMaxOverlap(),1e-5,1e-10,"Checking whether max overlap was preserved by domain transition");
        helpers::check(j->getSlidingSpring().X,8.7489e-06,1e-10,"Checking whether sliding spring was preserved by domain transition");
    }

};

int main()
{
    ForceLawsMPI2Test dpm;
    dpm.solve();
    //return 0;
}
