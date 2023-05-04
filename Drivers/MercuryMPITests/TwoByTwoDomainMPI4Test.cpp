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

#include <thread>
#include "Mercury3D.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
using constants::pi;
using mathsFunc::cubic;

/*!
 * A two particle pairs collide in the MPI transition zone
 */
class TwoByTwoMPIDomainMPI4Test : public Mercury3D {
public:

    //sets up name, domain, parallel decomposition, time step, repulsive-adhesive species
    TwoByTwoMPIDomainMPI4Test() {
        setName("TwoByTwoMPIDomainMPI4Test");
        setMax(Vec3D(3,3,3));
        setMin(-getMax());
        setNumberOfDomains({2,2,1});
        setTimeStep(1e-4);
        setTimeMax(4e-4);
        setSaveCount(1);
        eneFile.setFileType(FileType::NO_FILE);
        fStatFile.setFileType(FileType::NO_FILE);
        restartFile.writeFirstAndLastTimeStep();
        setXBallsAdditionalArguments("-solidf");
        
        //define species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        species->setDensity(6.0/pi);
        double mass  = species->getMassFromRadius(0.5);
        species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(1, 1, 1, mass);
        species->setSlidingFrictionCoefficient(1);
        speciesHandler.copyAndAddObject(species);
    }

    //sets up two particle pairs at equilibirium overlap at edge of domain, moving towards each other
    void setupInitialConditions() override
    {
        // place particle #0 on processor 1, give it a velocity so it moves to processor 3
        // place particle #1 on processor 1 (or 3)
        // place such that #0 is a ghost on processors 0 and 2, but #1 is not
        // in that case, a warning occurs:
        // [Process: 2]: In Domain::processReceivedInteractionData: otherParticle (id 1) is nullptr, the interaction data with pGhost (id 0) is not copied. Two particles possibly moved into domain simultaneously. nt = 1
        // [Process: 0]: In Domain::processReceivedInteractionData: otherParticle (id 1) is nullptr, the interaction data with pGhost (id 0) is not copied. Two particles possibly moved into domain simultaneously. nt = 1
        SphericalParticle particle;
        particle.setSpecies(speciesHandler.getLastObject());
        particle.setRadius(0.5);
        particle.setVelocity(Vec3D(0,100,0));
        particle.setPosition(Vec3D(0.51,-0.025,0));
        particleHandler.copyAndAddObject(particle);
        particle.setVelocity(Vec3D(0,0,0));
        particle.setPosition(Vec3D(1.49,-0.005,0));
        particleHandler.copyAndAddObject(particle);
    }
    
    void actionsAfterTimeStep() override {
        logger(INFO,"#time % #real % #ghost % #contacts %", getNumberOfTimeSteps(),
               particleHandler.getNumberOfRealObjectsLocal(),
               particleHandler.getNumberOfObjects()-particleHandler.getNumberOfRealObjectsLocal(),
               interactionHandler.getNumberOfObjects());
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        for (auto c : interactionHandler) {
            auto i = dynamic_cast<SlidingFrictionInteraction*>(c);
            logger(INFO, "contact % % %", i->getP()->getId(), i->getI()->getId(), i->getTangentialOverlap());
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
};

int main()
{
    TwoByTwoMPIDomainMPI4Test problem;
    problem.solve();
    for (auto c : problem.interactionHandler) {
        auto i = dynamic_cast<SlidingFrictionInteraction*>(c);
        helpers::check(i->getTangentialOverlap(), -0.0499896, 1e-6, "Tangential spring");
    }
    return 0;
}
