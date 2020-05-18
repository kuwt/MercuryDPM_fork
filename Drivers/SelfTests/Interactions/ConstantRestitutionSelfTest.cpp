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
#include "Species/LinearPlasticViscoelasticSpecies.h"


struct ConstantRestitutionSelfTest : public Mercury3D {
    void actionsBeforeTimeStep() override {
        if (interactionHandler.getSize()!=0) {
            collisionTime = getTime();
            maxOverlap = std::max(maxOverlap,interactionHandler.getLastObject()->getOverlap());
        }
    }
    Mdouble collisionTime = 0;
    Mdouble maxOverlap = 0;
};

int main() 
{
    ConstantRestitutionSelfTest dpm;
    dpm.setName("ConstantRestitutionSelfTest");
    dpm.setSaveCount(2);

    auto species = dpm.speciesHandler.copyAndAddObject(LinearPlasticViscoelasticSpecies());
    species->setDensity(1000);
    species->setConstantRestitution(true);
    Mdouble radius = 2.5e-3;
    const Mdouble stiffness = 24067/species->getMassFromRadius(radius);
    //species->setPlasticParameters(0.2*stiffness, stiffness, 0.873*stiffness, 0.05);
    species->setPlasticParameters(stiffness, stiffness, 0.0, 0.05);
    species->setRestitutionCoefficient(0.45, 1.0);
    dpm.setTimeStep(0.02*species->getCollisionTime(1.0));
    dpm.setTimeMax(60.0*dpm.getTimeStep());

    radius = 1e-4;
    Mdouble relativeVelocity = 1;

    SphericalParticle particle;
    particle.setSpecies(species);
    particle.setRadius(radius);
    particle.setPosition(Vec3D(-radius,0,0));
    particle.setVelocity(Vec3D(0.5*relativeVelocity,0,0));
    auto particle1 = dpm.particleHandler.copyAndAddObject(particle);

    particle.setPosition(-particle.getPosition());
    particle.setVelocity(-particle.getVelocity());
    auto particle2 = dpm.particleHandler.copyAndAddObject(particle);

    dpm.setDomain(radius*Vec3D(-2,-1,-1),radius*Vec3D(2,1,1));

    dpm.solve();

    //Values that are now independent of the particle radius
    logger(INFO,"Max. Overlap %",dpm.maxOverlap);
    logger(INFO,"Collision time %",dpm.collisionTime);
    logger(INFO,"Restitution %",particle2->getVelocity().X/relativeVelocity);
    return 0;
}

