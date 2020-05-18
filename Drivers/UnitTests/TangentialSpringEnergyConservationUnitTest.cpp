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

#include <iostream>
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "DPMBase.h"
#include "Walls/InfiniteWall.h"

class TangentialSpringEnergyConservationUnitTest : public DPMBase
{
public:
    TangentialSpringEnergyConservationUnitTest()
    {
        //set the species properties, with a very large sliding friction coefficient
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        species->setDensity(6.0 / constants::pi);
        species->setCollisionTimeAndRestitutionCoefficient(1e-3, .2, 1);
        species->setSlidingStiffness(species->getStiffness() * 2 / 7);
        species->setSlidingDissipation(0);
        species->setSlidingFrictionCoefficient(0.5);
        setName("TangentialSpringEnergyConservationUnitTest");
        setMin({-1, -1, -1});
        setMax({1, 1, 1});
        
        setTimeMax(10);
        setTimeStep(species->getCollisionTime(1) / 50);
        setParticlesWriteVTK(true);
        setWallsWriteVTK(FileType::ONE_FILE);
        setSaveCount(0.1 / getTimeStep());
    }
    
    void setupInitialConditions() override
    {
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setPosition({0, 0, 0});
        p.setVelocity({0, 0, 0});
        p.setAngularVelocity({0.1, 0, 0});
        p.setRadius(0.5);
        particleHandler.copyAndAddObject(p);
        
        p.setPosition({0.9999, 0, 0});
        p.fixParticle();
        particleHandler.copyAndAddObject(p);
        p.setPosition({-.9999,0,0});
        particleHandler.copyAndAddObject(p);
        p.setPosition({0,.9999,0});
        particleHandler.copyAndAddObject(p);
        p.setPosition({0,-.9999,0});
        particleHandler.copyAndAddObject(p);
        p.setPosition({0,0,.9999});
        particleHandler.copyAndAddObject(p);
        p.setPosition({0,0,-.9999});
        particleHandler.copyAndAddObject(p);
        
    }
    
    void printTime() const override
    {
        DPMBase::printTime();
        std::cout << " total energy: " << std::setprecision(10) << getRotationalEnergy() + getElasticEnergy() << '\n' << '\n';
    }
};

int main()
{
    TangentialSpringEnergyConservationUnitTest test;
    test.solve();
    return 0;
}
