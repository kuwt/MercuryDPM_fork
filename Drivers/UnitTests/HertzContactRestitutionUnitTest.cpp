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
#include "Math/ExtendedMath.h"
#include "Species/HertzianViscoelasticMindlinSpecies.h"
class HertzContactRestitutionUnitTest : public DPMBase
{
public:
    void setupInitialConditions() override
    {
        HertzianViscoelasticMindlinSpecies species;
        species.setEffectiveElasticModulusAndRestitutionCoefficient(20000, 0.8);
        speciesHandler.copyAndAddObject(species);
        
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        SphericalParticle p1 = *(p0.copy());
        
        p0.setPosition(Vec3D(0, 0.0, 0.0));
        p1.setPosition(Vec3D(3, 0.0, 0.0));
        p0.setVelocity(Vec3D(1, 0.0, 0.0));
        p1.setVelocity(Vec3D(-1, 0.0, 0.0));
        particleHandler.copyAndAddObject(p0);
        particleHandler.copyAndAddObject(p1);
        
        setTimeStep(species.getCollisionTime(1, 1, constants::pi / 6) / 50);
        setTimeMax(1.5);
        setMin(-1, -1, -1);
        setMax(4, 1, 1);
    }
    
    void test()
    {
        setName("HertzContactRestitutionUnitTest");
        setSaveCount(500);
        solve();
        logger.assert_always(mathsFunc::isEqual(particleHandler.getObject(1)->getVelocity().X, 0.8, 1e-2),
                             "Particle 1 has wrong end velocity: % (should be 0.8) ",
                             particleHandler.getObject(1)->getVelocity().X);
        logger.assert_always(mathsFunc::isEqual(particleHandler.getObject(0)->getVelocity().X, -0.8, 1e-2),
                             "Particle 0 has wrong end velocity: % (should be -0.8) ",
                             particleHandler.getObject(0)->getVelocity().X);

    }
};

int main(int argc, char* argv[])
{
    HertzContactRestitutionUnitTest problem;
    problem.test();
    return 0;
}
