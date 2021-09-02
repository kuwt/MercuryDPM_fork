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

#include <Species/HertzianViscoelasticMindlinSpecies.h>
#include "Mercury3D.h"
#include "Particles/SuperQuadricParticle.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Math/ExtendedMath.h"

class SlidingSpheresUnitTest: public Mercury3D
{
public:
    void setupInitialConditions() override
    {
        //from FreeFallHertzMindlinUnitTest
        Mdouble poissonRatio = 0.28;
        Mdouble density = 6.0/constants::pi/100;
        Mdouble elasticModulus = 4.11e9;
        Mdouble restitution = 0.8;
        
        HertzianViscoelasticMindlinSpecies species;
        species.setDensity(density);
        species.setEffectiveElasticModulusAndRestitutionCoefficient(elasticModulus, restitution);
        //https://en.wikipedia.org/wiki/Shear_modulus#References
        species.setEffectiveShearModulus(0.5 * elasticModulus / (1 + poissonRatio));
        species.setSlidingFrictionCoefficient(1.0);
        speciesHandler.copyAndAddObject(species);
        
        SuperQuadricParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setAxesAndExponents(1.0,1.0,1.0,1.0,1.0);
        p0.setInertia();
        SuperQuadricParticle p1 = *(p0.copy());
        p0.setPosition(Vec3D(1.2, 0.0, 0.0));
        p1.setPosition(Vec3D(2.9, 0.0, 2.0));
        p0.setVelocity(Vec3D(0.0, 0.0, 1.0));
        p1.setVelocity(Vec3D(0.0, 0.0, -1.0));
        particleHandler.copyAndAddObject(p0);
        particleHandler.copyAndAddObject(p1);
        logger.assert_always(mathsFunc::isEqual(p0.getMaxInteractionRadius(), 1.0, 1e-10),
                             "interaction radius p0 equals % but should be 1.0", p0.getMaxInteractionRadius());
        
        logger.assert_always(mathsFunc::isEqual(p1.getMaxInteractionRadius(), 1.0, 1e-10),
                             "interaction radius p1 equals % but should be 1.0", p1.getMaxInteractionRadius());
        
        
        SphericalParticle pSphere;
        pSphere.setSpecies(speciesHandler.getObject(0));
        pSphere.setRadius(1);
        pSphere.setInertia();
        SphericalParticle pSphere1 = pSphere;
        pSphere.setPosition(Vec3D(11.2, 0.0, 0.0));
        pSphere1.setPosition(Vec3D(12.9, 0.0, 2.0));
        pSphere.setVelocity(Vec3D(0.0, 0.0, 1.0));
        pSphere1.setVelocity(Vec3D(0.0, 0.0, -1.0));
        particleHandler.copyAndAddObject(pSphere);
        particleHandler.copyAndAddObject(pSphere1);
        
        
        setTimeStep(species.getCollisionTime(1, 1, density) / 50);
        setTimeMax(2);
        setMin(-10, -10, -10);
        setMax(20, 10, 10);
    }
    
    //check contact time
    void actionsAfterTimeStep() override
    {
        ///\todo should getTimeStep be getNTimeStep?
        BaseInteraction* interaction = particleHandler.getObject(0)->
                getInteractionWith(particleHandler.getObject(1), getTimeStep(), &interactionHandler);
        if (interaction!= nullptr) contactHasOccured = true;
    }
    
    void test()
    {
        setName("SlidingSpheresUnitTest");
        setSaveCount(5000);
        solve();
        logger.assert_always(contactHasOccured, "The contact has not occurred");
    }

private:
    bool contactHasOccured;
};

int main(int argc, char* argv[])
{
    SlidingSpheresUnitTest problem;
    // comment next line to turn on file output
    problem.setFileType(FileType::NO_FILE);
    problem.test();
    return 0;
}
