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
#include "Species/HertzianViscoelasticFrictionSpecies.h"

class Hertzian2DUnitTest : public Mercury2D
{
    public:

        Hertzian2DUnitTest()
        {

            setName("Hertzian2DSelfTest");
            setDomain({-1,-1,-1},{1,1,1});
            setTimeMax(4);
            setTimeStep(5e-5);
            setSaveCount(200);

            auto spec = new HertzianViscoelasticFrictionSpecies();
            spec->setDensity(1);
            spec->setEffectiveElasticModulusAndRestitutionCoefficient(6e2, 0.80);
            spec = speciesHandler.copyAndAddObject(spec);

            /* Collision between fixed and movable */
            auto pf = new SphericalParticle;
            pf->setSpecies(spec);
            pf->setRadius(0.1);
            pf->setPosition(Vec3D(0,0,0));
            pf->fixParticle();
            particleHandler.copyAndAddObject(pf);

            auto pm = new SphericalParticle;
            pm->setSpecies(spec);
            pm->setRadius(0.1);
            pm->setPosition(Vec3D(1,0,0));
            pm->setVelocity(Vec3D(-1,0,0));
            particleHandler.copyAndAddObject(pm);

            /* Collision between two movables */
            pm->setPosition(Vec3D(0,1,0));
            pm->setVelocity(Vec3D(1,0,0));
            particleHandler.copyAndAddObject(pm);
            
            pm->setPosition(Vec3D(1,1,0));
            pm->setVelocity(Vec3D(0,0,0));
            particleHandler.copyAndAddObject(pm);


        }

        ~Hertzian2DUnitTest() override = default;
};

int main() 
{
    auto problem = new Hertzian2DUnitTest;
    problem->solve();
    return 0;
}

