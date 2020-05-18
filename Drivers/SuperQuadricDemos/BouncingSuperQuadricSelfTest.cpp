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
#include <Walls/InfiniteWall.h>
#include "Mercury3D.h"
#include "Particles/SuperQuadricParticle.h"
#include "Species/LinearViscoelasticSpecies.h"
#include "Math/ExtendedMath.h"

class BouncingSuperQuadric : public Mercury3D
{
public:
    BouncingSuperQuadric()
    {
        LinearViscoelasticSpecies species;
        species.setCollisionTimeAndRestitutionCoefficient(0.005, 0.8, 1);
        species.setDensity(constants::pi / 6);
        speciesHandler.copyAndAddObject(species);
    
        setTimeStep(1e-4);
    }
    
    void setupInitialConditions() override
    {
        
        SuperQuadricParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setAxesAndExponents(1.0,1.0,1.0,1.0,1.0);
        p0.setInertia();
        
        p0.setPosition(Vec3D(0.0, 0.0, 2.0));
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        p0.setOrientationViaNormal(Vec3D(0.0,1.0,1.0));
        particleHandler.copyAndAddObject(p0);
        logger.assert_always(mathsFunc::isEqual(p0.getMaxInteractionRadius(), 1.0, 1e-10),
                             "interaction radius p0 equals % but should be 1.0", p0.getMaxInteractionRadius());
        
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(1);
        p.setPosition(Vec3D(4.0, 0.0, 2.0));
        p.setVelocity({0.0, 0.0, 0.0});
        particleHandler.copyAndAddObject(p);
        
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,0.0));
        wallHandler.copyAndAddObject(w0);
        
        setTimeMax(3);
        setMin(-1, -1, -1);
        setMax(4, 1, 1);
        setGravity({0,0,-1});
    }
    
    //check contact time
    void actionsAfterTimeStep() override
    {
        logger.assert_always(particleHandler.getObject(0)->getPosition().Z > 0.5 - 1e-2, "Particle fell through wall");
        logger.assert_always(std::abs(particleHandler.getObject(0)->getPosition().Z -
                                              particleHandler.getObject(1)->getPosition().Z) < 1e-0,
                             "BaseParticle and SuperQuad should have same Z-position at time %", getTime());
        logger.assert_always(std::abs(particleHandler.getObject(0)->getVelocity().Z -
                                      particleHandler.getObject(1)->getVelocity().Z) < 1e-0,
                             "BaseParticle and SuperQuad should have the same Z-velocity at time %, but got % and %",
                             getTime(), particleHandler.getObject(1)->getVelocity().Z,
                             particleHandler.getObject(0)->getVelocity().Z);
    }
    
    void test()
    {
        setName("BouncingSuperQuadricSelfTest");
        setSuperquadricParticlesWriteVTK(true);
        setSaveCount(500);
        solve();
    }

};

int main(int argc, char* argv[])
{
    BouncingSuperQuadric problem;
    // comment next line to turn on file output
    problem.setFileType(FileType::NO_FILE);
    problem.setMax(2,2,2);
    problem.setMin(0,0,1);
    problem.setNumberOfDomains({1,1,NUMBER_OF_PROCESSORS});
    problem.test();
    return 0;
}
