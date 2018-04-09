//Copyright (c) 2013-2017, The MercuryDPM Developers Team. All rights reserved.
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
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Species/Species.h"
#include "Species/HertzianViscoelasticMindlinSpecies.h"
#include "Species/HertzianViscoelasticSlidingFrictionSpecies.h"

/// In this file two particles are symmetrically placed in a bi-axial box are allowed to jump around under gravity. It tests walls gravity and symmetry.

class TwoParticleCollisionWithMindlinSelfTest : public Mercury2D
{


    void setupInitialConditions()
    {
        setMax({0.01,0.01,0.0});
        setGravity({0.0,-9.8,0.0});

        BaseParticle p0, p1;
        p0.setSpecies(speciesHandler.getObject(0));
        p1.setSpecies(speciesHandler.getObject(0));

        p0.setPosition(Vec3D(0.006, 0.0059, 0.0));
        p1.setPosition(Vec3D(0.004, 0.005, 0.0));

        p0.setVelocity(Vec3D(-0.1, 0.0, 0.0));
        p1.setVelocity(Vec3D(0.1, 0.0, 0.0));

        p0.setRadius(0.0005);
        p1.setRadius(0.0005);
        particleHandler.copyAndAddObject(p0);
        particleHandler.copyAndAddObject(p1);

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

public:
    Mdouble scaleUp = 4; //by how much the particle size is multiplied in test simulations (set to 1 for real simulations)
    Mdouble radius = 1.5e-3*scaleUp;
};

int main(int argc UNUSED, char* argv[] UNUSED)
{

   //The monosized 3 mm perfect-spherical tungsten particles are simulated.
    //Poisson Ratio is 0.28,
    //bulk density is 19250 Kg/m^3,
    //Young modulus 4.11E11 Pa,
    //restitution coefficient is 0.95,
    //Hertz-Mindlin model is used to describe the contact force between neighboring particles.
    Mdouble poissonRatio = 0.28;
    Mdouble density = 19250;
    Mdouble elasticModulus = 4.11e11; //100 times too soft
    Mdouble restitution = 1.0;

    {
        //simulate two particles with a Hertz-Mindlin contact force
        TwoParticleCollisionWithMindlinSelfTest tpcwm;
        tpcwm.setName("TwoParticleCollisionWithMindlinSelfTest");

        HertzianViscoelasticMindlinSpecies species;
        species.setDensity(density);
        species.setElasticModulusAndRestitutionCoefficient(elasticModulus, restitution);
        //https://en.wikipedia.org/wiki/Shear_modulus#References
        species.setShearModulus(0.5 * elasticModulus / (1 + poissonRatio));
        species.setSlidingFrictionCoefficient(0.1);
        tpcwm.speciesHandler.copyAndAddObject(species);

        tpcwm.setTimeMax(0.25);
        tpcwm.setSaveCount(1e8);
        tpcwm.fStatFile.setSaveCount(10);
        Mdouble relativeVelocity = 1;
        Mdouble tc = species.getCollisionTime(2.0 * tpcwm.radius, species.getDensity(), relativeVelocity);
        logger(INFO, "Collision time %", tc);
        tpcwm.setTimeStep(0.02 * tc);
        tpcwm.solve();
    }

    {
        //now do the same with the linear law
        TwoParticleCollisionWithMindlinSelfTest tpc;
        tpc.setName("TwoParticleCollisionWithSlidingFrictionSelfTest");

        HertzianViscoelasticSlidingFrictionSpecies species;
        species.setDensity(density);
        species.setElasticModulusAndRestitutionCoefficient(elasticModulus, restitution);
        //https://en.wikipedia.org/wiki/Shear_modulus#References
        species.setSlidingStiffness(10000000);
        species.setSlidingFrictionCoefficient(0.1);
        tpc.speciesHandler.copyAndAddObject(species);

        tpc.setTimeMax(0.25);
        tpc.setSaveCount(1e8);
        tpc.fStatFile.setSaveCount(10);
        Mdouble relativeVelocity = 1;
        Mdouble tc = species.getCollisionTime(2.0 * tpc.radius, species.getDensity(), relativeVelocity);
        logger(INFO, "Collision time %", tc);
        tpc.setTimeStep(0.02 * tc);
        tpc.solve();
    }

    helpers::writeToFile("TwoParticleCollisionWithMindlinSelfTest.gnu",
                         "p 'TwoParticleCollisionWithMindlinSelfTest.fstat' u 9:10 w p,0.1*x");
    helpers::writeToFile("TwoParticleCollisionWithSlidingFrictionSelfTest.gnu",
                         "p 'TwoParticleCollisionWithSlidingFrictionSelfTest.fstat' u 9:10 w p,0.1*x");
    logger(INFO,"type 'gnuplot', then 'load TwoParticleCollisionWithMindlinSelfTest.gnu' to check force law");
    logger(INFO,"type 'gnuplot', then 'load TwoParticleCollisionWithSlidingFrictionSelfTest.gnu' to check force law");
}
