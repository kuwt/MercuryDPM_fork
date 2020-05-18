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

#include <Walls/InfiniteWall.h>
#include <CG/TimeAveragedCG.h>
#include <random>
#include "Mercury3D.h"
#include "CG/CG.h"
#include "Species/LinearViscoelasticSpecies.h"

class Packing : public Mercury3D
{
public:

    /**
     * Defines the problem setup
     */
    void setupInitialConditions() override {
        //define a particle species
        LinearViscoelasticSpecies species;
        species.setDensity(6.0 / constants::pi);
        species.setCollisionTimeAndRestitutionCoefficient(50e-4, 0.1, 1.0);
        speciesHandler.copyAndAddObject(LinearViscoelasticSpecies());

        // set number and radius of particles
        unsigned nParticles = 10000;
        Mdouble meanRadius = 1;
        Mdouble stdRadius = .1;
        //static std::random_device rd;
        //static std::mt19937 gen(rd());
        static std::mt19937 gen(0);
        static std::normal_distribution<> d(meanRadius, stdRadius);

        setDomain({0,0,0},{1,1,1});

        particleHandler.setStorageCapacity(nParticles);
        SphericalParticle p(speciesHandler.getObject(0));
        for (unsigned i = 0; i < nParticles; ++i)
        {
            p.setRadius(std::fmin(2*meanRadius,std::fmax(0,d(gen))));
            //p.setRadius(d(gen));
            p.setPosition(p.getPosition()+p.getRadius()*Vec3D(0,0,2));
            particleHandler.copyAndAddObject(p);
        }
    }
};

int main()
{
    logger(INFO,"Testing if the particle size momenta are correctly computed by the cg implementation.");
    Packing dpm;
    dpm.setName("CGMomentSelfTest");
    dpm.setTimeStep(1e-4);
    dpm.setTimeMax(0);
    auto cg = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::O>());
    dpm.solve();
    //cg->getPoint(0).outputStandardisedParticleSizeMomenta(std::cout);
    auto momenta  = cg->getPoint(0).getStandardisedParticleSizeMomenta();

    //Checks if moments are computed correctly
    helpers::check(momenta[0],10000,1e-5, "Particle Number");
    helpers::check(momenta[1],1.001,1e-5, "Mean Radius");
    helpers::check(momenta[2],0.0101141,1e-5, "Standard deviation");
    helpers::check(momenta[3],-0.0198724,1e-5, "Skewness");
    helpers::check(momenta[4],2.91938,1e-5, "Kurtosis");
    helpers::check(momenta[5],-0.0762574,1e-5, "5-th moment");
    return 0;
}

