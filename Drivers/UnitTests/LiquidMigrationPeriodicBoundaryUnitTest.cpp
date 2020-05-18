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
#include "DPMBase.h"
#include "Particles/LiquidFilmParticle.h"
#include "Boundaries/PeriodicBoundary.h"
#include <Species/LinearViscoelasticFrictionLiquidMigrationWilletSpecies.h>
using mathsFunc::isEqual;

class LiquidMigrationPeriodicBoundaryInteraction : public DPMBase {
public:

    Mdouble getLiquidFilmVolume() const {
        Mdouble liquidFilmVolume = 0;
        for (auto& p : particleHandler) {
            auto l = dynamic_cast<LiquidFilmParticle*>(p);
            logger.assert_always(l!= nullptr,"Error in printTime");
            liquidFilmVolume += l->getLiquidVolume();
        }
        return liquidFilmVolume;
    };

    Mdouble getLiquidBridgeVolume() const {
        Mdouble liquidBridgeVolume = 0;
        for (auto& i : interactionHandler) {
            LiquidMigrationWilletInteraction* l = dynamic_cast<LiquidMigrationWilletInteraction*>(i);
            logger.assert_always(l!= nullptr,"Error in printTime");
            liquidBridgeVolume += l->getLiquidBridgeVolume();
            //logger(INFO,"o %",l->getOverlap());
        }
        return liquidBridgeVolume;
    };

    Mdouble getLiquidVolume() const {
        return getLiquidFilmVolume() + getLiquidBridgeVolume();
    };

    void printTime() const override {
        logger(INFO,"t % LF % LB % T %",getNumberOfTimeSteps(),getLiquidFilmVolume(), getLiquidBridgeVolume(), getLiquidVolume());
    };
};

int main()
{
    logger(INFO,"This unit test checks that liquid bridge migration is not affected by the slight inaccuracy created "
     "by having ghost particles. A pair of particles is created that touches right at the periodic boundary.");

    logger(INFO,"Testing formation");

    LiquidMigrationPeriodicBoundaryInteraction dpm;

    //set dpm properties
    {
        dpm.setName("LiquidMigrationPeriodicBoundaryInteraction");
        dpm.setTimeStep(1e-4);
        dpm.setTimeMax(5e-4);
        dpm.setMin({0,-.5,-.5});
        dpm.setMax({5,.5,.5});
        dpm.setSaveCount(1);
        dpm.setGravity({0,0,0});
    }

    //add species
    {
        LinearViscoelasticFrictionLiquidMigrationWilletSpecies s;
        s.setDensity(6. / constants::pi);
        s.setCollisionTimeAndRestitutionCoefficient(50*dpm.getTimeStep(), .9, 1);
        s.setLiquidBridgeVolumeMax(1.1);
        dpm.speciesHandler.copyAndAddObject(s);
    }

    //add boundary
    {
        PeriodicBoundary b;
        b.set({1,0,0},-.5,4.5);
        dpm.boundaryHandler.copyAndAddObject(b);
    }

    //add particles
    {
        LiquidFilmParticle p;
        p.setRadius(.5);
        p.setSpecies(dpm.speciesHandler.getLastObject());
        p.setLiquidVolume(.5);

        p.setPosition({4,0,0});
        dpm.particleHandler.copyAndAddObject(p);

        //p.setPosition({0,0,0});
        p.setPosition({2*pow(2,-52),0,0}); //touches on left boundary, but not on right
        dpm.particleHandler.copyAndAddObject(p);
    }

    //solve should cause a bridge to form
    Mdouble initialLiquidVolume = dpm.getLiquidVolume();
    dpm.printTime();
    dpm.solve();

    logger.assert_always(isEqual(dpm.getLiquidVolume(), initialLiquidVolume, 1e-14),
                         "Liquid Volume not conserved (%!=%)", dpm.getLiquidVolume(), initialLiquidVolume);

    logger(INFO,"Testing rupture");

    //move particles, causing a rupture
    {
        BaseParticle* p = dpm.particleHandler.getObject(0);
        p->setPosition({3,0,0});
    }


    //solve should cause a bridge to rupture
    dpm.setRestarted(true);
    dpm.setTimeMax(2.0*dpm.getTimeMax());
    dpm.printTime();
    dpm.solve();

    logger.assert_always(isEqual(dpm.getLiquidVolume(), initialLiquidVolume, 1e-14),
                         "Liquid Volume not conserved (%!=%)", dpm.getLiquidVolume(), initialLiquidVolume);



    return 0;
}
