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
#include <Species/LinearViscoelasticFrictionLiquidMigrationWilletSpecies.h>
#include <Particles/LiquidFilmParticle.h>
#include <CG/Fields/LiquidMigrationFields.h>
#include "Mercury3D.h"
#include "CG/CG.h"
using namespace constants;

class TwoParticles : public Mercury3D
{
public:
    void setupInitialConditions() override {
        setName("CGLiquidsUnitTest");
        setTimeStep(1e-4);
        setFileType(FileType::NO_FILE);

        //define a particle species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionLiquidMigrationWilletSpecies());
        species->setDensity(0.75 / constants::pi); // such that mass = 1
        species->setStiffnessAndRestitutionCoefficient(2e5, 0.1, 1.0);
        species->setLiquidBridgeVolumeMax(1);

        //set domain size
        setDomain({0,0,0},{1,1,1});

        //define two particles
        LiquidFilmParticle p;
        p.setSpecies(species);
        p.setRadius(.25);
        p.setLiquidVolume(1);
        p.setPosition({.5, .5, .3});
        particleHandler.copyAndAddObject(p);
        p.setPosition({.5, .5, .7});
        particleHandler.copyAndAddObject(p);
    }
};

int main()
{
    logger(INFO," Simulates a particle-particle collision.\n"
                " Checks cg \n"
                "   - at one instant in time\n"
                "   - standard fields (density and stress)\n"
                "   - for O, Z, and XYZ coordinates\n"
                "   - for both Gauss and Lucy kernel functions\n");

    TwoParticles dpm;
    dpm.setTimeMax(0);
    auto o = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::O,CGFunctions::Gauss,CGFields::LiquidMigrationFields>());
    auto z = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::Z,CGFunctions::Gauss,CGFields::LiquidMigrationFields>());
    z->setWidth(0.16);
    z->setN(100);
    auto xz = dpm.cgHandler.copyAndAddObject(CG<CGCoordinates::XZ,CGFunctions::Gauss,CGFields::LiquidMigrationFields>());
    xz->setWidth(0.16);
    xz->setN(100);
    dpm.solve();

    //A few checks
    logger(INFO,"Checking a few cg parameters");

    //Checks density for unresolved cg, rho = M/V = 2/16 = 0.125
    helpers::check(o->getPoint(0).getLiquidFilmVolume(),1,1e-15, "Average liquidFilmVolume");
    helpers::check(o->getPoint(0).getLiquidBridgeVolume(),1,1e-15, "Average liquidBridgeVolume");
    helpers::check(z->evaluateAverage().getLiquidFilmVolume(),1,3e-2, "Average liquidFilmVolume for z-resolved stats");
    helpers::check(z->evaluateAverage().getLiquidBridgeVolume(),1,1e-2, "Average liquidBridgeVolume for z-resolved stats");
    helpers::check(xz->evaluateAverage().getLiquidFilmVolume(),1,3e-2, "Average liquidFilmVolume for xz-resolved stats");
    helpers::check(xz->evaluateAverage().getLiquidBridgeVolume(),1,1e-2, "Average liquidBridgeVolume for xz-resolved stats");

    return 0;
}

