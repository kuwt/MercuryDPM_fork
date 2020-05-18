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
#include "Particles/LiquidFilmParticle.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include "Species/Species.h"
#include "Species/LinearViscoelasticFrictionLiquidMigrationWilletSpecies.h"
#include <iomanip>

/// In this file two particles of different species are symmetrically placed in a bi-axial box are allowed to jump around under gravity. It tests walls gravity and symmetry.
class DPM : public Mercury2D
{

    void setupInitialConditions() override {
        setName("LiquidMigrationTwoSpeciesSelfTest");
        setSystemDimensions(3);
        setXBallsAdditionalArguments("-v0 -solid -3dturn 1");
        
        //define 1st particle species
        auto species0 = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionLiquidMigrationWilletSpecies());
        species0->setStiffness(2e5);
        species0->setDensity(6.0/constants::pi);
        //define 2nd particle species/ wall species
        auto species1 = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionLiquidMigrationWilletSpecies());
        species1->setStiffness(2e5);
        species1->setDensity(6.0/constants::pi);
        species1->setStiffness(2e5);
        //define mixed species (the one that will be used
        //define the mixed species such that it has a finite interaction distance, but the particle species' don't
        auto species01 = speciesHandler.getMixedObject(species0,species1);
        species01->setStiffness(2e5);
        species01->setLiquidBridgeVolumeMax(mathsFunc::cubic(1e-3));
        species01->setDistributionCoefficient(1.0);
        species01->setSurfaceTension(1.0);
        species01->setContactAngle(0.0);
        
        logger(INFO,"species0 ID % maxID %",species0->getInteractionDistance(),species0->getMaxInteractionDistance());
        logger(INFO,"species1 ID % maxID %",species1->getInteractionDistance(),species1->getMaxInteractionDistance());
        logger(INFO,"species01 ID %",species01->getInteractionDistance());
        
        
        // simulate a particle-particle and a particle-wall collision
        setTimeMax(0.115);
        //setTimeMax(0.0561);
        setSaveCount(1);
        setTimeStep(1e-4);
        setGravity(Vec3D(0, 0, 0));
        setMax({.5,.5,.505});
        setMin(-getMax());
        
        LiquidFilmParticle p;
        p.setSpecies(species0);
        p.setLiquidVolume(mathsFunc::cubic(1e-3));
        p.setRadius(0.5);
        auto q = particleHandler.copyAndAddObject(p);
        logger(INFO,"particle IR %",q->getMaxInteractionRadius());
        
        p.setSpecies(species1);
        p.setPosition(Vec3D(0,0,getZMin()-p.getRadius()));
        p.setVelocity(Vec3D(0,0,.1));
        particleHandler.copyAndAddObject(p);
        
        wallHandler.clear();
        InfiniteWall w;
        w.setSpecies(species1);
        w.set(Vec3D(0, 0, 1), Vec3D(0, 0, getZMax()));
        wallHandler.copyAndAddObject(w);
    }
};

int main()
{
    DPM dpm;
    dpm.solve();
}
