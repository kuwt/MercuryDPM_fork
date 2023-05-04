//
// Created by WangX3 on 22-2-2022.
//
//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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
#include <random>
#include "Mercury2D.h"
#include "Walls/InfiniteWall.h"
#include <Boundaries/PeriodicBoundary.h>
#include "Species/LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Particles/LiquidFilmParticle.h"
#include "DPMBase.h"
#include <iostream>
#include <fstream>

class Splash_LS : public Mercury2D {

public:
    double radius, vel_x, vel_y, vel_z;

    void setupInitialConditions () override{

            //setName();//set the name of restart file to be read
            readRestartFile("Splash2D_LS_bed_000");
            setRestarted(false);
            setName("Splash2D_changeLV005");

            interactionHandler.clear();//clear all the interactions and detect again
            auto originalSpecies = dynamic_cast<LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies*>(speciesHandler.getObject(0));
            const Mdouble mass = originalSpecies->getMassFromRadius(radius);
            originalSpecies->setLiquidBridgeVolumeMax(mass * 0.05 * 2/1000.0);

            //change the liquid volume for each particle
            for (BaseParticle* bp : particleHandler) {
            auto lfp = dynamic_cast<LiquidFilmParticle*>(bp);//cast the baseparticle to LF type
            logger.assert_debug(lfp != nullptr, "Your particles need to be of type LiquidFilmParticle for this code to work.");
            Mdouble mass_lfp = 4.0/3.0 * constants::pi * pow(lfp->getRadius(), 3.0) * originalSpecies->getDensity();
            lfp->setLiquidVolume(0.05 * mass_lfp/1000.0);
            }

            //add the impact particle
            LiquidFilmParticle particle;
            particle.setSpecies(speciesHandler.getObject(0));
            particle.setRadius(radius);
            particle.setPosition(Vec3D(2*radius, getYMax()*2.0/3.0, 0.0));
            particle.setVelocity(Vec3D(vel_x, vel_y, vel_z));
            Mdouble mass_impact = 4.0/3.0 * constants::pi * pow(particle.getRadius(), 3.0) * originalSpecies->getDensity();
            particle.setLiquidVolume(0.05 * mass_impact/1000.0);
            particleHandler.copyAndAddObject(particle);



        //const Mdouble collisionTime = species->getCollisionTime(mass);
        setSaveCount(100);
        setTimeMax(0.2);
        setTimeStep(1e-5);
        }


    void printTime() const override {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime() << std::endl;
    }

//public:
//    LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies* species;

};

int main(int argc UNUSED, char *argv[] UNUSED){
//    std::string fileName ("Splash_LS_bed");
//    if (argc>1)
//    {
//        fileName = argv[1];
//        std::cout << "restarting from " << fileName << std::endl;
//    }


    Splash_LS sp0;

    sp0.setMax({0.0085,0.008,0.002});
    sp0.setGravity({0.0,-9.8,0.0});

    sp0.radius = 0.00025;

    sp0.vel_x = 0.3;
    sp0.vel_y = -0.2;
    sp0.vel_z = 0;

    sp0.setFileType(FileType::ONE_FILE);
    sp0.setParticlesWriteVTK(true);
    sp0.wallHandler.setWriteVTK(FileType::ONE_FILE);
    sp0.solve();
}

