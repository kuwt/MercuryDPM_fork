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
#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include <Boundaries/PeriodicBoundary.h>
#include "Species/LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
#include "Particles/LiquidFilmParticle.h"
#include "DPMBase.h"
#include <iostream>
#include <fstream>

class Splash_LS : public Mercury3D {

public:
    double radius, vel_x, vel_y, vel_z;

    void setupInitialConditions () override{

            //setName();//set the name of restart file to be read
            readRestartFile("Splash_bed_LS_002");
            setRestarted(false);
            setName("Splash_LS_LB_002");

            speciesHandler.clear();

            auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies());
            species->setDensity(2650);
            species->setStiffness(1500);
            species->setDissipation(0.002);
            const Mdouble mass = species->getMassFromRadius(radius);
            species->setRestitutionCoefficient(0.7, mass);
            species->setSlidingStiffness(0.4*species->getStiffness());
            species->setSlidingDissipation(0.4*species->getDissipation());
            species->setSlidingFrictionCoefficient(0.1);
            species->setSlidingFrictionCoefficientStatic(0.1);
            Mdouble mass_d50 = 4.0/3.0 * constants::pi * pow(radius, 3.0) * species->getDensity();
            species->setLiquidBridgeVolumeMax(0.02 * 2 * mass_d50/1000.0);
            species->setLiquidBridgeVolumeMin(0.0);
            species->setDistributionCoefficient(0.0);
            species->setSurfaceTension(0.07);
            species->setContactAngle(constants::pi/18);
            species->setViscosity(0.00085);


            for (BaseParticle* p : particleHandler)
                p->setSpecies(species);


            //add the impact particle
            LiquidFilmParticle particle;
            particle.setSpecies(speciesHandler.getObject(0));
            particle.setRadius(radius);
            particle.setPosition(Vec3D(2*radius, getYMax()/2, getZMax()/3*2));
            particle.setVelocity(Vec3D(vel_x, vel_y, vel_z));
            Mdouble mass_impact = 4.0/3.0 * constants::pi * pow(particle.getRadius(), 3.0) * species->getDensity();
            particle.setLiquidVolume(0.02 * mass_impact/1000.0);
            particleHandler.copyAndAddObject(particle);



        const Mdouble collisionTime = species->getCollisionTime(mass);
        setSaveCount(500);
        setTimeMax(0.2);
        setTimeStep(collisionTime/50);
        }


    void printTime() const override {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime() << std::endl;
    }

};

int main(int argc UNUSED, char *argv[] UNUSED){
    Splash_LS sp0;

    sp0.setMax({0.007,0.003,0.005});
    sp0.setGravity({0.0,0.0,-9.8});

    sp0.radius = 0.00025;

    sp0.vel_x = 0.7;
    sp0.vel_y = 0;
    sp0.vel_z = -0.4;

    sp0.setFileType(FileType::ONE_FILE);
    sp0.setParticlesWriteVTK(true);
    sp0.wallHandler.setWriteVTK(FileType::ONE_FILE);
    sp0.solve();
}

