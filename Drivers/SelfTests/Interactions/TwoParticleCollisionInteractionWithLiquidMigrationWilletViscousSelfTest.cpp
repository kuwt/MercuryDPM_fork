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
#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include <Species/LinearViscoelasticFrictionLiquidMigrationWilletViscousSpecies.h>
#include "Particles/LiquidFilmParticle.h"
#include "DPMBase.h"

/// In this file two particles are symmetrically placed in a domain with opposite initial velocities, to test the liquid bridge force after colliding.
class TwoParticleCollisionInteraction: public Mercury3D {

    public:
        double Radius, x,y,z, vel_x, vel_y, vel_z;


        void setupInitialConditions () override {

            auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionLiquidMigrationWilletViscousSpecies());
            species->setDensity(2500);
            Mdouble mass = species->getMassFromRadius(Radius);
            species->setStiffness(1500);
            species->setDissipation(0.002);
            species->setSlidingStiffness(0.4*species->getStiffness());
            species->setSlidingDissipation(0.4*species->getDissipation());
            species->setSlidingFrictionCoefficient(0.3);
            species->setRollingStiffness(0.4*species->getStiffness());//1.2e4);
            species->setRollingDissipation(0.4*species->getDissipation());//6.3e-2);
            species->setRollingFrictionCoefficient(0.5);
            species->setLiquidBridgeVolumeMax(1.06E-09);
            species->setLiquidBridgeVolumeMin(0.0);
            species->setDistributionCoefficient(1);
            species->setSurfaceTension(0.0728);
            species->setContactAngle(35*constants::pi/180);
            species->setViscosity(0.001);


            const Mdouble collisionTime = species->getCollisionTime(mass);
            setTimeMax(getXMax()/4/fabs(vel_x));
            setTimeStep(collisionTime/50);
            logger(INFO,"timestep=%",getTimeStep());
            setSaveCount(100);

            LiquidFilmParticle p0;
            p0.setRadius(Radius);
            p0.setSpecies(speciesHandler.getObject(0));
            p0.setPosition(Vec3D(getXMax()/2,getYMax()/2,getZMax()/2));
            p0.setVelocity(Vec3D(-vel_x, vel_y, vel_z));
            p0.setLiquidVolume(1.06E-09);
            particleHandler.copyAndAddObject(p0);
            p0.setRadius(Radius);
            p0.setSpecies(speciesHandler.getObject(0));
            p0.setPosition(Vec3D(getXMax()/4,getYMax()/2,getZMax()/2));
            p0.setVelocity(Vec3D(vel_x, vel_y, vel_z));
            p0.setLiquidVolume(1.06E-09);
            particleHandler.copyAndAddObject(p0);

            InfiniteWall w0;
            w0.setSpecies(speciesHandler.getObject(0));
            w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, getZMin()));
            wallHandler.copyAndAddObject(w0);
        };


        void printTime() const override {
            std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime() << std::endl;
        }

    };

    int main(int argc UNUSED, char *argv[] UNUSED){
        TwoParticleCollisionInteraction wpw;

        wpw.setName("TwoParticleCollisionInteractionWithLiquidMigrationWilletViscousSelfTest");
        wpw.setMax({0.01,0.01,0.02});
        wpw.setGravity({0.0,0.0,0.0});

        wpw.Radius = 1.74e-3*0.5;//d=1.74mm

        wpw.vel_x = 0.4;
        wpw.vel_y = 0;
        wpw.vel_z = 0;

        wpw.setParticlesWriteVTK(true);
        wpw.setWallsWriteVTK(FileType::ONE_FILE);
        wpw.solve(argc,argv);
        return 0;
    }