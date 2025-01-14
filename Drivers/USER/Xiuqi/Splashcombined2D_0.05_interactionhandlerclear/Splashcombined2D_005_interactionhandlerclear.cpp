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
#include <Species/LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies.h>
#include "Particles/LiquidFilmParticle.h"
#include "DPMBase.h"

class Splashcombined2D : public Mercury2D {

public:
    double radius_m, radius_dev, vel_x, vel_y, vel_z;
    unsigned int N;

    void setupInitialConditions () override {


        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies());
        species->setDensity(2650);
        species->setStiffness(1500);
        species->setDissipation(0.002);
        const Mdouble mass = species->getMassFromRadius(radius_m);
        species->setRestitutionCoefficient(0.7, mass);
        species->setSlidingStiffness(0.4*species->getStiffness());
        species->setSlidingDissipation(0.4*species->getDissipation());
        species->setSlidingFrictionCoefficient(0.1);
        species->setSlidingFrictionCoefficientStatic(0.1);
        species->setLiquidBridgeVolumeMax(0.0);
        species->setLiquidBridgeVolumeMin(0.0);
        species->setDistributionCoefficient(0.0);
        species->setSurfaceTension(0.07);
        species->setContactAngle(constants::pi/18);
        species->setViscosity(0.00085);

        //const Mdouble collisionTime = species->getCollisionTime(mass);
        setSaveCount(100);
        setTimeMax(0.2);
        setTimeStep(1e-5);


        std::default_random_engine gen;
        std::normal_distribution<double> dis(radius_m, radius_dev);
        int N1=static_cast<int>(sqrt(N))+1;
        LiquidFilmParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        for (int i=0;i<N;i++) {
            int ix = i % N1;
            int iy = i / N1;
            // set particle position
            double x = (getXMax() - getXMin()) * (ix + 1) / (N1 + 1);
            double y = (getYMax() - getYMin()) * (iy + 1) / (N1 + 1);
            p0.setPosition(Vec3D(x, y, 0.0));
            // set random velocities for the particle
            p0.setVelocity(Vec3D(random.getRandomNumber(-0.0005,0.0005),  random.getRandomNumber(-0.0005,0.0005), 0.0));
            p0.setRadius(dis(gen));
            Mdouble mass_p0 = 4.0/3.0 * constants::pi * pow(p0.getRadius(), 3.0) * species->getDensity();
            p0.setLiquidVolume(mass_p0 * 0.05/1000.0);
            particleHandler.copyAndAddObject(p0);
        }


        PeriodicBoundary pb;
        pb.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(pb);

        InfiniteWall w;
        w.setSpecies(speciesHandler.getObject(0));
        w.set(Vec3D(0,1,0), Vec3D(0.0, getYMax(), 0.0));
        wallHandler.copyAndAddObject(w);
        w.set(Vec3D(0,-1,0), Vec3D(0.0, getYMin(), 0.0));
        wallHandler.copyAndAddObject(w);
    };

    void actionsAfterTimeStep() override {
        if (getTime() < 0.1 && getTime()+getTimeStep()>=0.1) {
            interactionHandler.clear();//clear all the interaction
            auto originalSpecies = dynamic_cast<LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies*>(speciesHandler.getObject(0));
            const Mdouble mass = originalSpecies->getMassFromRadius(radius_m);
            originalSpecies->setLiquidBridgeVolumeMax(mass * 0.05 * 2/1000.0);
        }
        if (getTime() < 0.12 && getTime()+getTimeStep()>=0.12) {
            //add the impact particle
            LiquidFilmParticle particle;
            auto originalSpecies = dynamic_cast<LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies*>(speciesHandler.getObject(0));
            particle.setSpecies(speciesHandler.getObject(0));
            particle.setRadius(radius_m);
            particle.setPosition(Vec3D(2*radius_m, getYMax()*2.0/3.0, 0.0));
            particle.setVelocity(Vec3D(vel_x, vel_y, vel_z));
            Mdouble mass_impact = 4.0/3.0 * constants::pi * pow(particle.getRadius(), 3.0) * originalSpecies->getDensity();
            particle.setLiquidVolume(0.05 * mass_impact/1000.0);
            particleHandler.copyAndAddObject(particle);
        }
    }

    void printTime() const override {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime() << std::endl;
    }


};

int main(int argc UNUSED, char *argv[] UNUSED){

    Splashcombined2D sp0;
    sp0.setName("Splashcombined2D_005_interactionhandlerclear");
    sp0.setMax({0.0085,0.008,0.002});
    sp0.setGravity({0.0,-9.8,0.0});
    sp0.radius_m = 0.00025;
    sp0.radius_dev = 0.00003;//lower than 3D case
    sp0.N = 100;//particle number

    sp0.vel_x = 0.3;
    sp0.vel_y = -0.2;
    sp0.vel_z = 0;

    sp0.setFileType(FileType::ONE_FILE);
    sp0.setParticlesWriteVTK(true);
    sp0.wallHandler.setWriteVTK(FileType::ONE_FILE);
    sp0.setInteractionsWriteVTK(true);
    sp0.solve();

}
