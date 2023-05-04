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
#include <Species/LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies.h>
#include "Particles/LiquidFilmParticle.h"
#include "DPMBase.h"

class Splash_LS_bed_liquidbridge : public Mercury3D {

public:
    double radius_m, radius_dev;
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
        Mdouble mass_d50 = 4.0/3.0 * constants::pi * pow(radius_m, 3.0) * species->getDensity();
        species->setLiquidBridgeVolumeMax(0.05 * mass_d50/1000.0);
        species->setLiquidBridgeVolumeMin(0.0);
        species->setDistributionCoefficient(0.0);
        species->setSurfaceTension(0.07);
        species->setContactAngle(constants::pi/18);
        species->setViscosity(0.00085);

        const Mdouble collisionTime = species->getCollisionTime(mass);
        setSaveCount(500);
        setTimeMax(0.2);
        setTimeStep(collisionTime/50);

        const unsigned int N1 = static_cast<unsigned int>(pow(N, 0.33)) + 1;
        //particleHandler.clear();
        LiquidFilmParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));


        //PSD with normal distribution
        std::default_random_engine gen;
        //int seed = 3; 也可以设置一个种子
        //std::default_random_engine gen(seed);
        std::normal_distribution<double> dis(radius_m, radius_dev);//radius with gaussian distribution

        for (unsigned int i = 0; i < N; ++i) {
            const unsigned int ix = (i % N1);
            const unsigned int iz = static_cast<unsigned int>( static_cast<double>(i) / N1 / N1);
            const unsigned int iy = (i - ix - N1 * N1 * iz) / N1;
            //std::cout<< ix << iy << iz << std::endl;
            // set particle position
            const double x = (getXMax() - getXMin()) * (ix + 1) / (N1 + 1);
            const double y = (getYMax() - getYMin()) * (iy + 1) / (N1 + 1);
            const double z = (getZMax() - getZMin()) * (iz + 1) / (N1 + 1);
            p0.setPosition(Vec3D(x, y, z));
            //std::cout<< x << y << z << std::endl;
            // set random velocities for the particle
            p0.setVelocity(Vec3D(random.getRandomNumber(-0.0005,0.0005),  random.getRandomNumber(-0.0005,0.0005), 0));
            //p0.setVelocity(Vec3D(0,0,0));
            //p0.setRadius(random.getRandomNumber(0.0003, 0.0007));
            p0.setRadius(dis(gen));//generate particle with radius in normal distribution
            Mdouble mass_p0 = 4.0/3.0 * constants::pi * pow(p0.getRadius(), 3.0) * species->getDensity();
            p0.setLiquidVolume(mass_p0 * 0.05);
            particleHandler.copyAndAddObject(p0);
        }


        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
//        w0.set(Vec3D(1.0, 0.0, 0.0), Vec3D(getXMax(), 0.0, 0.0));
//        wallHandler.copyAndAddObject(w0);
//        w0.set(Vec3D(-1.0, 0.0, 0.0), Vec3D(getXMin(), 0.0, 0.0));
//        wallHandler.copyAndAddObject(w0);
        PeriodicBoundary pb;
        pb.set(Vec3D(1.0, 0.0, 0.0), getXMin(), getXMax());
        boundaryHandler.copyAndAddObject(pb);

        w0.set(Vec3D(0.0, 1.0, 0.0), Vec3D(0.0, getYMax(), 0.0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0, getYMin(), 0.0));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, getZMin()));
        wallHandler.copyAndAddObject(w0);
        w0.set(Vec3D(0.0, 0.0, 1.0), Vec3D(0.0, 0.0, getZMax()));
        wallHandler.copyAndAddObject(w0);
    };


    void printTime() const override {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime() << std::endl;
    }

};

int main(int argc UNUSED, char *argv[] UNUSED){

    Splash_LS_bed_liquidbridge sp0;
    sp0.setName("Splash_LS_bed_liquidbridge");
    sp0.setMax({0.007,0.003,0.005});
    sp0.setGravity({0.0,0.0,-9.8});
    sp0.radius_m = 0.00025;
    sp0.radius_dev = 0.00005;
    sp0.N = 300;//particle number

    sp0.setFileType(FileType::ONE_FILE);
    sp0.setParticlesWriteVTK(true);
    sp0.wallHandler.setWriteVTK(FileType::ONE_FILE);
    sp0.solve();

}
