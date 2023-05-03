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
#include <random>
#include "Mercury2D.h"
#include "Walls/InfiniteWall.h"
#include <Boundaries/PeriodicBoundary.h>
#include <Species/LinearViscoelasticSlidingFrictionLiquidMigrationLSSpecies.h>
#include "Particles/LiquidFilmParticle.h"
#include "DPMBase.h"

class Drop2D : public Mercury2D {

public:
    double radius;

    void setupInitialConditions () override {


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
        species->setLiquidBridgeVolumeMax(0.05*2*mass/1000.0);
        species->setLiquidBridgeVolumeMin(0.0);
        species->setDistributionCoefficient(0.0);
        species->setSurfaceTension(0.07);
        species->setContactAngle(constants::pi/18);
        species->setViscosity(0.00085);

        //const Mdouble collisionTime = species->getCollisionTime(mass);
        setSaveCount(100);
        setTimeMax(0.1);
        setTimeStep(1e-5);


        LiquidFilmParticle p0,p1;
        p0.setSpecies(speciesHandler.getObject(0));
        p1.setSpecies(speciesHandler.getObject(0));
        p0.setPosition(Vec3D(getXMax()/2.0,getYMax()/2.0, 0.0));
        p1.setPosition(Vec3D(getXMax()/2.0,getYMax()/2.0-2.2*radius, 0.0));
        p0.setVelocity(Vec3D(0.0, -0.05, 0.0));
        p1.setVelocity(Vec3D(0.0, -0.05, 0.0));
        p0.setRadius(radius);
        p1.setRadius(radius);
        //Mdouble mass_p = 4.0/3.0 * constants::pi * pow(radius, 3.0) * species->getDensity();
        p0.setLiquidVolume(mass * 0.05/1000.0);
        p1.setLiquidVolume(mass * 0.05/1000.0);
        particleHandler.copyAndAddObject(p0);
        particleHandler.copyAndAddObject(p1);


        InfiniteWall w;
        w.setSpecies(speciesHandler.getObject(0));
        w.set(Vec3D(1,0,0), Vec3D(getXMax(), 0.0, 0.0));
        wallHandler.copyAndAddObject(w);
        w.set(Vec3D(-1,0,0), Vec3D(getXMin(), 0.0, 0.0));
        wallHandler.copyAndAddObject(w);
        w.set(Vec3D(0,1,0), Vec3D(0.0, getYMax(), 0.0));
        wallHandler.copyAndAddObject(w);
        w.set(Vec3D(0,-1,0), Vec3D(0.0, getYMin(), 0.0));
        wallHandler.copyAndAddObject(w);
    };


    void printTime() const override {
        std::cout << "t=" << std::setprecision(3) << std::left << std::setw(6) << getTime() << std::endl;
    }

};

int main(int argc UNUSED, char *argv[] UNUSED){

    Drop2D sp0;
    sp0.setName("DropTest2D_005");
    sp0.setMax({0.005,0.005,0.005});
    sp0.setGravity({0.0,-9.8,0.0});
    sp0.radius = 0.00025;

    sp0.setFileType(FileType::ONE_FILE);
    sp0.setParticlesWriteVTK(true);
    sp0.wallHandler.setWriteVTK(FileType::ONE_FILE);
    sp0.solve();

}
