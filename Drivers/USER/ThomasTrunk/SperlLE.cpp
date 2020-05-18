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

#include <Boundaries/PeriodicBoundary.h>
#include <MercuryTime.h>
#include <Boundaries/StressStrainControlBoundary.h>
#include "Mercury3D.h"
#include "CG/CG.h"
#include "Boundaries/LeesEdwardsBoundary.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
using constants::pi;
using helpers::to_string;

/**
 * Simulates homogeneous shear for fixed shear-rate and volume fraction.
 *
 * The contact parameters and non-dimensionalisation are the same as for the chute flow simulations:
 *  - Particle diameter: 1 (mono)
 *  - Particle mass: 1
 *  - Stiffness: k=2e5 (equivalent to t_c=1/200)
 *  - Restitution: 0.7
 *  - Friction: 0.5
 *  - Tangential stiffness: 2/7*k
 *  - Tangential dissipation: 2/7*disp
 *  - Pressure: 10-200 (100)
 *  - Inertial number: 1e-3 to 1e-1 (1e-2) => Shear rate: (I/d)*sqrt(p/rho)=0.1
 *  - System size N=1000.
 */
class LE : public Mercury3D
{
public:

    explicit LE(Mdouble pressure, Mdouble shearRate) {
        //set parameters
        Mdouble collisionTime = 2e-2;
        Mdouble restitution = 0.7;
        Mdouble friction = 0.5;
        Mdouble N = 1000;

        //define species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        species->setDensity(6.0 / pi);
        species->setSlidingFrictionCoefficient(0.5);
        species->setCollisionTimeAndNormalAndTangentialRestitutionCoefficient(collisionTime,restitution,restitution,1);

        setTimeStep(0.1*collisionTime);
        setTimeMax(10000);
        setSaveCount(500);
        setName("SperlLE_P"+to_string(pressure,4)+"_S"+to_string(shearRate,4));
        logger(INFO,"Name %",getName());
        setXBallsAdditionalArguments("-v0 -solidf -3dturn 1");

        Mdouble initialVolumeFraction = 0.6;
        Mdouble L = std::cbrt(pi / 6. * static_cast<double>(N) / initialVolumeFraction);
        setDomain({0,0,0},{L,L,L});

        Matrix3D stress, strainRate, gain;
        stress.YY = pressure; gain.YY = 1e-4;
        strainRate.XY = shearRate;
        StressStrainControlBoundary boundary;
        boundary.setHandler(&boundaryHandler);
        boundary.set(stress, strainRate, gain, true);
        boundaryHandler.copyAndAddObject(boundary);

        //define common particle properties
        SphericalParticle p(species);
        p.setRadius(0.5);

        //insert particles
        Vec3D position;
        while (particleHandler.getNumberOfObjects() < N) {
            position.X = random.getRandomNumber(getXMin(), getXMax());
            position.Y = random.getRandomNumber(getYMin(), getYMax());
            position.Z = random.getRandomNumber(getYMin(), getYMax());
            p.setPosition(position);
            particleHandler.copyAndAddObject(p);
        }

        logger(INFO, "%: % particles in 0<x,y,z<%", getName(), N, L);
    }

    void actionsAfterTimeStep() {
        static bool stressControl = true;
        if (stressControl && getTime()>100) {
            logger(INFO,"Turning off stress control at t=100");
            stressControl = false;
        }
    }

    void printTime() const override
    {
        static Time2Finish time(getTime(),getTimeMax());
        Matrix3D stress = getTotalStress();
        auto boundary = dynamic_cast<const StressStrainControlBoundary*>(boundaryHandler.getLastObject());
        logger.assert(boundary,"No StressStrainControlBoundary");
        Mdouble strainRate = boundary->getStrainRate().YY;
        Mdouble shearRate = boundary->getStrainRate().XY;
        Mdouble solidFraction = particleHandler.getVolume()/getTotalVolume();
        Mdouble temperature = getKineticEnergy();
        Mdouble inertialNumber = shearRate/sqrt(stress.YY/solidFraction);
        logger(INFO, "time % confiningStress % strainRate % shearrate % solidFraction % temperature % inertialNumber %", getTime(), stress.YY, strainRate, shearRate, solidFraction, temperature, inertialNumber);
    }

};

int main(int argc, char* argv[])
{
    Mdouble pressure = helpers::readFromCommandLine(argc,argv,"-p",100.0);
    Mdouble shearRate = helpers::readFromCommandLine(argc,argv,"-s",0.1);
    LE leesEdwards(pressure, shearRate);
    leesEdwards.solve();
    return 0;
}
