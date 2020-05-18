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
#include "Mercury3D.h"
#include "CG/CG.h"
#include "Boundaries/LeesEdwardsBoundary.h"
#include "Species/LinearViscoelasticSlidingFrictionSpecies.h"
using constants::pi;

/**
 * Simulates homogeneous shear for fixed shear-rate and volume fraction.
 *
 * The contact parameters and non-dimensionalisation are the same as for the chute flow simulations:
 *  - Particle diameter: 1 (mono)
 *  - Particle mass: 1
 *  - Stiffness: 2e5
 *  - Dissipation: 50
 *  - Friction: 0.5
 *  - Tangential stiffness: 2/7*2e5
 *  - Tangential dissipation: 50
 *
 *  System size is N=1000.
 *
 *  Note, the shear rate and volume fraction relates to certain inclinations and pressures of the corresponding chute flows:
 *  - Inclinations: \theta \in [21,28] deg
 *  - Pressures: P \in [0,40]
 *  - Volume fraction: \phi = 0.610 - exp((\theta-46.2)/7.02)
 *  - Shear rate: \dot\gamma = sqrt(p) * sqrt(pi/6) * 0.617 * (mu-19.67)/(39.89-mu), mu=\tan\theta
 *
 *  So we need to study the following volume fractions and shear rates:
 *  - Volume fraction: 0.535 - 0.582
 *  - Shear rates: 0.165 - 1.618
 *
 * @param argc
 * @param argv
 * @return
 */
class LeesEdwardsInit : public Mercury3D
{
public:

    explicit LeesEdwardsInit(Mdouble volumeFraction) {
        //define species
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        species->setDensity(6.0 / pi);
        species->setStiffness(2e5);
        species->setDissipation(50);
        species->setSlidingStiffness(2. / 7. * 2e5);
        species->setSlidingDissipation(50);
        species->setSlidingFrictionCoefficient(0.5);

        //set gravity, time step
        setGravity(Vec3D(0, 0, 0)); // such that gravity is one
        setTimeStep(1e-4);
        setTimeMax(1e20);
        setSaveCount(std::numeric_limits<unsigned>::max());
        eneFile.setSaveCount(200);
        dataFile.setSaveCount(200);
        fStatFile.setFileType(FileType::NO_FILE);
        setName("LE_V"+helpers::to_string(volumeFraction,4));
        setXBallsAdditionalArguments("-v0 -solidf -3dturn 1");
        //logger(INFO, "File name %", getName());

        //set domain size
        const unsigned N = mathsFunc::cubic(10);
        const Mdouble L = pow(pi / 6. * static_cast<double>(N) / volumeFraction, 1. / 3.);
        setXMin(0);
        setYMin(0);
        setZMin(0);
        setXMax(L);
        setYMax(L);
        setZMax(L);

        //define leesEdwardsBoundary
        LeesEdwardsBoundary leesEdwardsBoundary;
        leesEdwardsBoundary.set(
                [](double time UNUSED) { return 0; },
                [](double time UNUSED) { return 0; },
                getXMin(), getXMax(), getYMin(), getYMax());
        boundaryHandler.copyAndAddObject(leesEdwardsBoundary);

        PeriodicBoundary periodicBoundary;
        periodicBoundary.set({0, 0, 1}, getZMin(), getZMax());
        boundaryHandler.copyAndAddObject(periodicBoundary);

        //define common particle properties
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(0.5);

        Vec3D position;
        while (particleHandler.getNumberOfObjects() < N) {
            position.X = random.getRandomNumber(getXMin(), getXMax());
            position.Y = random.getRandomNumber(getYMin(), getYMax());
            position.Z = random.getRandomNumber(getYMin(), getYMax());
            p.setPosition(position);
            particleHandler.copyAndAddObject(p);
        }

        logger(INFO, "%: % particles in 0<x,y,z<%", getName(), N, L);
        solve();
    }

    void computeExternalForces(BaseParticle* particle) override
    {
        particle->addForce(-50*particle->getVelocity());
    }

    void actionsBeforeTimeStep() override
    {
        if (getNumberOfTimeSteps()%100!=0) return;

        const Mdouble ene =  getKineticEnergy();
        if (ene<1.0) {
            setName(getName()+"_ini");
            setTimeMax(getTime());
        }
        //logger(INFO, "T % E %", getTime(), getKineticEnergy());
    }

    void printTime() const override
    {
        static Time2Finish time(getTime(),getTimeMax());
        logger(INFO, "T % Ek % E % t %", getTime(), getKineticEnergy(), getKineticEnergy() + getElasticEnergy(),time.getTime2Finish(getTime()));
    }

};

class LeesEdwards : public Mercury3D
{
public:

    LeesEdwards(Mdouble volumeFraction, Mdouble shearRate) {
        setName("LE_V"+helpers::to_string(volumeFraction,4)+"_ini");
        while (!readRestartFile()) {
            logger(INFO, "Creating initial conditions %", getName());
            LeesEdwardsInit leesEdwardsInit(volumeFraction);
        }
        //logger(INFO, "Read initial conditions %", getName());
        setName("LE_V"+helpers::to_string(volumeFraction,4)+"_S"+helpers::to_string(shearRate,4));

        //define leesEdwardsBoundary
        const Mdouble L = getXMax()-getXMin();
        const Mdouble V = L * shearRate;
        auto* leesEdwardsBoundary = dynamic_cast<LeesEdwardsBoundary*>(boundaryHandler.getObject(0));
        leesEdwardsBoundary->set(
                [V](double time) { return time * V; },
                [V](double time UNUSED) { return V; },
                getXMin(), getXMax(), getYMin(), getYMax());

        setTimeMax(50);
        dataFile.setSaveCount(1000);
        logger(INFO, "%: Volume fraction %, shear rate %", getName(), volumeFraction, shearRate);

        cg = cgHandler.copyAndAddObject(CG<CGCoordinates::O>());
        cg->statFile.setSaveCount(50);
        solve();
    }

    void printTime() const override
    {
        static Time2Finish time(getTime(),getTimeMax());
        Matrix3D stress = cg->getPoint(0).getContactStress();
        Mdouble P = stress.trace()/3.0;
        Mdouble T = stress.deviator();
        logger(INFO, "T % t % P % T % mu %", getTime(), time.getTime2Finish(getTime()),P,T,atan2(T,P)*180./pi);
    }

    CG<CGCoordinates::O>* cg;
};

int main(int argc, char* argv[])
{
    /*
     *  Note, the shear rate and volume fraction relates to certain inclinations and pressures of the corresponding chute flows:
     *  - Inclinations: \theta \in [21,28] deg
     *  - Pressures: P \in [0,40]
     *  - Volume fraction: \phi = 0.610 - exp((\theta-46.2)/7.02)
     *  - Shear rate: \dot\gamma = sqrt(p) * sqrt(pi/6) * 0.617 * (mu-19.67)/(39.89-mu), mu=\tan\theta
     *
     *  So we need to study the following volume fractions and shear rates:
     *  - Volume fraction: 0.535 - 0.582
     *  - Shear rates: 0.165 - 1.618
     */
    if (argc<2) logger(ERROR,"argument needed:\n LE $inclination");
    const Mdouble inclination = atof(argv[1]);
    const Mdouble volumeFraction = 0.610 - exp((inclination-46.2)/7.02);
    const Mdouble mu = tan(inclination*pi/180.);
    const Mdouble mu1 = tan(19.67*pi/180.);
    const Mdouble mu2 = tan(39.89*pi/180.);

    std::array<Mdouble,4> pressures = {40, 20, 10, 5};
    for (const Mdouble pressure : pressures) {
        const Mdouble shearRate = std::sqrt(pressure * pi / 6) * 0.617 * (mu - mu1) / (mu2 - mu);
        logger(INFO, "Inclination %, pressure %", inclination, pressure);
        logger(INFO, "Volume fraction %, shear rate %", volumeFraction, shearRate);
        LeesEdwards leesEdwards(volumeFraction, shearRate);
        leesEdwards.solve();
    }
}
