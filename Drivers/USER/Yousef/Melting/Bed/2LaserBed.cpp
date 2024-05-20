//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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
#include "Species/MeltableSpecies.h"
#include "Particles/ThermalParticle.h"
#include "Boundaries/PeriodicBoundary.h"
#include <random>



class LaserBed: public Mercury3D
{
private:
    // set main simulation variables
    Mdouble timeStep = 2.7e-6; //2.7e-6 -> MS1e8 //2.7e-7 -> MS1e6 //2.7e-8 -> MS1e4 //2.7e-10 -> MS1
    Mdouble scaleMass = 1e8;
    //
    unsigned int laserON = 1;
    //
    unsigned int wallsBoundary = 1;
    //
    bool writeVTK = true;
    Mdouble MaxSimTime = 10.0;
    unsigned int numberOfCounts = 1e4;
    //
    int thermalExpansionOnOff = 0;
    //
    bool addNewLayer = false;
    //
    // control system size
    int spatialScale = 2;//4
    int volumeScale = 6; //3 //6
    // DPM parameters
    Mdouble gravityValue = -9.81/sqrt(scaleMass); //scaleMass;
    //
    Mdouble XMaxDomain = spatialScale*500e-6/2;
    Mdouble YMaxDomain = spatialScale*500e-6/2;
    Mdouble ZMinDomain = 400e-6;//200e-6;
    Mdouble ZMaxDomain = 100e-6;
    //
    Mdouble MatDensity = 1050 * scaleMass;
    Mdouble elasticModulus = 0.8*2.95e9; //8e7;
    Mdouble possionsRation = 0.4;
    Mdouble damping = 0.9; // * / sqrt(scaleMass);
    //
    Mdouble latentHeat = 56.4e3;
    Mdouble heatCapacity = 1200;
    //
    Mdouble surfaceTension = 34.3e-3 * 2e-3 * sqrt(scaleMass);
    Mdouble refVis = 2e-3 * sqrt(scaleMass);
    // cond , conv, rad
    Mdouble thermalcond = 0.12 * scaleMass;
    Mdouble thermalConvCoff = 150.0 * scaleMass;
    Mdouble emmisivity = 0.9 * scaleMass;
    //
    Mdouble aTemp = 155 + 273.15;//428.15
    Mdouble mTemp = 178 + 273.15;//451.15;
    Mdouble deltaT = 20;//+273.15;// (this value matched initially with Exp, but very high rate) //5.0;
    Mdouble liquidHeatCapacity = 1.2 * heatCapacity;//2.0*heatCapacity;
    Mdouble vaporizationHeatCapacity = 2.5 * heatCapacity;//4.0*heatCapacity;
    Mdouble vaporizationLatentHeat = 20.0 * latentHeat;
    Mdouble vaporizationTemp = 3.5 * mTemp;
    //
    Mdouble thermalExpansionCoeff = 90 * 1e-6; // K-1
    //
    Mdouble particlesTemp = 155 + 273.15;//413.15;
    //
    Mdouble heatingTime = 0.002;
    Mdouble heatSource = 0.0 * scaleMass;
    Mdouble coolingTime = 0.01;
    //
    Mdouble absorptivity = 0.33 * scaleMass;
    Mdouble laserPower = 5;//1;
    Mdouble beamRadius = 250e-6;//125e-6;
    Mdouble powderBedPorosity = 0.55;
    Mdouble initialLaserPX = beamRadius;//-beamRadius;
    Mdouble laserSpeedX = 0.0;//0.1;//-0.1;//100e-6;
    Mdouble initialLaserPY = beamRadius;//-beamRadius;//beamRadius * 2.0 + (beamRadius / 2.0); // second pass -> beamRadius*2.0;// first scan -> beamRadius;
    Mdouble laserSpeedY = 0.0;//0.1;//0.0;//0.1;
    //
    // Domain
    Mdouble ToolLength = spatialScale*500e-6/2;//0.00025; //45e-1 * scale * 2.0; // unit m
    Mdouble ToolWidth = spatialScale*500e-6/2;//0.00025;//30e-1 * scale * 2.0;
    Mdouble ToolHight = 0.01;//40e-1 * scale * 4.0; //20e-1 * scale * 4.0;
    Mdouble ToolThickness = 50e-6;//5e-1 * scale;
    Mdouble Gap = 30e-6;
    //
    //the volume to be added:
    Mdouble addVolume = volumeScale*Gap*ToolLength*ToolWidth;
    //particles:
    //
    Mdouble meanRadius1 = 27e-6;// in m
    Mdouble MinRadius = 20e-6;//22.6e-6/2.0;  //19.5e-6/2.0;//in m
    Mdouble MaxRadius = 75e-6;//153.4/2.0;    //174.5e-6/2.0;//in m
    Mdouble LogmeanRadius = -10.3693;
    Mdouble LogstdRadius = 0.2758;
    //
    MeltableSpecies *species;
    //
public:
    //-------------------------------------------------- --------------------------------------------------
    LaserBed() {
        //
        setName("InitialBed");
        //setName("BedFixAdd");
        //
        //setName("LaserBed");
        //setName("LaserBed2");
        //
        readRestartFile();
        setRestarted(false);
        //
        //removeParticles();
        //
        setName("LaserBed");
        //setName("LaserBed2");
        //setName("LaserBed3");
        species = dynamic_cast<MeltableSpecies *>(speciesHandler.getObject(0));
        //
        setSystemDimensions(3);
        setGravity(Vec3D(0.0,0.0,gravityValue));
        //setXMin(0.0);
        //setYMin(0.0);
        setZMin(-ZMinDomain);
        setXMax(XMaxDomain);
        setYMax(YMaxDomain);
        setZMax(ZMaxDomain);
        setTimeMax(MaxSimTime);
        //
        species->setDensity(MatDensity);
        species->setElasticModulus(elasticModulus);
        species->setPoissonRatio(possionsRation);
        species->setDissipation(damping);
        //
        species->setLatentHeat(latentHeat);
        species->setSolidHeatCapacity(heatCapacity);
        species->setThermalConductivityCoefficient(thermalcond);
        species->setThermalConvectionCoefficient(thermalConvCoff);
        species->setAmbientTemperature(aTemp);
        species->setMeltingTemperature(mTemp);
        species->setMaterialEmissivity(emmisivity);
        //
        species->setSurfaceTension(surfaceTension);
        species->setRefViscosity(refVis);
        //
        species->setDeltaT(deltaT);
        species->setLiquidHeatCapacity(liquidHeatCapacity);
        species->setVaporizationHeatCapacity(vaporizationHeatCapacity);
        species->setVaporizationLatentHeat(vaporizationLatentHeat);
        species->setVaporizationTemperature(vaporizationTemp);
        //
        species->setThermalExpansionCoefficient(thermalExpansionCoeff);
        species->setThermalExpansionOnOff(thermalExpansionOnOff);
        //
        species->setHeatingTime(heatingTime);
        species->setCoolingTime(coolingTime);
        switch (laserON)
        {
            case 1:
            {
                // Laser
                species->setMaterialAbsorptivity(absorptivity);
                species->setLaserPower(laserPower);
                species->setBeamSpotRadius(beamRadius);
                species->setPowderBedPorosity(powderBedPorosity);
                species->setInitialLaserPositionX(initialLaserPX);
                species->setLaserSpeedX(laserSpeedX);
                species->setInitialLaserPositionY(initialLaserPY);
                species->setLaserSpeedY(laserSpeedY);
            }
            case 0:
            {
                species->setHeatInput(heatSource);
            }
        }
        //
        setSaveCount(numberOfCounts);
        setParticlesWriteVTK(writeVTK);
        setWallsWriteVTK(writeVTK);
        //
        //species->setMaxTimeStep(MinRadius,MatDensity,elasticModulus);
        setTimeStep(timeStep);
        //
        if (addNewLayer)
        {
            //
            std::mt19937 gen;
            gen.seed(100);
            std::lognormal_distribution<> d(LogmeanRadius, LogstdRadius);
            ThermalParticle p;
            p.setSpecies(speciesHandler.getObject(0));
            p.setRadius(meanRadius1);
            p.setTemperature(particlesTemp);
            Mdouble fillHeight = 0.0;
            while (addVolume > 0) {
                Mdouble x = random.getRandomNumber(ToolThickness, ToolLength);
                Mdouble y = random.getRandomNumber(ToolThickness, (ToolWidth - ToolThickness));
                Mdouble z = random.getRandomNumber(0.0, fillHeight);
                p.setPosition({x, y, z});
                // check if particle can be inserted
                if (checkParticleForInteraction(p)) {
                    particleHandler.copyAndAddObject(p);
                    addVolume -= p.getVolume();
                    do {
                        p.setRadius(d(gen));
                    } while (p.getRadius() < MinRadius || p.getRadius() > MaxRadius);
                    if (particleHandler.getNumberOfObjects() % 100 == 0) std::cout << '.' << std::flush;
                    if (particleHandler.getNumberOfObjects() % 1000 == 0) std::cout << ' ';
                    if (particleHandler.getNumberOfObjects() % 10000 == 0) std::cout << addVolume << '\n';
                } else {
                    fillHeight += 10e-6 *meanRadius1;
                }
            }
            logger(INFO, " Inserted % particles", particleHandler.getNumberOfObjects());
            //
        }
        //
    };
    //
/*    void removeParticles()
    {
        for (int i = 0; i <= particleHandler.getNumberOfObjects()-1; ++i) {
            Vec3D pP = particleHandler.getObject(i)->getPosition();
            if (getTime()==0.001 & pP.Z >= 0.0)
                particleHandler.removeObject(i);
        }
        logger(INFO, " Removed % particles", particleHandler.getNumberOfObjects());
    }*/
    //void actionsAfterTimeStep() override{}
    //
    //void setupInitialConditions() override {}
    //
    //
    Mdouble getSurfaceHeight() const
    {
        double height = -ZMinDomain; //0;
        for (const BaseParticle* p : particleHandler)
        {
            double newHeight = p->getPosition().Z + p->getRadius();
            if (height<newHeight) height = newHeight;
        }
        return height+ZMinDomain;
    }
    //
    void writeEneTimeStep(std::ostream& os) const override
    {
        Mdouble massTotalParticles = particleHandler.getMass();
        Mdouble volumeTotalParticles = particleHandler.getVolume();
        //
        Mdouble sysBaseArea = 0.0;
        //
        switch(wallsBoundary)
        {
            case 1:
            {
                sysBaseArea = (ToolLength-ToolThickness) * (ToolWidth - (2.0*ToolThickness));
            }break;
            case 0:
            {
                sysBaseArea = getXMax()*getYMax();
            }
        }
        //
        Mdouble totalSysVolumeParticlesAndVoids = sysBaseArea * getSurfaceHeight();
        //
        Mdouble bulkDensity = massTotalParticles/totalSysVolumeParticlesAndVoids;
        Mdouble relativeDensity = bulkDensity/MatDensity;
        //
        //write header line
        if (eneFile.getCounter() == 1) os << "1time 2KineticEnergy 3massTotalParticles 4SurfaceHeight 5sysBaseArea 6totalSysVolumeParticlesAndVoids"
                                             "7volumeTotalParticles 8bulkDensity 9relativeDensity\n";
        //
        // ADD relative velocity?
        //
        os  << getTime() << '\t'
        << getKineticEnergy() << '\t'
        << massTotalParticles << '\t'
        << getSurfaceHeight() << '\t'
        << sysBaseArea << '\t'
        << totalSysVolumeParticlesAndVoids << '\t'
        << volumeTotalParticles << '\t'
        << bulkDensity << '\t'
        << relativeDensity << '\n';
    }
    //
    void writeFstatHeader(std::ostream& os) const override
    {
        //
        for (const auto &iH: interactionHandler) {
            auto i0 = dynamic_cast<MeltableInteraction *>(iH);
            logger.assert_debug(i0 != nullptr, "i0 not set");
            os << i0->getId() << " "
            << i0->getForce() << " "
            << i0->getRelativeVelocity() << " "
            << i0->getOverlap() << " "
            << i0->getSolidOverlap() << " "
            << i0->getBondingOverlap() << " "
            << i0->getContactRadiusMeltable() << " "
            << i0->getElasticForce() << " "
            << i0->getDampingForce() << " "
            << i0->getBondingForce() << " "
            << i0->getViscousForce() << " "
            << i0->getTensionForce() << " "
            << i0->getSolidContactRadius() << " "
            << i0->getBondingContactRadius() << std::endl;
        }

    }
    //
    void printTime() const override
    {
        logger(INFO,"Time % TimeMax % kineticEnergy %", getTime(), getTimeMax(),getKineticEnergy());
    }
    //
};
//
int main(int argc UNUSED, char *argv[] UNUSED)
{
    LaserBed problem;
    problem.solve();
    return 0;
}
