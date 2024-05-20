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
    Mdouble timeStep = 5e-6; //2e-7; // 1e-6;
    Mdouble scaleMass = 1e8; // 1e6;
    //
    unsigned int laserON = 1;
    Mdouble laserPenetrationZ = -100e-6;
    //
    bool writeVTK = true;
    //
    unsigned int wallsBoundary = 1; // 0: PBC -> update laser penetration depth and
    //
    // each layer simulation time
    Mdouble setSimTime = 8.0; //30.0; //0.5; //1.0;
    //waiting time to avoid overlap, keep fixed?
    Mdouble waitingInterval = 2.0; //0.0001;
    // max sim time for each layer
    Mdouble MaxSimTime = setSimTime + waitingInterval;
    //
    Mdouble heatingDuration = 0.002; //0.006;
    Mdouble startHeatinTime = 0.001; //0.002;
    // waiting time interval between adding layers
    Mdouble waitingTime = setSimTime; //0.5; //1.0; //30;
    Mdouble numberOfLayers = 3;//10.0;
    //:::::Mdouble finalTime = waitingTime*numberOfLayers::::::
    //
    bool addingNewLayers = true;
    Mdouble minCond = setSimTime;
    Mdouble maxCond = MaxSimTime ;//getTimeMax();
    //
    int thermalExpansionOnOff = 0;
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
    Mdouble heatingTime = startHeatinTime;
    Mdouble heatSource = 0.0 * scaleMass;
    Mdouble coolingTime = heatingTime+heatingDuration;
    //
    Mdouble absorptivity = 0.33 * scaleMass;
    Mdouble laserPower = 8.0; //16.0
    Mdouble beamRadius = 250e-6;
    Mdouble powderBedPorosity = 0.55;
    Mdouble initialLaserPX = beamRadius;//-beamRadius;
    Mdouble laserSpeedX = 0.0; //0.1;//0.0;//-0.1
    Mdouble initialLaserPY = beamRadius;//-beamRadius;//beamRadius * 2.0 + (beamRadius / 2.0); // second pass -> beamRadius*2.0;// first scan -> beamRadius;
    Mdouble laserSpeedY = 0.0; //0.1;//0.0;//0.1;
    //
    unsigned int numberOfCounts = 1e3; //1e4;
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
    int counter_ = 1;
    double laserZ_ = laserPenetrationZ;
public:
    //-------------------------------------------------- --------------------------------------------------
    LaserBed() {
        //
        setName("InitialBed");
        readRestartFile();
        setRestarted(false);
        //
        //removeParticles();
        //
        setName("LaserBed");
        //
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
                species->setLaserPenetrationZCoor(laserPenetrationZ);
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
        // HERE
        //
    };
    //
    void actionsAfterTimeStep() override
    {
        // using one condition to add all layers:
        if ( addingNewLayers & (maxCond-waitingInterval) < getTime()) //& getTime() < maxCond
        {
            minCond = MaxSimTime + waitingInterval + (counter_-1)*waitingTime;
            maxCond = MaxSimTime + waitingInterval + counter_*waitingTime;
            setTimeMax(maxCond);
            logger(INFO, " setting new max time %, minCond %, maxCond %", getTimeMax(),minCond,maxCond);
            //
            logger(INFO, "ADDING PARTICLES");
            addParticles();
            //
            applyLaser();
            laserZ_ += 100e-6; //75.0e-6; //50.0e-6;
            //
            counter_+=1;
            logger(INFO, "next Layer %", counter_);
            //
        }
        else if (counter_ > numberOfLayers)
            addingNewLayers = false;
    }
    //
    void addParticles ()
    {
        Mdouble addVolumeNew = 1*Gap*ToolLength*ToolWidth;
        Mdouble startFillHeight = -150e-6;
        //
        std::mt19937 gen;
        gen.seed(100);
        std::lognormal_distribution<> d(LogmeanRadius, LogstdRadius);
        ThermalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(meanRadius1);
        p.setTemperature(particlesTemp);
        Mdouble fillHeight = startFillHeight;
        while (addVolumeNew > 0) {
            Mdouble x = random.getRandomNumber(ToolThickness, ToolLength);
            Mdouble y = random.getRandomNumber(ToolThickness, (ToolWidth - ToolThickness));
            Mdouble z = random.getRandomNumber(startFillHeight, fillHeight);
            p.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p)) {
                particleHandler.copyAndAddObject(p);
                addVolumeNew -= p.getVolume();
                do {
                    p.setRadius(d(gen));
                } while (p.getRadius() < MinRadius || p.getRadius() > MaxRadius);
                if (particleHandler.getNumberOfObjects() % 100 == 0) std::cout << '.' << std::flush;
                if (particleHandler.getNumberOfObjects() % 1000 == 0) std::cout << ' ';
                if (particleHandler.getNumberOfObjects() % 10000 == 0) std::cout << addVolumeNew << '\n';
            } else {
                fillHeight += 10e-6 *meanRadius1;
            }
        }
        logger(INFO, " Inserted % particles", particleHandler.getNumberOfObjects());
        //
    }
    //
    void applyLaser()
    {
        Mdouble newHeatingTime = startHeatinTime+minCond; //getTimeMax();
        Mdouble newCoolingTime = newHeatingTime+heatingDuration;
        logger(INFO, "Applying Laser at Zplan %, laserOn at time %, laserOff at time %", laserZ_, newHeatingTime, newCoolingTime);
        species->setHeatingTime(newHeatingTime);
        species->setCoolingTime(newCoolingTime);
        species->setLaserPenetrationZCoor(laserZ_);
        //
/*        Mdouble newLaserXPosition = -beamRadius*counter_;
        Mdouble newLaserYPosition = -beamRadius*counter_;
        species->setInitialLaserPositionX(newLaserXPosition);
        species->setInitialLaserPositionY(newLaserYPosition);*/
        //OR: re set initial laser X and Y, by setXLaser and setYLaser to 0.0 , add function in Thermal Particle and add new species Var for that
    }

    /*
     * Returns the height of the highest particle
     */
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
};
//
int main(int argc UNUSED, char *argv[] UNUSED)
{
    LaserBed problem;
    problem.solve();
    return 0;
}
