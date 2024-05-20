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



class ApplyLaserOnLayer: public Mercury3D
{
private:
    // set main simulation variables
    Mdouble timeStep = 1e-6;
    Mdouble scaleMass = 1e8;
    //
    unsigned int laserON = 1;
    //
    bool writeVTK = true;
    //
    unsigned int wallsBoundary = 1; // 0: PBC -> update laser penetration depth and
    //
    // ----------------------------------------------------------
    // each layer simulation time
    Mdouble setSimTime = 0.1; //3.0;
    //waiting time to avoid overlap, keep fixed?
    Mdouble waitingInterval = 0.1; //1.0;
    // max sim time for each layer
    Mdouble MaxSimTime = 0.1; //1.0; //setSimTime + waitingInterval;
    //
/*    Mdouble heatingDuration = 0.002; //0.006;
    Mdouble startHeatinTime = 0.001; //0.002;*/
    // waiting time interval between adding layers
    Mdouble waitingTime = 3.0; //setSimTime;
    Mdouble numberOfLayers = 5;
    //:::::Mdouble finalTime = waitingTime*numberOfLayers::::::
    //
    bool addingNewLayers = true;
    Mdouble minCond = 0.0; //setSimTime;
    Mdouble maxCond = MaxSimTime ;
    //
    //----------------------------------------------------------
    //
    int thermalExpansionOnOff = 0;
    //
    // control system size
    int spatialScale = 10; //2 //4
    Mdouble volumeScale = 2.5; //5; //3 //6
    // DPM parameters
    Mdouble gravityValue = -9.81/sqrt(scaleMass); //scaleMass;
    //
    Mdouble XMaxDomain = spatialScale*100e-6; //spatialScale*500e-6;
    Mdouble YMaxDomain = spatialScale*100e-6; //spatialScale*500e-6;
    Mdouble ZMinDomain = 120e-6; //400e-6;//200e-6;
    Mdouble ZMaxDomain = 100e-6;

    Mdouble ToolLength = spatialScale*100e-6; //spatialScale*500e-6; // unit m
    Mdouble ToolWidth = spatialScale*100e-6; //spatialScale*500e-6;
    Mdouble ToolHight = 0.01;//40e-1 * scale * 4.0; //20e-1 * scale * 4.0;
    Mdouble ToolThickness = 50e-6;//5e-1 * scale;
    Mdouble Gap = 30e-6; //10e-6; //30e-6;
    //
    Mdouble MatDensity = 1050 * scaleMass;
    Mdouble elasticModulus = 0.8*2.95e9; //8e7;
    Mdouble possionsRation = 0.4;
    Mdouble damping = 0.9; // * / sqrt(scaleMass);
    //
    Mdouble latentHeat = 102e3; //Sintratec PA12 -> 102e3 //TU/E -> 56.4e3;
    Mdouble heatCapacity = 2900; //Sintratec PA12 -> 2900 //TU/E ->1200;
    //

    // Test Effect?
    Mdouble surfaceTension = 34.3e-3 * sqrt(scaleMass); //* 1e-5
    Mdouble refVis = 1e-3 * sqrt(scaleMass); //2e-5

    // cond , conv, rad
    Mdouble thermalcond = 0.12 * scaleMass;
    Mdouble thermalConvCoff = 150.0 * scaleMass;
    Mdouble emmisivity = 0.9 * scaleMass;
    //
    Mdouble aTemp = 170 + 273.15;// -> 443.15 -- 140 + 273.15;// -> 413.15
    Mdouble particlesTemp =  aTemp;
    Mdouble mTemp = 185 + 273.15;//458.15

    Mdouble deltaT = 20;//+273.15;// (this value matched initially with Exp, but very high rate) //5.0;
    Mdouble liquidHeatCapacity = 3350; //1.2 * heatCapacity;//2.0*heatCapacity;

    Mdouble vaporizationHeatCapacity = 2.5 * heatCapacity;//4.0*heatCapacity;
    Mdouble vaporizationLatentHeat = 20.0 * latentHeat;
    Mdouble vaporizationTemp = 3.5 * mTemp;
    //
    Mdouble thermalExpansionCoeff = 90 * 1e-6; // K-1
    //
    bool laserOnOff = true;
    Mdouble absorptivity = 0.33 * scaleMass;
    Mdouble laserPower = 2.3; // (W) -> Sintratec Kit laser power
    Mdouble beamRadius = 250e-6;
    Mdouble powderBedPorosity = 0.55;
    Mdouble initialLaserPX =  -2.0*beamRadius; //0.0; //2.0*beamRadius; //(getXMax()+getXMin())/2.0 ;//beamRadius;
    Mdouble laserSpeedX = 0.0;
    Mdouble initialLaserPY = -2.0*beamRadius; // 0.0; //2.0*beamRadius; //(getYMax()+getXMin())/2.0; 2.0);
    Mdouble laserSpeedY = 0.0; //-0.5; //500 mm/s

    Mdouble laserPenetrationZ = -100e-6;
    int counter_ = 1;
    double laserZ_ = laserPenetrationZ;

    // ----------------------  Laser passes in Y-dir ----------------------
    int counterLaserPass = 1; //2;

    Mdouble hatchDistance = beamRadius;//2.0*meanRadius; //250e-6;

    Mdouble passDistance = YMaxDomain;
    Mdouble laserSpeed = 0.25; //0.1;  //unit: 0.5 m/s = 500 mm/s
    Mdouble timeLaserPass = passDistance/laserSpeed;

    bool activeLaser = true;
    Mdouble numberOfLaserPasses = (XMaxDomain/hatchDistance); //-1;
    //
    Mdouble timeLaserAfterFirstLayer = 0.0;
    // ------------------------------------------------------
    //
    unsigned int numberOfCounts = 1e3; //1e3;
    //
    //the volume to be added:
    Mdouble addVolume = volumeScale*Gap*ToolLength*ToolWidth;
    //particles:
    //
    // d50_p1 = 55.239016 ,d10_p1 = 23.218276, d90_p1 = 131.763261, dMax = 153.380932, dMin = 19.893925
    Mdouble meanRadius1 = 55e-6/2.0; //27e-6; // in m
    Mdouble MinRadius = 20e-6/2.0; //13.5e-6; //22.6e-6/2.0; //in m
    Mdouble MaxRadius = 153e-6/2.0; //72e-6;  //153.4/2.0; in m
    Mdouble LogmeanRadius = -10.4970; //-10.41;  //-10.3693;
    Mdouble LogstdRadius = 0.6446; //0.2758;
    //
    MeltableSpecies *species;
    //
    Mdouble heatingTime = 0.0; //startHeatinTime;
    Mdouble heatSource = 0.0 * scaleMass;
    Mdouble coolingTime = MaxSimTime; //heatingTime+heatingDuration;
    //
    // startFillHeight of new powder layer
    Mdouble startFillHeight = -120e-6;

public:
    //-------------------------------------------------- --------------------------------------------------
    ApplyLaserOnLayer() {
        //
        setName("InitialLayer");
        readRestartFile();
        setRestarted(false);
        //
        //removeParticles();
        //
        setName("ApplyLaserMultipleLayers");
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
        setTimeMax(setSimTime); //MaxSimTime);
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
        //species->setHeatingTime(heatingTime);
        //species->setCoolingTime(coolingTime);
        switch (laserON)
        {
            case 1:
            {
                // Laser
                species->setTurnLaserOnOff(laserOnOff);
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
                species->setHeatingTime(heatingTime);
                species->setCoolingTime(coolingTime);
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
        //
        if (getTime()>= 0.0 & getTime() < 1.0*getTimeStep())
        {
            species->setHighestParticleZForLaser(getSurfaceHeight());
            logger(INFO, " -firstLayer-: setting surface height for laser = %",species->getHighestParticleZForLaser());
        }
        //
        if (getTime() >= timeLaserPass + timeLaserAfterFirstLayer & getTime() < getTimeStep() + timeLaserPass + timeLaserAfterFirstLayer)
        {
            species->setHighestParticleZForLaser(getSurfaceHeight());
            logger(INFO, " -nextLayers-: setting surface height for laser = %",species->getHighestParticleZForLaser());
        }
        //
        // using one condition to add all layers:
        if ( addingNewLayers & (maxCond-waitingInterval) < getTime())
        {
            minCond = (counter_-1)*waitingTime; //MaxSimTime + waitingInterval + (counter_-1)*waitingTime;
            maxCond = counter_*waitingTime; //MaxSimTime + waitingInterval + counter_*waitingTime;
            setTimeMax(maxCond);
            logger(INFO, " setting new max time %, minCond %, maxCond %", getTimeMax(),minCond,maxCond);
            //
            timeLaserAfterFirstLayer = (minCond+maxCond)/2.0;
            //
            startFillHeight *= counter_;
            logger(INFO, "ADDING PARTICLES starting from Z = %", startFillHeight);
            addParticles();
            //
            laserZ_ -= 100e-6; //+= 100e-6;
            logger(INFO, "apply laser at new laser pass time = % and counter is rest to = %", timeLaserAfterFirstLayer + (timeLaserPass * counterLaserPass), counterLaserPass);
            //
            counter_+=1;
            logger(INFO, "next Layer %", counter_);
            //
            activeLaser = true;
        }
        else if (counter_ > numberOfLayers)
            addingNewLayers = false;
/*        else // this is for the first laser passes i.e. first layer
            applyLaser();*/
        applyLaser();
    }
    //
    void addParticles ()
    {
        Mdouble addVolumeNew = volumeScale*Gap*ToolLength*ToolWidth; //1*Gap*ToolLength*ToolWidth;
        //Mdouble startFillHeight = -150e-6;
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
        //
        species->setInitialLaserPositionX(-beamRadius); //0.0);
        species->setInitialLaserPositionY(-beamRadius); //0.0);
        //
        if (getTime() >= (timeLaserPass*counterLaserPass) + timeLaserAfterFirstLayer & getTime() < (timeLaserPass * counterLaserPass) + timeStep + timeLaserAfterFirstLayer & activeLaser)
        {
            //
            logger(INFO, "-------- Laser pass number = %,-------- Number of laser passes = % --------",
                   counterLaserPass,numberOfLaserPasses);
            //
            for (const auto &q: particleHandler)
            {
                auto p0 = dynamic_cast<ThermalParticle *>(q);
                logger.assert_debug(p0 != nullptr, "Thermal Particles required");
                //
                p0->setXLaser(counterLaserPass*hatchDistance);
                if (counterLaserPass%2 != 0)
                    p0->setYLaser(passDistance);
                else
                    p0->setYLaser(beamRadius); // same as initialYLaser
            }
            //
            species->setLaserSpeedX(0.0);
            if (counterLaserPass%2 != 0)
                species->setLaserSpeedY(-laserSpeed);
            else
                species->setLaserSpeedY(+laserSpeed);
            //
            counterLaserPass += 1;
            logger(INFO, "-------- Next pass = %",counterLaserPass);
            logger(INFO, "-------- Next time laser pass = %",timeLaserPass*counterLaserPass);
        }
        else if (counterLaserPass > numberOfLaserPasses)
        {
            logger(INFO, "-------- LASER OFF, laser pass time = % and counter = %", timeLaserPass,counterLaserPass);
            activeLaser = false;
            counterLaserPass = 1;
            //timeLaserAfterFirstLayer = maxCond;
            species->setLaserPenetrationZCoor(laserZ_);
            logger(INFO, "-------- LASER OFF,  new laser pass time = % , reset counter = % and penetration coordinate Z = %",
                   timeLaserAfterFirstLayer + (timeLaserPass * counterLaserPass), counterLaserPass, species->getLaserPenetrationZCoor());
        }
        //
/*        Mdouble newHeatingTime = startHeatinTime+minCond; //getTimeMax();
        Mdouble newCoolingTime = newHeatingTime+heatingDuration;
        logger(INFO, "Applying Laser at Zplan %, laserOn at time %, laserOff at time %", laserZ_, newHeatingTime, newCoolingTime);
        species->setHeatingTime(newHeatingTime);
        species->setCoolingTime(newCoolingTime);
        species->setLaserPenetrationZCoor(laserZ_);*/
        //
    }

    /*
     * Returns the height of the highest particle
     */
    Mdouble getSurfaceHeight() const
    {
        double height = -ZMinDomain; //0;
        for (const BaseParticle* p : particleHandler)
        {
            double newHeight = p->getPosition().Z;// + p->getRadius();
            if (height<newHeight) height = newHeight;
        }
        return height;//+ZMinDomain;
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
    ApplyLaserOnLayer problem;
    problem.solve();
    return 0;
}
