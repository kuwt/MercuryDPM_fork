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
#include <Walls/InfiniteWall.h>
#include <Walls/IntersectionOfWalls.h>
#include "DPMBase.h"


class TestLaser: public Mercury3D
{
private:
    // set main simulation variables
    Mdouble timeStep = 1e-6; //2.9e-4 -> MS1e8 //4.6e-5 -> MS1e6 //7.4e-6 -> MS1e4 //1.8e-7 -> MS1
    Mdouble scaleMass = 1;
    //
    unsigned int laserON = 1;
    //
    unsigned int wallsBoundary = 2;
    //
    bool writeVTK = true;
    Mdouble maxSimTime = 0.05; //0.2; //2.0; //20.0;
    unsigned int numberOfCounts = 1e3;//e4;
    //
    //particles:
    //
    Mdouble meanRadius = 50e-6; //27e-6;// in m
    Mdouble MinRadius = 20e-6;//22.6e-6/2.0;  //19.5e-6/2.0;//in m
    Mdouble MaxRadius = 75e-6;//153.4/2.0;    //174.5e-6/2.0;//in m
    //
    int numberOfParticlesInX = 5;
    int numberOfParticlesInY = 4; //10;
    int numberOfParticlesInZ = 5;
    Mdouble x = meanRadius;
    Mdouble z = meanRadius;
    Mdouble overlap = 0.0;
    //
    //Mdouble heatingTime = 0.0;
    //Mdouble heatSource = 0.0 * scaleMass;
    //Mdouble coolingTime = 2.0;

    Mdouble absorptivity = 0.33;
    Mdouble laserPower = 5;//1;
    Mdouble beamRadius = 100e-6; // 250e-6; //125e-6;
    Mdouble powderBedPorosity = 0.55;
    Mdouble initialLaserPX = beamRadius; //-beamRadius;
    Mdouble laserSpeedX = 0.0; //0.1;//-0.1;//100e-6;
    Mdouble initialLaserPY = beamRadius;//-beamRadius;//beamRadius * 2.0 + (beamRadius / 2.0); // second pass -> beamRadius*2.0;// first scan -> beamRadius;
    Mdouble laserSpeedY = 0.1; //0.0; //0.1;

    Mdouble laserPenetrationZ = 250e-6; //-250e-6;

    Mdouble heatingTime = 0.0;
    Mdouble heatSource = 0.0 * scaleMass;
    Mdouble coolingTime = maxSimTime; //2.0;
    //
    int thermalExpansionOnOff = 0;
    //
    // DPM parameters
    Mdouble gravityAcc = -9.81/sqrt(scaleMass); //sqrt(scaleMass); //scaleMass; //sqrt(scaleMass);
    //
    //
    // control system size
    int spatialScale = 2;//4
    int volumeScale = 6; //3 //6
    //
    Mdouble XMaxDomain = numberOfParticlesInX*2.0 *meanRadius + 2.0*meanRadius; //spatialScale*500e-6/2;
    Mdouble YMaxDomain = numberOfParticlesInY*2.0 *meanRadius + 2.0*meanRadius; //spatialScale*500e-6/2;
    Mdouble ZMinDomain = 550e-6;//200e-6;
    Mdouble ZMaxDomain = 100e-6;
    // Domain
    Mdouble ToolLength = spatialScale*500e-6/2;//0.00025; //45e-1 * scale * 2.0; // unit m
    Mdouble ToolWidth = spatialScale*500e-6/2;//0.00025;//30e-1 * scale * 2.0;
    Mdouble ToolHight = 0.01;//40e-1 * scale * 4.0; //20e-1 * scale * 4.0;
    Mdouble ToolThickness = 50e-6;//5e-1 * scale;
    Mdouble Gap = 30e-6;
    //
    Mdouble MatDensity = 1050 * scaleMass;
    Mdouble elasticModulus = 1e5; //0.8*2.95e9; //8e7;
    Mdouble possionsRation = 0.4;
    Mdouble damping = 0.8; // * / sqrt(scaleMass);
    //
    Mdouble latentHeat = 56.4e3;
    Mdouble heatCapacity = 1200;
    //
    Mdouble surfaceTension = 0.0;// 34.3e-3 * sqrt(scaleMass);
    Mdouble refVis = 0.0;//2e-3 * sqrt(scaleMass);
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
    Mdouble particleInitialTemp = 155 + 273.15;//413.15;
    //
    //
    MeltableSpecies *species;
    //

    // ----------------------  Laser passes in Y-dir ----------------------
    int counterLaserPass = 2;

    Mdouble hatchDistance = beamRadius;//2.0*meanRadius; //250e-6;

    Mdouble passDistance = YMaxDomain;
    Mdouble laserSpeed = 0.1;
    Mdouble timeLaserPass = passDistance/laserSpeed;

    bool activeLaser = true;
    Mdouble numberOfLaserPasses = (XMaxDomain/hatchDistance)-1;
    // ------------------------------------------------------
    //
public:
    //
    void setupInitialConditions() override {
        //
        setName("TestLaser");
        setSystemDimensions(3);
        setGravity(Vec3D(0.0, 0.0, gravityAcc));
        setZMin(-ZMinDomain);
        setXMax(XMaxDomain);
        setYMax(YMaxDomain);
        setZMax(ZMaxDomain);
        setTimeMax(maxSimTime);
        //
        //auto species = speciesHandler.copyAndAddObject(MeltableSpecies());
        species = speciesHandler.copyAndAddObject(MeltableSpecies());
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
        //
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
        //species->setMaxTimeStep(MinRadius,MatDensity,elasticModulus);
        setTimeStep(timeStep);
        //
        //restartFile.setSaveCount(restartSaveCount);
        //fStatFile.setSaveCount(10000);
        dataFile.setFileType(FileType::ONE_FILE);
        //
        //restartFile.writeFirstAndLastTimeStep();
        restartFile.setFileType(FileType::ONE_FILE);
        //
        fStatFile.setFileType(FileType::ONE_FILE);
        //eneFile.setFileType(FileType::ONE_FILE);//NO_FILE);
        //Data Output
        setSaveCount(numberOfCounts);
        setParticlesWriteVTK(writeVTK);
        setWallsWriteVTK(writeVTK);
        //setWallsWriteVTK(FileType::MULTIPLE_FILES);
        //
        //Bottom Wall
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0)); //w0.setSpecies(speciesHandler.getObject(1));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, -meanRadius)); //-ZMinDomain));
        wallHandler.copyAndAddObject(w0);
        //
        switch(wallsBoundary) {
            case 1:
            {
                //
                Mdouble damThickness = ToolLength + ToolThickness;//(5.0e-1 * scale);
                //
                //side wall
                Vec3D minPoint1_tool = Vec3D(0.0, 0.0, -ZMinDomain);
                Vec3D maxPoint1_tool = Vec3D(ToolLength, ToolThickness, ToolHight); //Vec3D(getXMax(),0.05*getYMax(),getZMax())
                //back wall
                Vec3D minPoint2_tool = Vec3D(0.0, 0.0, -ZMinDomain); //Gap);
                Vec3D maxPoint2_tool = Vec3D(ToolThickness, ToolWidth, ToolHight); //Vec3D(ToolThickness,(ToolWidth-ToolThickness),ToolHight);
                //side wall
                Vec3D minPoint3_tool = Vec3D(0.0, (ToolWidth - ToolThickness), -ZMinDomain);
                Vec3D maxPoint3_tool = Vec3D(ToolLength, ToolWidth, ToolHight); // Vec3D(getXMax(),getYMax(),getZMax())
                //
                //Front Dam
                Vec3D minPoint_dam = Vec3D(ToolLength, 0.0, -ZMinDomain);
                Vec3D maxPoint_dam = Vec3D(damThickness, ToolWidth, ToolHight);
                //
                Mdouble topWallZ = 100e-6;
                Vec3D minPointTop = Vec3D(0.0,0.0,topWallZ);
                Vec3D maxPointTop = Vec3D(ToolLength,ToolWidth,topWallZ+ToolThickness);
                //
                IntersectionOfWalls w1;
                w1.setSpecies(speciesHandler.getObject(0));
                w1.addObject(Vec3D(1.0, 0.0, 0.0), minPoint1_tool);
                w1.addObject(Vec3D(0.0, 1.0, 0.0), minPoint1_tool);
                w1.addObject(Vec3D(0.0, 0.0, 1.0), minPoint1_tool);
                w1.addObject(Vec3D(-1.0, 0.0, 0.0), maxPoint1_tool);
                w1.addObject(Vec3D(0.0, -1.0, 0.0), maxPoint1_tool);
                w1.addObject(Vec3D(0.0, 0.0, -1.0), maxPoint1_tool);
                wallHandler.copyAndAddObject(w1);
                //
                IntersectionOfWalls w2;
                w2.setSpecies(speciesHandler.getObject(0));
                w2.addObject(Vec3D(1.0, 0.0, 0.0), minPoint2_tool);
                w2.addObject(Vec3D(0.0, 1.0, 0.0), minPoint2_tool);
                w2.addObject(Vec3D(0.0, 0.0, 1.0), minPoint2_tool);
                w2.addObject(Vec3D(-1.0, 0.0, 0.0), maxPoint2_tool);
                w2.addObject(Vec3D(0.0, -1.0, 0.0), maxPoint2_tool);
                w2.addObject(Vec3D(0.0, 0.0, -1.0), maxPoint2_tool);
                wallHandler.copyAndAddObject(w2);
                //
                IntersectionOfWalls w3;
                w3.setSpecies(speciesHandler.getObject(0));
                w3.addObject(Vec3D(1.0, 0.0, 0.0), minPoint3_tool);
                w3.addObject(Vec3D(0.0, 1.0, 0.0), minPoint3_tool);
                w3.addObject(Vec3D(0.0, 0.0, 1.0), minPoint3_tool);
                w3.addObject(Vec3D(-1.0, 0.0, 0.0), maxPoint3_tool);
                w3.addObject(Vec3D(0.0, -1.0, 0.0), maxPoint3_tool);
                w3.addObject(Vec3D(0.0, 0.0, -1.0), maxPoint3_tool);
                wallHandler.copyAndAddObject(w3);
                //
                IntersectionOfWalls w4;
                w4.setSpecies(speciesHandler.getObject(0));
                w4.addObject(Vec3D(1.0, 0.0, 0.0), minPoint_dam);
                w4.addObject(Vec3D(0.0, 1.0, 0.0), minPoint_dam);
                w4.addObject(Vec3D(0.0, 0.0, 1.0), minPoint_dam);
                w4.addObject(Vec3D(-1.0, 0.0, 0.0), maxPoint_dam);
                w4.addObject(Vec3D(0.0, -1.0, 0.0), maxPoint_dam);
                w4.addObject(Vec3D(0.0, 0.0, -1.0), maxPoint_dam);
                wallHandler.copyAndAddObject(w4);
                //
                /*                IntersectionOfWalls w5;
                                w5.setSpecies(speciesHandler.getObject(0));
                                w5.addObject(Vec3D(1.0, 0.0, 0.0), minPointTop);
                                w5.addObject(Vec3D(0.0, 1.0, 0.0), minPointTop);
                                w5.addObject(Vec3D(0.0, 0.0, 1.0), minPointTop);
                                w5.addObject(Vec3D(-1.0, 0.0, 0.0), maxPointTop);
                                w5.addObject(Vec3D(0.0, -1.0, 0.0), maxPointTop);
                                w5.addObject(Vec3D(0.0, 0.0, -1.0), maxPointTop);
                                wallHandler.copyAndAddObject(w5);*/
                //
            }break;
                //
            case 0:
            {
                PeriodicBoundary by1;
                by1.set(Vec3D(0, 1, 0), 0.0, ToolWidth);
                boundaryHandler.copyAndAddObject(by1);
                PeriodicBoundary bx1;
                bx1.set(Vec3D(1, 0, 0), 0.0, ToolLength);
                boundaryHandler.copyAndAddObject(bx1);
            }break;
            case 2:
            {
                logger(INFO, "No walls, no PBC");
            }
        }
        //
/*                ThermalParticle p0;
                p0.setSpecies(speciesHandler.getObject(0));
                p0.setTemperature(particleInitialTemp);
                p0.setRadius(pRadius1); //particle-1 radius 50 um
                p0.setPosition(Vec3D(pRadius1, 0.0, pRadius1));
                p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
                particleHandler.copyAndAddObject(p0); // 1st particle created
                //
                Mdouble overlap = 1.0;
                //
                ThermalParticle p1;
                p1.setSpecies(speciesHandler.getObject(0));
                p1.setTemperature(particleInitialTemp);
                p1.setRadius(pRadius2); // particle-2 radius
                p1.setPosition(Vec3D(pRadius1,(pRadius1*2.0)-(overlap*pRadius1),pRadius2));
                p1.setVelocity(Vec3D(0.0, 0.0, 0.0));
                particleHandler.copyAndAddObject(p1); // 2nd particle created*/
        //
        // creat monodisperse particles in a grid
        ThermalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setTemperature(particleInitialTemp);
        p0.setRadius(meanRadius);
        p0.setVelocity(Vec3D(0.0, 0.0, 0.0));
        for (int i=0; i <= numberOfParticlesInY;i++) {
            for (int j=0; j <= numberOfParticlesInX;j++) {
                for (int k=0; k <= numberOfParticlesInZ; k++) {
                    logger(INFO, " Inserting particles %", particleHandler.getNumberOfObjects());
                    //p0.setRadius(random.getRandomNumber(MinParticleRadius, MaxParticleRadius));
                    //p0.setPosition(Vec3D(j*2.0*meanRadius-(overlap*meanRadius), i*2.0*meanRadius-(overlap*meanRadius),k*2.0*meanRadius));
                    p0.setPosition(Vec3D(j*2.0*meanRadius+meanRadius, i*2.0*meanRadius+meanRadius,k*2.0*meanRadius));
                    particleHandler.copyAndAddObject(p0);
                }
            }
        }
//
    }
    //
    void printTime() const override
    {
        logger(INFO,"Time % TimeMax %", getTime(), getTimeMax());
    }
    //
    void actionsAfterTimeStep() override
    {
        //
        species->setInitialLaserPositionX(0.0);
        species->setInitialLaserPositionY(0.0);
        //
        if (getTime() >= timeLaserPass & getTime() < timeLaserPass+timeStep & activeLaser)
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
            logger(INFO, "-------- time laser pass before = %",timeLaserPass);
            timeLaserPass += timeLaserPass;
            logger(INFO, "-------- time laser pass AFTER = %",timeLaserPass);

            counterLaserPass += 1;
            logger(INFO, "-------- Next pass = %",counterLaserPass);
        }
        else if (counterLaserPass > numberOfLaserPasses & getTime() >= timeLaserPass & getTime() < timeLaserPass+timeStep )
        {
            activeLaser = false;
            logger(INFO, "---------------------- Finished laser passes ---------------------------------");
        }
        // ---------------------------------------------------------------------------------
       /* // SECOND LASER PASS
        Mdouble initialLaserPX2 = getXMin();
        Mdouble laserSpeedX2 = 0.1;
        Mdouble initialLaserPY2 = getYMax()/2.0;
        Mdouble laserSpeedY2 = 0.0;

        Mdouble timeLaserPass2 = initialLaserPY2/laserSpeedY; //0.01 (s)
        //
        // THIRD LASER PASS
        Mdouble initialLaserPX3 = ((getXMax()+getXMin())/2.0)+meanRadius;
        Mdouble laserSpeedX3 = 0.0;
        Mdouble initialLaserPY3 = getYMax()/2.0;
        Mdouble laserSpeedY3 = -0.1;

        Mdouble timeLaserPass3 = timeLaserPass2 + (initialLaserPX3/laserSpeedX2);
        //
        //Mdouble hatchDis = 250e-6;
        //Mdouble firstXCoor = sqrt(2) * hatchDis;
        // odd
        //std::vector<double> xCoor = {firstXCoor,3.0*firstXCoor,5.0*firstXCoor,7.0*firstXCoor,0.2,0.1,0.05};
        //std::vector<double> yCoor = {0.4,0.3,0.2,0.1,0.05,0.005};
        // even

        // 4th LASER PASS
        Mdouble hatchDistance = meanRadius*3.0;

        Mdouble initialLaserPX4 = sqrt(2) * hatchDistance;
        Mdouble laserSpeedX4 = -0.1;
        Mdouble initialLaserPY4 = 0.0;
        Mdouble laserSpeedY4 = 0.1;

        Mdouble timeLaserPass4 = timeLaserPass3 + ((getYMax()+getYMin())/2.0)/0.1; //100*timeStep; //(initialLaserPX4/(-laserSpeedX3));
        //
        // 5th LASER PASS

        Mdouble initialLaserPX5 = 0.0;
        Mdouble laserSpeedX5 = +0.1;
        Mdouble initialLaserPY5 = 2.0 * sqrt(2) * hatchDistance;
        Mdouble laserSpeedY5 = -0.1;

        Mdouble timeLaserPass5 = timeLaserPass4 + ((getXMax()+getXMin())/2.0)/0.1;
        //
        //
        //
        // ------- 2
        if (getTime() >= timeLaserPass2 && getTime() < timeLaserPass2+timeStep)
        {
            //
            for (const auto &q: particleHandler)
            {
                auto p0 = dynamic_cast<ThermalParticle *>(q);
                logger.assert_debug(p0 != nullptr, "Thermal Particles required");
                //
                p0->setXLaser(initialLaserPX2);
                p0->setYLaser(initialLaserPY2);
            }
            //
            species->setLaserSpeedX(laserSpeedX2);
            species->setLaserSpeedY(laserSpeedY2);
            //
        }
        //
        // ------- 3
        else if (getTime() >= timeLaserPass3 && getTime() < timeLaserPass3+timeStep)
        {
            //
            for (const auto &q: particleHandler)
            {
                auto p0 = dynamic_cast<ThermalParticle *>(q);
                logger.assert_debug(p0 != nullptr, "Thermal Particles required");
                //
                p0->setXLaser(initialLaserPX3);
                p0->setYLaser(initialLaserPY3);
            }
            //
            species->setLaserSpeedX(laserSpeedX3);
            species->setLaserSpeedY(laserSpeedY3);
        }
        //
        // ------- 4
        else if (getTime() >= timeLaserPass4 && getTime() < timeLaserPass4+timeStep)
        {
            //
            for (const auto &q: particleHandler)
            {
                auto p0 = dynamic_cast<ThermalParticle *>(q);
                logger.assert_debug(p0 != nullptr, "Thermal Particles required");
                //
                p0->setXLaser(initialLaserPX4);
                p0->setYLaser(initialLaserPY4);
            }
            //
            species->setLaserSpeedX(laserSpeedX4);
            species->setLaserSpeedY(laserSpeedY4);
        }
        //
        // ------- 5
        else if (getTime() >= timeLaserPass5 && getTime() < timeLaserPass5+timeStep)
        {
            //
            for (const auto &q: particleHandler)
            {
                auto p0 = dynamic_cast<ThermalParticle *>(q);
                logger.assert_debug(p0 != nullptr, "Thermal Particles required");
                //
                p0->setXLaser(initialLaserPX5);
                p0->setYLaser(initialLaserPY5);
            }
            //
            species->setLaserSpeedX(laserSpeedX5);
            species->setLaserSpeedY(laserSpeedY5);
        }*/
        //
    }
    //
    // Write custom output to the ene file
    void writeEneTimeStep(std::ostream &os) const override
    {
        //write header line
        if (eneFile.getCounter() == 1) os << "time relativeVelocity tensionForce visForce bondingForce dampingForce elasticForce "
                                             "overlap solidOverlap bondingOverlap meltRateP meltRateI "
                                             "contactRadius solidContactRadius bondingContactRadius "
                                             "visCoeff visTimeStep setTimeStep Temp MeltLayerThickness\n";
        //write interaction properties at every time step
        auto i = interactionHandler.getObject(0); //getLastObject(); // getObject(0);
        logger.assert_debug(i,"No interaction exists");
        auto mi = dynamic_cast<const MeltableInteraction *>(i);
        logger.assert_debug(mi,"Interaction is not of type MeltableInteraction");
        //
        auto p = particleHandler.getLastObject();
        logger.assert_debug(p,"No particle");
        auto pi = dynamic_cast<const ThermalParticle *>(p);
        logger.assert_debug(pi,"particle is not of type ThermalParticle");
        //
        os <<  getTime() << ' '
           << mi->getNormalRelativeVelocity() << ' '
           << mi->getTensionForce() << ' '
           << mi->getViscousForce() << ' '
           << mi->getBondingForce() << ' '
           << mi->getDampingForce() << ' '
           << mi->getElasticForce() << ' '
           << mi->getOverlap() << ' '
           << mi->getSolidOverlap() << ' '
           << mi->getBondingOverlap() << ' '
           << mi->getMoltenLayerThicknessRateP() << ' '
           << mi->getMoltenLayerThicknessRateI() << ' '
           << mi->getContactRadiusMeltable() << ' '
           << mi->getSolidContactRadius() << ' '
           << mi->getBondingContactRadius() << ' '
           << mi->getVisCoeff() << ' '
           //<< mi->getTimeStepVis() << ' '
           << mi->getHandler()->getDPMBase()->getTimeStep() << ' '
           << pi->getTemperature() << ' '
           << pi-> getMoltenLayerThickness() << ' '
           << pi->getMass() << ' '
           << pi->getInvMass() << ' '
           << helpers::getEffectiveMass(pi->getMass(),pi->getMass()) << ' '
           << pi->getRadius() << std::endl;
    }
};


int main(int argc UNUSED, char* argv[] UNUSED)
{
    TestLaser problem;
    problem.setXBallsAdditionalArguments("-solidf -v0 -s .85");
    problem.solve();
    return 0;
}
