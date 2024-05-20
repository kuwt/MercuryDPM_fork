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
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include <random>


class InitialBed: public Mercury3D
        {
        private:
            // set main simulation variables
            Mdouble timeStep = 5e-5; //2e-5; //2.9e-4 -> MS1e8 //4.6e-5 -> MS1e6 //7.4e-6 -> MS1e4 //1.8e-7 -> MS1
            Mdouble scaleMass = 1e8;
            //
            unsigned int laserON = 0;
            //
            unsigned int wallsBoundary = 1; //0->PBC , 1->wall
            //
            bool writeVTK = true;

            Mdouble MaxSimTime = 10.0 + 1.0; //2.0; //20.0;
            Mdouble timeToRemoveOutofBondParticles = 10.5;

            unsigned int numberOfCounts = 1e3;//e1;
            //
            unsigned int randomParticlesPacking = 0;
            unsigned int seedNo = 100;
            //
            int thermalExpansionOnOff = 0;
            //
            // control system size
            int spatialScale = 10; //2;//4
            Mdouble volumeScale = 2.5; //2.5 -> height 100 um //5 -> ;
            // DPM parameters
            Mdouble gravityValue = -9.81/sqrt(scaleMass);
            // Domain:
            Mdouble XMaxDomain = spatialScale*100e-6; //spatialScale*500e-6;
            Mdouble YMaxDomain = spatialScale*100e-6; //spatialScale*500e-6;
            Mdouble ZMinDomain = 120e-6;//200e-6;
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
            Mdouble damping = 0.9;
            //
            Mdouble latentHeat = 102e3; //Sintratec PA12 -> 102e3 //TU/E -> 56.4e3;
            Mdouble heatCapacity = 2900; //Sintratec PA12 -> 2900 //TU/E ->1200;
            //
            Mdouble surfaceTension =  34.3e-3 * 1e-5 * sqrt(scaleMass); //<- from HSM Cal. //34.3e-3 * sqrt(scaleMass);
            Mdouble refVis = 2e-5 * sqrt(scaleMass); //<- from HSM Cal. //2e-3 * sqrt(scaleMass);
            // cond , conv, rad
            Mdouble thermalcond = 0.12 * scaleMass;
            Mdouble thermalConvCoff = 150.0 * scaleMass;
            Mdouble emmisivity = 0.9 * scaleMass;
            //
            Mdouble aTemp = 170 + 273.15; //<- SLS  // TU/E -> 155 + 273.15;//428.15
            Mdouble mTemp = 185 + 273.15; //<- SLS  // TU/E -> 178 + 273.15;//451.15;

            Mdouble deltaT = 20;//+273.15;// (this value matched initially with Exp, but very high rate) //5.0;
            Mdouble liquidHeatCapacity = 3350; //1.2 * heatCapacity;//2.0*heatCapacity;

            Mdouble vaporizationHeatCapacity = 2.5 * heatCapacity;//4.0*heatCapacity;
            Mdouble vaporizationLatentHeat = 20.0 * latentHeat;
            Mdouble vaporizationTemp = 3.5 * mTemp;
            //
            Mdouble thermalExpansionCoeff = 90 * 1e-6; // K-1
            //
            Mdouble particlesTemp = aTemp; //<- 170 + 273.15; //155 + 273.15 = 413.15;
            //
            Mdouble heatingTime = 0.0;
            Mdouble heatSource = 0.0 * scaleMass;
            Mdouble coolingTime = 0.0;
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
        public:
            void setupInitialConditions() override
            {
                //
                setName("InitialLayer");
                setSystemDimensions(3);
                setGravity(Vec3D(0.0,0.0,gravityValue));
                setZMin(-ZMinDomain);
                setXMax(XMaxDomain);
                setYMax(YMaxDomain);
                setZMax(ZMaxDomain);
                setTimeMax(MaxSimTime);
                //
                auto species = speciesHandler.copyAndAddObject(MeltableSpecies());
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
                //Data Output
                setSaveCount(numberOfCounts);
                setParticlesWriteVTK(writeVTK);
                setWallsWriteVTK(writeVTK);
                //
                //Bottom Wall
                InfiniteWall w0;
                w0.setSpecies(speciesHandler.getObject(0)); //w0.setSpecies(speciesHandler.getObject(1));
                w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, -ZMinDomain));
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
/*                        Mdouble topWallZ = 100e-6;
                        Vec3D minPointTop = Vec3D(0.0,0.0,topWallZ);
                        Vec3D maxPointTop = Vec3D(ToolLength,ToolWidth,topWallZ+ToolThickness);*/
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
                    }
                }
                //
                //TOP Wall
                /*        InfiniteWall wTop;
                        wTop.setSpecies(speciesHandler.getObject(0)); //w0.setSpecies(speciesHandler.getObject(1));
                        wTop.set(Vec3D(0.0, 0.0, 1.0), Vec3D(0.0, 0.0, 0.0));
                        wallHandler.copyAndAddObject(wTop);*/
                //
                // Set up random number generator
                //
                /*        std::mt19937 gen;
                        switch (randomParticlesPacking) {
                            case 0:
                            {
                                //std::mt19937 gen;
                                gen.seed(seedNo);
                            }break;
                            case 1:
                            {
                                std::random_device rd;
                                std::mt19937 gen(rd());
                            }

                        }*/
                //
                std::mt19937 gen;
                gen.seed(seedNo);
                std::lognormal_distribution<> d(LogmeanRadius, LogstdRadius);
                //add particles until the volume to be added is zero
                //logger(INFO,"Adding particles ...");
                ThermalParticle p;
                p.setSpecies(speciesHandler.getObject(0));
                p.setRadius(meanRadius1);
                p.setTemperature(particlesTemp);
                Mdouble fillHeight = -ZMinDomain; //Gap; //0.0;
                while (addVolume > 0) {
                    Mdouble x = random.getRandomNumber(ToolThickness, ToolLength);
                    Mdouble y = random.getRandomNumber(ToolThickness, (ToolWidth-ToolThickness));
                    Mdouble z = random.getRandomNumber(-ZMinDomain, fillHeight); //Gap, fillHeight);
                    p.setPosition({x, y, z});
                    // check if particle can be inserted
                    if (checkParticleForInteraction(p)) {
                        particleHandler.copyAndAddObject(p);
                        addVolume -= p.getVolume();
                        do {
                            p.setRadius(d(gen));
                        } while (p.getRadius() < MinRadius || p.getRadius() > MaxRadius);
                        if (particleHandler.getNumberOfObjects() % 10 == 0) std::cout << '.' << std::flush;
                        if (particleHandler.getNumberOfObjects() % 100 == 0) std::cout << ' ';
                        if (particleHandler.getNumberOfObjects() % 1000 == 0) std::cout << addVolume << '\n';
                    } else {
                        fillHeight += 5e-12; //1e-11; //1e-12;
                        //30e-6*meanRadius1;//0.01 * scale * meanRadius1; //increase fill height (slowly to insert particles as low as possible)
                    }
                }
                logger(INFO, " Inserted % particles", particleHandler.getNumberOfObjects());
            }
            //
            void actionsAfterTimeStep() override
            {
                if (getTime()>= timeToRemoveOutofBondParticles) {
                    for (int i = 0; i < particleHandler.getNumberOfObjects(); ++i) {
                        Vec3D particlePosition = particleHandler.getObject(i)->getPosition();
                        if (particlePosition.Z > 0.0) {
                            particleHandler.removeObject(i);
                        }
                    }
                }
            }
            //
            //override continueSolve function such that the code stops when the packing is relaxed (Ekin<1e-16)
            /*   bool continueSolve() const override
                {
                    static unsigned int counter = 0;
                    if (++counter>100)
                    {
                        counter=0;
                        if (getKineticEnergy()<1e-16) // *getElasticEnergy()
                            return false;
                    }
                    return true;
                }*/
            //
            /*    void actionsAfterTimeStep() override
                {
                    // Remove top wall
                    static bool wallRemoved = false;
                    if (wallRemoved == false && getTime() == getTimeMax())
                    {
                        //logger(INFO,"walls removed");
                        wallHandler.removeObject(6);
                        wallRemoved = true;
                    }
                }*/
            //
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
            // No need to export these info here:
/*            void writeFstatHeader(std::ostream& os) const override
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

            }*/
            //
            void printTime() const override
            {
                logger(INFO,"Time % TimeMax % kineticEnergy %", getTime(), getTimeMax(),getKineticEnergy());
            }
        };
//
int main(int argc UNUSED, char *argv[] UNUSED)
{
    InitialBed problem;
    problem.solve();
    return 0;
}
