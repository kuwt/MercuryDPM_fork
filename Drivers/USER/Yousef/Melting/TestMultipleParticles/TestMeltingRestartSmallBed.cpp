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
#include "../Initialization.h"
#include <random>



class TestMeltingRestartSmallBed: public Initialization
{
private:
    //
    Mdouble timeStep = 1e-7;
    unsigned int addParticles = 0;
    unsigned int setNewSpecies = 0; // 0->off
    unsigned int logNormalDistParticles = 1;
    bool writeVTK = false;
    //
    Mdouble monoRadius = 15e-6;
    //
    Mdouble meanRadius1 = 27e-6;// in m
    Mdouble MinRadius = 20e-6;//22.6e-6/2.0;  //19.5e-6/2.0;//in m
    Mdouble MaxRadius = 75e-6; //monoRadius;//75e-6;    //153.4/2.0;    //174.5e-6/2.0;//in m
    Mdouble LogmeanRadius = -10.3693;
    Mdouble LogstdRadius = 0.2758;
    //
    Mdouble addVolume = MaxRadius * MaxRadius * MaxRadius*2;
    //
    MeltableSpecies *previousSpecies;
    //
public:
    TestMeltingRestartSmallBed()
    {
        setName("TestMeltingSmallParticlesBed");
        //setName("TestMeltingRestartSmallBed");
        readRestartFile();
        //setRestarted(false);
        setName("TestMeltingRestartSmallBed");
        //setName("TestMeltingRestartSmallBed2");
        previousSpecies = dynamic_cast<MeltableSpecies *>(speciesHandler.getObject(0));
        //
        setSystemDimensions(3);
        setGravity(Vec3D(0.0, 0.0, gravityAcc));
        setXMin(0.0);
        setYMin(0.0);//0.0
        setZMin(0.0);
        setXMax(MaxRadius*2); //1
        setYMax(MaxRadius*2); //1
        setZMax(MaxRadius*8); //3
        setTimeMax(maxSimTime);
        // change properties of previous specious
        // e.g. previousSpecies->setHeatInput(1e-3);
        //
        switch (setNewSpecies) {
            case 0:
            {
                previousSpecies->setDensity(MatDensity);
                previousSpecies->setElasticModulus(elasticModulus);
                previousSpecies->setPoissonRatio(possionsRation);
                previousSpecies->setDissipation(damping);
                previousSpecies->setLatentHeat(latentHeat);
                previousSpecies->setSolidHeatCapacity(heatCapacity);
                previousSpecies->setThermalConductivityCoefficient(thermalcond);
                previousSpecies->setThermalConvectionCoefficient(thermalConvCoff);
                previousSpecies->setAmbientTemperature(aTemp);
                previousSpecies->setMeltingTemperature(mTemp);
                previousSpecies->setMaterialEmissivity(emmisivity);
                previousSpecies->setHeatingTime(heatingTime);
                previousSpecies->setHeatInput(heatSource);
                previousSpecies->setCoolingTime(coolingTime);
                previousSpecies->setSurfaceTension(surfaceTension);
                previousSpecies->setRefViscosity(refVis);
                previousSpecies->setDeltaT(deltaT);
                previousSpecies->setLiquidHeatCapacity(liquidHeatCapacity);
                previousSpecies->setVaporizationHeatCapacity(vaporizationHeatCapacity);
                previousSpecies->setVaporizationLatentHeat(vaporizationLatentHeat);
                previousSpecies->setVaporizationTemperature(vaporizationTemp);
                previousSpecies->setThermalExpansionCoefficient(thermalExpansionCoeff);
                previousSpecies->setThermalExpansionOnOff(thermalExpansionOnOff);
            }break;
            case 1:
            {
                auto species = speciesHandler.copyAndAddObject(MeltableSpecies());
                species->setDensity(MatDensity);
                species->setElasticModulus(elasticModulus);
                species->setPoissonRatio(possionsRation);
                species->setDissipation(damping);
                species->setLatentHeat(latentHeat);
                species->setSolidHeatCapacity(heatCapacity);
                species->setThermalConductivityCoefficient(thermalcond);
                species->setThermalConvectionCoefficient(thermalConvCoff);
                species->setAmbientTemperature(aTemp);
                species->setMeltingTemperature(mTemp);
                species->setMaterialEmissivity(emmisivity);
                species->setHeatingTime(heatingTime);
                species->setHeatInput(heatSource);
                species->setCoolingTime(coolingTime);
                species->setSurfaceTension(surfaceTension);
                species->setRefViscosity(refVis);
                species->setDeltaT(deltaT);
                species->setLiquidHeatCapacity(liquidHeatCapacity);
                species->setVaporizationHeatCapacity(vaporizationHeatCapacity);
                species->setVaporizationLatentHeat(vaporizationLatentHeat);
                species->setVaporizationTemperature(vaporizationTemp);
                species->setThermalExpansionCoefficient(thermalExpansionCoeff);
                species->setThermalExpansionOnOff(thermalExpansionOnOff);
            }
        }

        //
        setTimeStep(timeStep); //species->getMaxTimeStep());//1e-8
        //
        //
        setSaveCount(saveCount);
        restartFile.setFileType(FileType::ONE_FILE);
        //restartFile.setSaveCount(restartSaveCount);
        //restartFile.writeFirstAndLastTimeStep();
        dataFile.setFileType(FileType::ONE_FILE);
        fStatFile.setFileType(FileType::ONE_FILE);
        //fStatFile.setSaveCount(10000);
        eneFile.setFileType(FileType::NO_FILE);
        //
        setParticlesWriteVTK(writeVTK);
        setWallsWriteVTK(writeVTK);
        //
        //
        switch (addParticles) {
            case 1: {
                switch (logNormalDistParticles) {
                    case 1: {
                        //set up random number generator
                        //std::random_device rd;
                        //std::mt19937 gen(rd());
                        std::mt19937 gen;
                        gen.seed(3786497);
                        //
                        std::lognormal_distribution<> d(LogmeanRadius, LogstdRadius);
                        //add particles until the volume to be added is zero
                        //logger(INFO,"Adding particles ...");
                        ThermalParticle p;
                        p.setSpecies(speciesHandler.getObject(setNewSpecies));
                        p.setRadius(meanRadius1);
                        p.setTemperature(particleInitialTemp);
                        Mdouble fillHeight = 0.0;
                        while (addVolume > 0) {
                            Mdouble x = random.getRandomNumber(getXMin(), getXMax());
                            Mdouble y = random.getRandomNumber(getYMin(), getYMax());
                            Mdouble z = random.getRandomNumber(getZMin(), fillHeight);
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
                                fillHeight += 30e-6 * meanRadius1;
                            }
                        }
                        logger(INFO, " Inserted % particles", particleHandler.getNumberOfObjects());
                    }
                        break;
                        //
                    case 0: {
                        ThermalParticle p;
                        p.setSpecies(speciesHandler.getObject(setNewSpecies));
                        p.setRadius(monoRadius);
                        p.setTemperature(particleInitialTemp);
                        Mdouble fillHeight = 0.0;
                        while (addVolume > 0) {
                            Mdouble x = random.getRandomNumber(getXMin(), getXMax());
                            Mdouble y = random.getRandomNumber(getYMin(), getYMax());
                            Mdouble z = random.getRandomNumber(getZMin(), fillHeight);
                            p.setPosition({x, y, z});
                            // check if particle can be inserted
                            if (checkParticleForInteraction(p)) {
                                particleHandler.copyAndAddObject(p);
                                addVolume -= p.getVolume();
                                //
                                if (particleHandler.getNumberOfObjects() % 100 == 0) std::cout << '.' << std::flush;
                                if (particleHandler.getNumberOfObjects() % 1000 == 0) std::cout << ' ';
                                if (particleHandler.getNumberOfObjects() % 10000 == 0) std::cout << addVolume << '\n';
                            } else {
                                fillHeight += 30e-6 * meanRadius1;
                            }
                        }
                        logger(INFO, " Inserted % particles", particleHandler.getNumberOfObjects());
                    }
                }
            }break;
            case 0: {
                logger(INFO, " NO Inserted particles, total Particles %", particleHandler.getNumberOfObjects());
            }
        }
        //
    };
    //
/*    void setupInitialConditions() override {
    }*/
//
// printing to screen
    void printTime() const override
    {
        logger(INFO,"Time % TimeMax % timeStep %", getTime(), getTimeMax(), getTimeStep());
        //
/*        for (const auto& q: particleHandler) {
            auto p0 = dynamic_cast<ThermalParticle *>(q);
            logger.assert_debug(p0 != nullptr, "Thermal Particles required");
            logger(INFO, "time % temperature % moltenLayer %", getTime(), p0->getTemperature(), p0->getMoltenLayerThickness());
        }
        //
        for (const auto &iH: interactionHandler) {
            auto i0 = dynamic_cast<MeltableInteraction *>(iH);
            logger.assert_debug(i0 != nullptr, "i0 not set");
            //double force = i0->getForce().y()/i0->getNormal().y();
            logger(INFO, "time % force % overlap % solidOverlap % bondingOverlap % contactRadius %", getTime(), i0->getForce(), i0->getOverlap(),
                   i0->getSolidOverlap(), i0->getBondingOverlap(), i0->getContactRadiusMeltable());
        }*/
    }
};


int main(int argc UNUSED, char* argv[] UNUSED)
{
    TestMeltingRestartSmallBed problem;
    problem.setXBallsAdditionalArguments("-solidf -v0 -s .85");
    problem.solve();
    return 0;
}
