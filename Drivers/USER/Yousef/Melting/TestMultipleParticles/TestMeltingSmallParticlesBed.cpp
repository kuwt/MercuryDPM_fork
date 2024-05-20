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
#include "../Initialization.h"
#include <random>
#include "Walls/IntersectionOfWalls.h"



class TestMeltingSmallParticlesBed: public Initialization
{
private:
    //
    Mdouble timeStep = 1e-9;
    //
    unsigned int wallsBoundary = 1; // if false decrease domain size, to allow particles to be top of each other
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
public:
    void setupInitialConditions() override {
        //
        setName("TestMeltingSmallParticlesBed");
        setSystemDimensions(3);
        setGravity(Vec3D(0.0, 0.0, gravityAcc));
        setXMin(0.0);
        setYMin(0.0);//0.0
        setZMin(0.0);
        setXMax(MaxRadius*2); //1
        setYMax(MaxRadius*2); //1
        setZMax(MaxRadius*8); //3
        setTimeMax(maxSimTime);
        //setHGridMaxLevels(3);
        //
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
        //
        //species->setMaxTimeStep(0.4*MinRadius,MatDensity,elasticModulus);
        setTimeStep(timeStep);
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
        // Sim setup:
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set({0.0, 0.0, -1.0}, {0.0, 0.0, 0.0});
        wallHandler.copyAndAddObject(w0);
        //
        Mdouble thicknessW = minRadius;
        //
        Vec3D minP1 = Vec3D(getXMin(),getYMin(),getZMin());
        Vec3D maxP1 = Vec3D(getXMin()+thicknessW,getYMax(),getZMax());
        //
        Vec3D minP2 = minP1;//Vec3D(getXMin(),getYMin(),getZMin());
        Vec3D maxP2 = Vec3D(getXMax(),getYMin()+thicknessW,getZMax());
        //
        Vec3D minP3 = Vec3D(getXMax()-thicknessW,getYMin(),getZMin());
        Vec3D maxP3 = Vec3D(getXMax(),getYMax(),getZMax());
        //
        Vec3D minP4 = Vec3D(getXMin(),getYMax()-thicknessW,getZMin());
        Vec3D maxP4 = Vec3D(getXMax(),getYMax(),getZMax());
        //
        switch (wallsBoundary)
        {
            case 1:
            {
                IntersectionOfWalls w1;
                w1.setSpecies(speciesHandler.getObject(0));
                w1.addObject(Vec3D(1.0, 0.0, 0.0), minP1);
                w1.addObject(Vec3D(0.0, 1.0, 0.0), minP1);
                w1.addObject(Vec3D(0.0, 0.0, 1.0), minP1);
                w1.addObject(Vec3D(-1.0, 0.0, 0.0), maxP1);
                w1.addObject(Vec3D(0.0, -1.0, 0.0), maxP1);
                w1.addObject(Vec3D(0.0, 0.0, -1.0), maxP1);
                wallHandler.copyAndAddObject(w1);
                //
                IntersectionOfWalls w2;
                w2.setSpecies(speciesHandler.getObject(0));
                w2.addObject(Vec3D(1.0, 0.0, 0.0), minP2);
                w2.addObject(Vec3D(0.0, 1.0, 0.0), minP2);
                w2.addObject(Vec3D(0.0, 0.0, 1.0), minP2);
                w2.addObject(Vec3D(-1.0, 0.0, 0.0), maxP2);
                w2.addObject(Vec3D(0.0, -1.0, 0.0), maxP2);
                w2.addObject(Vec3D(0.0, 0.0, -1.0), maxP2);
                wallHandler.copyAndAddObject(w2);
                //
                IntersectionOfWalls w3;
                w3.setSpecies(speciesHandler.getObject(0));
                w3.addObject(Vec3D(1.0, 0.0, 0.0), minP3);
                w3.addObject(Vec3D(0.0, 1.0, 0.0), minP3);
                w3.addObject(Vec3D(0.0, 0.0, 1.0), minP3);
                w3.addObject(Vec3D(-1.0, 0.0, 0.0), maxP3);
                w3.addObject(Vec3D(0.0, -1.0, 0.0), maxP3);
                w3.addObject(Vec3D(0.0, 0.0, -1.0), maxP3);
                wallHandler.copyAndAddObject(w3);
                //
                IntersectionOfWalls w4;
                w4.setSpecies(speciesHandler.getObject(0));
                w4.addObject(Vec3D(1.0, 0.0, 0.0), minP4);
                w4.addObject(Vec3D(0.0, 1.0, 0.0), minP4);
                w4.addObject(Vec3D(0.0, 0.0, 1.0), minP4);
                w4.addObject(Vec3D(-1.0, 0.0, 0.0), maxP4);
                w4.addObject(Vec3D(0.0, -1.0, 0.0), maxP4);
                w4.addObject(Vec3D(0.0, 0.0, -1.0), maxP4);
                wallHandler.copyAndAddObject(w4);
        }break;
            //
            case 0:
            {
                PeriodicBoundary by1;
                by1.set(Vec3D(0, 1, 0), 0.0, getYMax());
                boundaryHandler.copyAndAddObject(by1);
                PeriodicBoundary bx1;
                bx1.set(Vec3D(1, 0, 0), 0.0, getXMax());
                boundaryHandler.copyAndAddObject(bx1);
            }
        }
        //
        Mdouble addVolume = MaxRadius * MaxRadius * MaxRadius*2;
        //
        switch (logNormalDistParticles)
        {
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
                p.setSpecies(speciesHandler.getObject(0));
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
            }break;
            //
            case 0:
            {
                ThermalParticle p;
                p.setSpecies(speciesHandler.getObject(0));
                p.setRadius(monoRadius);
                p.setTemperature(particleInitialTemp);
                Mdouble fillHeight = 0.0;
                while (addVolume > 0) {
                    Mdouble x = random.getRandomNumber(getXMin(), getXMax());
                    Mdouble y = random.getRandomNumber(getYMin(), getYMax());
                    Mdouble z = random.getRandomNumber(getZMin(), fillHeight);
                    p.setPosition({x, y, z});
                    // check if particle can be inserted
                    if (checkParticleForInteraction(p))
                    {
                        particleHandler.copyAndAddObject(p);
                        addVolume -= p.getVolume();
                        //
                        if (particleHandler.getNumberOfObjects() % 100 == 0) std::cout << '.' << std::flush;
                        if (particleHandler.getNumberOfObjects() % 1000 == 0) std::cout << ' ';
                        if (particleHandler.getNumberOfObjects() % 10000 == 0) std::cout << addVolume << '\n';
                    }
                    else
                    {
                        fillHeight += 30e-6 *meanRadius1;
                    }
                }
                logger(INFO, " Inserted % particles", particleHandler.getNumberOfObjects());
            }
        }
        //
    }

/*    void writeFstatHeader(std::ostream& os) const override
    {
        for (const auto &q: particleHandler) {
            auto p0 = dynamic_cast<ThermalParticle *>(q);
            logger.assert_debug(p0 != nullptr, "Thermal Particles required");
            os << p0->getTemperature() << " " << p0->getMoltenLayerThickness()
               << " " << 00 << " " << 00 << " " << 00 <<  " " << 00 << " " << 00
               << " " << 00 << " " << 00 << " " << 00 << " " << 00 << " " << 00 << std::endl;
        }
        //
        for (const auto &iH: interactionHandler) {
            auto i0 = dynamic_cast<MeltableInteraction *>(iH);
            logger.assert_debug(i0 != nullptr, "i0 not set");
            //double force = i0->getForce().z()/i0->getNormal().z();
            //Mdouble force = i0->getAbsoluteNormalForce();
            os << i0->getForce() << " " << i0->getOverlap()
               << " " << i0->getSolidOverlap() << " " << i0->getBondingOverlap() << " " << i0->getContactRadiusMeltable()
               << " " << i0->getElasticForce() << " " << i0->getDampingForce()
               << " " << i0->getBondingForce() << " " << i0->getViscousForce() << " " << i0->getTensionForce() << std::endl;
        }
    }*/
//
// printing to screen
    void printTime() const override
    {
        logger(INFO,"Time % TimeMax %", getTime(), getTimeMax());
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
    TestMeltingSmallParticlesBed problem;
    problem.setXBallsAdditionalArguments("-solidf -v0 -s .85");
    problem.solve();
    return 0;
}
