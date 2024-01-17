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


class BedFixAdd: public Mercury3D
{
private:
    // set main simulation variables
    Mdouble timeStep = 1e-7;
    Mdouble scaleMass = 1;
    unsigned int laserON = 0;
    bool writeVTK = true;
    bool fixBaseParticles = false;
    bool deleteParticles = true;
    Mdouble MaxSimTime = 0.05;
    //
    int spatialScale = 2;//4
    int volumeScale = 3;//6
    //
    Mdouble gravityValue = -9.81;
    Mdouble XMaxDomain = spatialScale*500e-6/2;
    Mdouble YMaxDomain = spatialScale*500e-6/2;
    Mdouble ZMinDomain = 400e-6;//200e-6;
    Mdouble ZMaxDomain = 100e-6;
    //
    Mdouble MatDensity = 1050;//1400.0;
    Mdouble elasticModulus = 8e7;//8e7;//1e7;
    Mdouble possionsRation = 0.4;//0.2;
    Mdouble damping = 0.9;//5e-6*10.0;
    //
    Mdouble latentHeat = 56.4e3;//245e3;
    Mdouble heatCapacity = 1200;//2270;
    //
    Mdouble surfaceTension = 34.3e-3;//0.035;
    Mdouble refVis = 2e-3;//2e-3;
    // cond , conv, rad
    Mdouble thermalcond = 0.12;
    Mdouble thermalConvCoff = 150.0;
    Mdouble emmisivity = 0.9;
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
    int thermalExpansionOnOff = 0;
    //
    Mdouble particlesTemp = 155 + 273.15;//413.15;
    //
    Mdouble heatingTime = 0.0;//0.0075;
    Mdouble heatSource = 0.0;//2.5e-4;//deposition_problem.heat;//2.5e-4;
    //
    unsigned int numberOfCounts = 1e4;
    //
    // Domain
    Mdouble ToolLength = spatialScale*500e-6/2;//0.00025; //45e-1 * scale * 2.0; // unit cm
    Mdouble ToolWidth = spatialScale*500e-6/2;//0.00025;//30e-1 * scale * 2.0;
    Mdouble ToolHight = 0.01;//40e-1 * scale * 4.0; //20e-1 * scale * 4.0;
    Mdouble ToolThickness = 50e-6;
    Mdouble Gap = 30e-6;
    //Front Dam
    Mdouble damThickness = ToolLength + ToolThickness;//(5.0e-1 * scale);
    //
    //Front Dam
    Vec3D minPoint_dam = Vec3D(ToolLength, 0.0, -ZMinDomain);
    Vec3D maxPoint_dam = Vec3D(damThickness, ToolWidth, ToolHight);
    //
    //
    Mdouble meanRadius1 = 27e-6;// in m
    Mdouble MinRadius = 20e-6;//22.6e-6/2.0;  //19.5e-6/2.0;//in m
    Mdouble MaxRadius = 75e-6;//153.4/2.0;    //174.5e-6/2.0;//in m
    Mdouble LogmeanRadius = -10.3693;
    Mdouble LogstdRadius = 0.2758;
    //
    //the volume to be added:
    Mdouble addVolume = volumeScale*Gap*ToolLength*ToolWidth;//(20e-1 * scale)*ToolLength*ToolWidth;//(ToolHight/4.0)*ToolLength*ToolWidth; //0.0018 cm3
    //
    MeltableSpecies *species;
    //
public:
    //
    BedFixAdd() {
        setName("InitialBed");
        //
        readRestartFile();
        setRestarted(true);
        //
        if (fixBaseParticles)
        {
            for (int i = 0; i < particleHandler.getNumberOfObjects(); ++i) {
                Vec3D pP = particleHandler.getObject(i)->getPosition();
                if (pP.Z <= -160.0 * 1e-6)
                    particleHandler.getObject(i)->fixParticle();
            }
            logger(INFO, " Fixed % particles", particleHandler.getNumberOfObjects());
        }
        //
/*        if(deleteParticles)
        {
            for (int i = 0; i < particleHandler.getNumberOfObjects(); ++i) {
                Vec3D pP = particleHandler.getObject(i)->getPosition();
                if (pP.Z >= 0.0)
                    particleHandler.removeObject(i);
            }
            logger(INFO, " Removed % particles", particleHandler.getNumberOfObjects());
        }*/
        //
        setName("BedFixAdd");
        setSystemDimensions(3);
        setGravity(Vec3D(0.0,0.0,gravityValue));
        //setXMin(0.0);
        //setYMin(0.0);
        setZMin(-ZMinDomain);
        setXMax(XMaxDomain);
        setYMax(YMaxDomain);
        setZMax(ZMaxDomain);
        setTimeMax(MaxSimTime);
        species = dynamic_cast<MeltableSpecies *>(speciesHandler.getObject(0));
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
        //
        species->setHeatingTime(heatingTime);
        species->setHeatInput(heatSource);
        //
        //
        setSaveCount(numberOfCounts);
        setParticlesWriteVTK(writeVTK);
        setWallsWriteVTK(writeVTK);
        //
        //species->setMaxTimeStep(MinRadius,MatDensity,elasticModulus);
        setTimeStep(timeStep);
        //
        std::mt19937 gen;
        gen.seed(10);
        std::lognormal_distribution<> d(LogmeanRadius, LogstdRadius);
        ThermalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(meanRadius1);
        p.setTemperature(particlesTemp);
        Mdouble fillHeight = 0.0;
        while (addVolume > 0) {
            Mdouble x = random.getRandomNumber(ToolThickness, ToolLength);
            Mdouble y = random.getRandomNumber(ToolThickness, (ToolWidth - ToolThickness));
            Mdouble z = random.getRandomNumber(-100e-6, fillHeight);
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
                fillHeight += 30e-6 *meanRadius1;
            }
        }
        logger(INFO, " Inserted % particles", particleHandler.getNumberOfObjects());
        //
    };
    //
    //void setupInitialConditions() override {}
    //
};
//
int main(int argc UNUSED, char *argv[] UNUSED)
{
    BedFixAdd problem;
    problem.solve();
    return 0;
}
