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
#include "Species/MeltableFrictionSpecies.h"
//#include "Species/MeltableSpecies.h"
//#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Particles/ThermalParticle.h"
#include "Boundaries/PeriodicBoundary.h"
//#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include <random>
//#include <DPMBase.h>


class SinterInitialBed: public Mercury3D
{
public:
    //
    Mdouble bottomWallZ = -3e-4;
    Mdouble maxSimTime = 0.03;//6.0;
    unsigned int numberOfCounts = 10000;
    //
    // material properties and simulation parameters:
    Mdouble particleInitialTemp = 155 + 273.15;//20 + 273.15;//428.15
    //
    Mdouble MatDensity = 1050;//1400.0;
    Mdouble elasticModulus = 8e7;//8e7;//1e7;
    Mdouble possionsRation = 0.4;//0.2;
    Mdouble damping = 0.95;//0.9;//5e-6*10.0;
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
    int thermalExpansionOnOff = 1;
    //------------------------------------------------------------------------------------------------
    //particles:
    Mdouble minRadius = 40e-6;
    Mdouble maxRadius = 50e-6;
    Mdouble meanRadius1 = 0.5*(minRadius+maxRadius);
    //------------------------------------------------------------------------------------------------
    //
    //Mdouble scale = 5e-3; //1e-1; //copy to speed and simulation boundaries
    //
    void setupInitialConditions() override {

        //TOOL:
        Mdouble ToolLength = 25e-5;//0.00025; //45e-1 * scale * 2.0; // unit cm
        Mdouble ToolWidth = 25e-5;//0.00025;//30e-1 * scale * 2.0;
        Mdouble ToolHight = 0.01;//40e-1 * scale * 4.0; //20e-1 * scale * 4.0;
        Mdouble ToolThickness = 25e-6;//0.000025;//5e-1 * scale;
        Mdouble Gap = 0.1e-3;
        //Front Dam
        Mdouble damThickness = ToolLength + ToolThickness;//(5.0e-1 * scale);
        //TOOL
        //side wall
        Vec3D minPoint1_tool = Vec3D(0.0, 0.0, 0.0);
        Vec3D maxPoint1_tool = Vec3D(ToolLength, ToolThickness, ToolHight); //Vec3D(getXMax(),0.05*getYMax(),getZMax())
        //back wall
        Vec3D minPoint2_tool = Vec3D(0.0, 0.0, 0.0); //Gap);
        Vec3D maxPoint2_tool = Vec3D(ToolThickness, ToolWidth, ToolHight); //Vec3D(ToolThickness,(ToolWidth-ToolThickness),ToolHight);
        //side wall
        Vec3D minPoint3_tool = Vec3D(0.0, (ToolWidth - ToolThickness), 0.0);
        Vec3D maxPoint3_tool = Vec3D(ToolLength, ToolWidth, ToolHight); // Vec3D(getXMax(),getYMax(),getZMax())
        //
        //Front Dam
        Vec3D minPoint_dam = Vec3D(ToolLength, 0.0, 0.0);
        Vec3D maxPoint_dam = Vec3D(damThickness, ToolWidth, ToolHight);
        //Back Dam
        Vec3D minPoint_back_dam = Vec3D(0.0, 0.0, 0.0); //Vec3D(0.0,ToolThickness,0.0);
        Vec3D maxPoint_back_dam = Vec3D(ToolThickness, ToolWidth, Gap); //Vec3D(ToolThickness,(ToolWidth-ToolThickness),Gap);
        //
        //Bottom Wall
        //typedef ThermalInfiniteWall<InfiniteWall> thermalWall;
        //thermalWall w0;
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0)); //w0.setSpecies(speciesHandler.getObject(1));
        //w0.setSolidHeatCapacity(inf);
        //w0.setThermalConductivity(0.0);
        //w0.setEmmisivity(0.0);
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, bottomWallZ));
        wallHandler.copyAndAddObject(w0);
        //
/*        IntersectionOfWalls w1;
        w1.setSpecies(speciesHandler.getObject(0));
        w1.addObject(Vec3D(1.0,0.0,0.0),minPoint1_tool);
        w1.addObject(Vec3D(0.0,1.0,0.0),minPoint1_tool);
        w1.addObject(Vec3D(0.0,0.0,1.0),minPoint1_tool);
        w1.addObject(Vec3D(-1.0,0.0,0.0),maxPoint1_tool);
        w1.addObject(Vec3D(0.0,-1.0,0.0),maxPoint1_tool);
        w1.addObject(Vec3D(0.0,0.0,-1.0),maxPoint1_tool);
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
        w3.addObject(Vec3D(1.0,0.0,0.0),minPoint3_tool);
        w3.addObject(Vec3D(0.0,1.0,0.0),minPoint3_tool);
        w3.addObject(Vec3D(0.0,0.0,1.0),minPoint3_tool);
        w3.addObject(Vec3D(-1.0,0.0,0.0),maxPoint3_tool);
        w3.addObject(Vec3D(0.0,-1.0,0.0),maxPoint3_tool);
        w3.addObject(Vec3D(0.0,0.0,-1.0),maxPoint3_tool);
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
        */
        //
        //Back Dam
/*
        IntersectionOfWalls w3;
        w3.setSpecies(speciesHandler.getObject(0));
        w3.addObject(Vec3D(1.0, 0.0, 0.0), minPoint_back_dam);
        w3.addObject(Vec3D(0.0, 1.0, 0.0), minPoint_back_dam);
        w3.addObject(Vec3D(0.0, 0.0, 1.0), minPoint_back_dam);
        w3.addObject(Vec3D(-1.0, 0.0, 0.0), maxPoint_back_dam);
        w3.addObject(Vec3D(0.0, -1.0, 0.0), maxPoint_back_dam);
        w3.addObject(Vec3D(0.0, 0.0, -1.0), maxPoint_back_dam);
        wallHandler.copyAndAddObject(w3);
*/
        //------------------------------------------------------------------------------------------------
        PeriodicBoundary b0;
        b0.set(Vec3D(0, 1, 0), 0.0, ToolWidth);
        boundaryHandler.copyAndAddObject(b0);
        PeriodicBoundary b1;
        b1.set(Vec3D(1, 0, 0), 0.0, ToolLength);
        boundaryHandler.copyAndAddObject(b1);
        //------------------------------------------------------------------------------------------------
        //
        //the volume to be added:
        Mdouble addVolume = Gap*ToolLength*ToolWidth;
        //(20e-1 * scale)*ToolLength*ToolWidth;//(ToolHight/4.0)*ToolLength*ToolWidth; //0.0018 cm3
        //------------------------------------------------------------------------------------------------
        // //set up random number generator
        //std::random_device rd;
        //std::mt19937 gen(rd());
        //
        std::mt19937 gen;
        gen.seed(1);
        std::uniform_real_distribution<> d(minRadius, maxRadius);
        //
        //add particles until the volume to be added is zero
        //logger(INFO,"Adding particles ...");
        ThermalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(meanRadius1);
        p.setTemperature(particleInitialTemp);
        Mdouble fillHeight = 0.0;
        while (addVolume > 0) {
            Mdouble x = random.getRandomNumber(ToolThickness, ToolLength);
            Mdouble y = random.getRandomNumber(ToolThickness, (ToolWidth-ToolThickness));
            Mdouble z = random.getRandomNumber(bottomWallZ, fillHeight);
            p.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p)) {
                particleHandler.copyAndAddObject(p);
                addVolume -= p.getVolume();
                do {
                    p.setRadius(d(gen));
                } while (p.getRadius() < minRadius || p.getRadius() > maxRadius);
                if (particleHandler.getNumberOfObjects() % 100 == 0) std::cout << '.' << std::flush;
                if (particleHandler.getNumberOfObjects() % 1000 == 0) std::cout << ' ';
                if (particleHandler.getNumberOfObjects() % 10000 == 0) std::cout << addVolume << '\n';
            } else {
                fillHeight += 0.01*meanRadius1;//0.01 * scale * meanRadius1; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        logger(INFO, " Inserted % particles", particleHandler.getNumberOfObjects());
        //------------------------------------------------------------------------------------------------
    }
};
//
int main(int argc UNUSED, char *argv[] UNUSED)
{
    //logger(INFO,"Simple box for creating particles");
    SinterInitialBed problem;
    //
    Mdouble gravityValue = -9.81;
    Mdouble XMaxDomain = 25e-5;//0.00025;//problem.scale*45e-1+(2.0*problem.scale);
    Mdouble YMaxDomain = 25e-5;//0.00025;//problem.scale*30e-1+(2.0*problem.scale);
    Mdouble ZMaxDomain = 0.01;//problem.scale*40e-1*(4.0);//*problem.scale);
    //
    //
    //-----------------------------------------------------------------------------------------------------------------------------------------------------
/*    Mdouble SlidingFrictionCoeff = 0.25;//0.25;
    Mdouble RollingFrictionCoeff = 0.05;//0.05;*/
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    //
    problem.setName("SinteringInitialBedUniform");
    problem.setSystemDimensions(3);
    problem.setGravity(Vec3D(0.0,0.0,gravityValue));
    //deposition_problem.setXMin(0.0);
    //deposition_problem.setYMin(0.0);
    problem.setZMin(problem.bottomWallZ);
    problem.setXMax(XMaxDomain);
    problem.setYMax(YMaxDomain);
    problem.setZMax(ZMaxDomain);
    problem.setTimeMax(problem.maxSimTime);
    //problem.setHGridMaxLevels(2);
    //
    auto species = problem.speciesHandler.copyAndAddObject(MeltableFrictionSpecies());
    //auto species2 = deposition_problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    //auto wallSpecies = deposition_problem.speciesHandler.copyAndAddObject(species2);
    //auto wallParticleSpecies = deposition_problem.speciesHandler.getMixedObject(species2,wallSpecies);
    //
    species->setDensity(problem.MatDensity);
    species->setElasticModulus(problem.elasticModulus);
    species->setPoissonRatio(problem.possionsRation);
    species->setDissipation(problem.damping);
    //
    species->setLatentHeat(problem.latentHeat);
    species->setSolidHeatCapacity(problem.heatCapacity);
    species->setThermalConductivityCoefficient(problem.thermalcond);
    species->setThermalConvectionCoefficient(problem.thermalConvCoff);
    species->setAmbientTemperature(problem.aTemp);
    species->setMeltingTemperature(problem.mTemp);
    species->setMaterialEmissivity(problem.emmisivity);
    //
    species->setSurfaceTension(problem.surfaceTension);
    species->setRefViscosity(problem.refVis);
    //
    species->setDeltaT(problem.deltaT);
    species->setLiquidHeatCapacity(problem.liquidHeatCapacity);
    species->setVaporizationHeatCapacity(problem.vaporizationHeatCapacity);
    species->setVaporizationLatentHeat(problem.vaporizationLatentHeat);
    species->setVaporizationTemperature(problem.vaporizationTemp);
    //
    species->setThermalExpansionCoefficient(problem.thermalExpansionCoeff);
    species->setThermalExpansionOnOff(problem.thermalExpansionOnOff);
    //
/*    species->setHeatingTime(problem.heatingTime);
    species->setHeatInput(problem.heatSource);
    species->setCoolingTime(problem.coolingTime);*/
    //
    species->setMaxTimeStep(problem.minRadius,problem.MatDensity,problem.elasticModulus);
    //------------------------------------------------------------------------------Friction coeffiecients Values---------------------------------------------------------------------------
/*    species->setSlidingStiffness(2./7.*species->getEffectiveElasticModulus());
    species->setSlidingDissipation(2./7.*species->getViscosityCoefficient());
    species->setSlidingFrictionCoefficient(SlidingFrictionCoeff);
    //species->setSlidingFrictionCoefficientStatic(0.6);
    //
    species->setRollingStiffness(2./5.*species->getElasticModulus());
    species->setRollingDissipation(2./5.*species->getDissipation());
    species->setRollingFrictionCoefficient(RollingFrictionCoeff);*/
    //------------------------------------------------------------------------------Data Output--------------------------------------------------------------------------------------------
    problem.setSaveCount(problem.numberOfCounts);//50
    //problem.restartFile.setSaveCount(10000);
    //deposition_problem.dataFile.setFileType(FileType::NO_FILE);javascript:void(0)
    //deposition_problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.setParticlesWriteVTK(true);
    problem.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    //
    problem.setTimeStep(species->getMaxTimeStep());
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    problem.solve();
    return 0;
}
