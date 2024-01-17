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
#include "Species/MeltableSpecies.h"
//#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Particles/ThermalParticle.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include <random>
//#include <DPMBase.h>


class SinterInitialBed: public Mercury3D
{
public:
    //
    //Mdouble scale = 5e-3; //1e-1; //copy to speed and simulation boundaries
    //
    void setupInitialConditions() {

        //TOOL:
        Mdouble ToolLength = 0.00025; //45e-1 * scale * 2.0; // unit cm
        Mdouble ToolWidth = 0.00025;//30e-1 * scale * 2.0;
        Mdouble ToolHight = 0.01;//40e-1 * scale * 4.0; //20e-1 * scale * 4.0;
        Mdouble ToolThickness = 0.000025;//5e-1 * scale;
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
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, 0.0));
        wallHandler.copyAndAddObject(w0);
        //
        IntersectionOfWalls w1;
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

        //
/*
        PeriodicBoundary b0;
        b0.set(Vec3D(0, 1, 0), 0.0, ToolWidth);
        boundaryHandler.copyAndAddObject(b0);
        PeriodicBoundary b1;
        b1.set(Vec3D(1, 0, 0), 0.0, ToolLength);
        boundaryHandler.copyAndAddObject(b1);*/
        //
        //particles:
        Mdouble meanRadius1 = 0.000019; // V2_1 -> 18.8075e-4; //17.3370e-4; //17.572e-4;  //35.14e-4/2;    //17.59e-4;//35.18e-4/2;
        //Mdouble mode_log_normal = 34.674e-4 / 2.0;
        //Mdouble stdRadius1 = 1.4027; // 1.7824; //1.8549; //0.2 * meanRadius1;
        // log-normal distribution mean and std, calculated from data points via lognfit matlab function
        Mdouble LogmeanRadius = -10.8936 ;//-6.2884; //-6.4266;//-6.3575;//log(meanRadius1);
        Mdouble LogstdRadius = 0.3384;//0.3384; //0.6178; //log(stdRadius1);
        Mdouble MinRadius = 23.0e-6/2.0;//(13.183e-4 / 2.0); //10e-4; //diameter = 20//10microM = 10e-3 mm = 10e-4 cm
        Mdouble MaxRadius = 60.0e-6/2.0;//(79.433e-4 / 2.0); //27.5e-4; // diameter = 55 // 27.5microM = 27.5e-3 mm = 27.5e-4cm
        //
        //the volume to be added:
        Mdouble addVolume = Gap*ToolLength*ToolWidth;//(20e-1 * scale)*ToolLength*ToolWidth;//(ToolHight/4.0)*ToolLength*ToolWidth; //0.0018 cm3
        // //set up random number generator
        //std::random_device rd;
        //std::mt19937 gen(rd());
        //
        std::mt19937 gen;
        gen.seed(1);
        std::lognormal_distribution<> d(LogmeanRadius, LogstdRadius);
        //
        //add particles until the volume to be added is zero
        //logger(INFO,"Adding particles ...");
        ThermalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(meanRadius1);
        p.setTemperature(443.15);
        Mdouble fillHeight = 0.0;
        while (addVolume > 0) {
            Mdouble x = random.getRandomNumber(ToolThickness, ToolLength);
            Mdouble y = random.getRandomNumber(ToolThickness, (ToolWidth-ToolThickness));
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
                fillHeight += 0.01*meanRadius1;//0.01 * scale * meanRadius1; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        logger(INFO, " Inserted % particles", particleHandler.getNumberOfObjects());
    }
};
//
int main(int argc UNUSED, char *argv[] UNUSED)
{
    //logger(INFO,"Simple box for creating particles");
    SinterInitialBed problem;
    //
    Mdouble gravityValue = -9.81;
    Mdouble XMaxDomain = 0.00025;//problem.scale*45e-1+(2.0*problem.scale);
    Mdouble YMaxDomain = 0.00025;//problem.scale*30e-1+(2.0*problem.scale);
    Mdouble ZMaxDomain = 0.01;//problem.scale*40e-1*(4.0);//*problem.scale);
    //
    Mdouble MaxSimTime = 0.03;
    //
    Mdouble MatDensity = 1400.0;
    Mdouble elasticModulus = 8e7;
    Mdouble possionsRation = 0.2;
    Mdouble latentHeat = 245e3;
    Mdouble heatCapacity = 2270;
    Mdouble thermalcond = 0.36;
    Mdouble thermalConvCoff = 15.0;
    Mdouble aTemp = 443.15;//413.15;
    Mdouble mTemp = 451.15;
    Mdouble damping = 0.1;//5e-2;//5e-6;//0.5;//5e-6;
    Mdouble heatingTime = 0.0;//0.0075;
    Mdouble heatSource = 0.0;//2.5e-4;//deposition_problem.heat;//2.5e-4;
    Mdouble emmisivity = 0.9;
    //
    //
    //Mdouble restitutionCoeff = 0.1;
    //Mdouble tc = 90e-7;
    //Mdouble timeStep = 1e-7;//0.00005/50.0;//0.02*tc;
    //
    //-----------------------------------------------------------------------------------------------------------------------------------------------------
    Mdouble SlidingFrictionCoeff = 0.25;//0.25;
    Mdouble RollingFrictionCoeff = 0.05;//0.05;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    //
    problem.setName("SinteringInitialBed");
    problem.setSystemDimensions(3);
    problem.setGravity(Vec3D(0.0,0.0,gravityValue));
    //deposition_problem.setXMin(0.0);
    //deposition_problem.setYMin(0.0);
    //deposition_problem.setZMin(0.0);
    problem.setXMax(XMaxDomain);
    problem.setYMax(YMaxDomain);
    problem.setZMax(ZMaxDomain);
    problem.setTimeMax(MaxSimTime);
    problem.setHGridMaxLevels(2);
    //
    auto species = problem.speciesHandler.copyAndAddObject(MeltableFrictionSpecies());
    //auto species2 = deposition_problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    //auto wallSpecies = deposition_problem.speciesHandler.copyAndAddObject(species2);
    //auto wallParticleSpecies = deposition_problem.speciesHandler.getMixedObject(species2,wallSpecies);
    //
    species->setDensity(MatDensity);
    species->setElasticModulus(elasticModulus);
    species->setPoissonRatio(possionsRation);
    species->setLatentHeat(latentHeat);
    species->setSolidHeatCapacity(heatCapacity);
    species->setThermalConductivityCoefficient(thermalcond);
    species->setThermalConvectionCoefficient(thermalConvCoff);
    species->setDissipation(damping);
    species->setAmbientTemperature(aTemp);
    species->setMeltingTemperature(mTemp);
    species->setMaterialEmissivity(emmisivity);
    //
    species->setHeatingTime(heatingTime);
    species->setHeatInput(heatSource);
    //
    species->setMaxTimeStep(22.0e-6/2.0,MatDensity,elasticModulus);
    //Mdouble mass = species->getMassFromRadius(13.183e-4/2);
    //------------------------------------------------------------------------------Friction coeffiecients Values---------------------------------------------------------------------------
    species->setSlidingStiffness(2./7.* species->getElasticModulus());
    species->setSlidingDissipation(2./7.* species->getDissipation());
    species->setSlidingFrictionCoefficient(SlidingFrictionCoeff);
    //species->setSlidingFrictionCoefficientStatic(0.6);
    //
    species->setRollingStiffness(2./5.* species->getElasticModulus());
    species->setRollingDissipation(2./5.* species->getDissipation());
    species->setRollingFrictionCoefficient(RollingFrictionCoeff);
    //------------------------------------------------------------------------------Data Output--------------------------------------------------------------------------------------------
    problem.setSaveCount(10000);//50
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
