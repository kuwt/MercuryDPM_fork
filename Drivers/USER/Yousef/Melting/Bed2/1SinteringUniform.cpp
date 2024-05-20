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


class Sintering: public Mercury3D
{
public:
    //
    Mdouble bottomWallZ = -3e-4;
    bool laserON = true;
    Mdouble maxSimTime = 6.0;//6.0;
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
    Mdouble heatingTime = 0.075;
    Mdouble coolingTime = 0.075 + 0.001;
    //
    Mdouble heatSource = 0.0;//384e-6/1e-3;
    //
    //bool resetLaser = true;
    Mdouble absorptivity = 0.33;
    Mdouble laserPower = 5;//1;
    Mdouble beamRadius = 250e-6;//125e-6;
    Mdouble powderBedPorosity = 0.55;
    Mdouble initialLaserPX = beamRadius;//-beamRadius;
    Mdouble laserSpeedX = 0.0;//0.1;//-0.1;//100e-6;
    Mdouble initialLaserPY = beamRadius;//-beamRadius;//beamRadius * 2.0 + (beamRadius / 2.0); // second pass -> beamRadius*2.0;// first scan -> beamRadius;
    Mdouble laserSpeedY = 0.0;//0.1;//0.0;//0.1;
    //------------------------------------------------------------------------------------------------
    //particles:
    Mdouble minRadius = 40e-6;
    Mdouble maxRadius = 50e-6;
    Mdouble meanRadius1 = 0.5*(minRadius+maxRadius);
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------
    //
    //Mdouble scale = 5e-3; //1e-1; //copy to speed and simulation boundaries
    //
    //Mdouble startHeating = 0.035;//0.2;
    //Mdouble stopHeating = 0.0075;//0.25;
    //
    Sintering ()
    {
        setName("SinteringInitialBedUniform");
        readRestartFile();
        //setRestarted(false);
        setName("SinteringUniform");
        species = dynamic_cast<MeltableFrictionSpecies *>(speciesHandler.getObject(0));
    };
    //
    void setupInitialConditions() override {

        //TOOL:
/*        Mdouble ToolLength = 0.025;//45e-1 * scale * 2.0; // unit cm
        Mdouble ToolWidth = 0.025;//30e-1 * scale * 2.0;
        Mdouble ToolHight = 0.08;//20e-1 * scale * 4.0;
        Mdouble ToolThickness = 0.0025;//5e-1 * scale;
        Mdouble Gap = 0.1e-1;
        //Front Dam
        Mdouble damThickness = ToolLength + ToolThickness;//(5.0e-1 * scale);*/
    }
    MeltableFrictionSpecies *species;
};
//
int main(int argc UNUSED, char *argv[] UNUSED)
{
    //logger(INFO,"Simple box for creating particles");
    Sintering problem;
    //
    Mdouble gravityValue = -9.81;
    Mdouble XMaxDomain = 25e-5;//0.00025;//problem.scale*45e-1+(2.0*problem.scale);
    Mdouble YMaxDomain = 25e-5;//0.00025;//problem.scale*30e-1+(2.0*problem.scale);
    Mdouble ZMaxDomain = 0.01;//problem.scale*40e-1*(4.0);//*problem.scale); Mdouble ZMaxDomain = problem.scale*20e-1+(4.0*problem.scale);
    //
    //-----------------------------------------------------------------------------------------------------------------------------------------------------
/*    Mdouble SlidingFrictionCoeff = 0.25;
    Mdouble RollingFrictionCoeff = 0.05;*/
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    //
    problem.setName("SinteringUniform");
    problem.setSystemDimensions(3);
    problem.setGravity(Vec3D(0.0,0.0,gravityValue));
    //problem.setXMin(0.0);
    //problem.setYMin(0.0);
    problem.setZMin(problem.bottomWallZ);
    problem.setXMax(XMaxDomain);
    problem.setYMax(YMaxDomain);
    problem.setZMax(ZMaxDomain);
    problem.setTimeMax(problem.maxSimTime);
    //problem.setHGridMaxLevels(2);
    //
    //auto species = problem.speciesHandler.copyAndAddObject(MeltableFrictionSpecies());
    //auto species2 = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    //auto wallSpecies = problem.speciesHandler.copyAndAddObject(species2);
    //auto wallParticleSpecies = problem.speciesHandler.getMixedObject(species2,wallSpecies);
    //
    problem.species->setDensity(problem.MatDensity);
    problem.species->setElasticModulus(problem.elasticModulus);
    problem.species->setPoissonRatio(problem.possionsRation);
    problem.species->setDissipation(problem.damping);
    //
    problem.species->setLatentHeat(problem.latentHeat);
    problem.species->setSolidHeatCapacity(problem.heatCapacity);
    problem.species->setThermalConductivityCoefficient(problem.thermalcond);
    problem.species->setThermalConvectionCoefficient(problem.thermalConvCoff);
    problem.species->setAmbientTemperature(problem.aTemp);
    problem.species->setMeltingTemperature(problem.mTemp);
    problem.species->setMaterialEmissivity(problem.emmisivity);
    //
    problem.species->setSurfaceTension(problem.surfaceTension);
    problem.species->setRefViscosity(problem.refVis);
    //
    problem.species->setDeltaT(problem.deltaT);
    problem.species->setLiquidHeatCapacity(problem.liquidHeatCapacity);
    problem.species->setVaporizationHeatCapacity(problem.vaporizationHeatCapacity);
    problem.species->setVaporizationLatentHeat(problem.vaporizationLatentHeat);
    problem.species->setVaporizationTemperature(problem.vaporizationTemp);
    //
    problem.species->setThermalExpansionCoefficient(problem.thermalExpansionCoeff);
    problem.species->setThermalExpansionOnOff(problem.thermalExpansionOnOff);
    //
/*    species->setHeatingTime(problem.heatingTime);
    species->setHeatInput(problem.heatSource);
    species->setCoolingTime(problem.coolingTime);*/
    if(problem.laserON)
    {
        // Laser
        problem.species->setHeatingTime(problem.heatingTime);
        problem.species->setCoolingTime(problem.coolingTime);
        //
        problem.species->setMaterialAbsorptivity(problem.absorptivity);
        problem.species->setLaserPower(problem.laserPower);
        problem.species->setBeamSpotRadius(problem.beamRadius);
        problem.species->setPowderBedPorosity(problem.powderBedPorosity);
        problem.species->setInitialLaserPositionX(problem.initialLaserPX);
        problem.species->setLaserSpeedX(problem.laserSpeedX);
        problem.species->setInitialLaserPositionY(problem.initialLaserPY);
        problem.species->setLaserSpeedY(problem.laserSpeedY);
    }
    else if (!problem.laserON)
    {
        problem.species->setHeatingTime(problem.heatingTime);
        problem.species->setHeatInput(problem.heatSource);
        problem.species->setCoolingTime(problem.coolingTime);
    }
    else
    {
        problem.species->setHeatInput(0.0);
    }
    //
    problem.species->setMaxTimeStep(problem.minRadius,problem.MatDensity,problem.elasticModulus);
    //------------------------------------------------------------------------------Friction coeffiecients Values---------------------------------------------------------------------------
/*    problem.species->setSlidingStiffness(2./7.*problem.species->getEffectiveElasticModulus());
    problem.species->setSlidingDissipation(2./7.*problem.species->getViscosityCoefficient());
    problem.species->setSlidingFrictionCoefficient(SlidingFrictionCoeff);
    //species->setSlidingFrictionCoefficientStatic(0.6);
    //
    problem.species->setRollingStiffness(2./5.*problem.species->getElasticModulus());
    problem.species->setRollingDissipation(2./5.*problem.species->getDissipation());
    problem.species->setRollingFrictionCoefficient(RollingFrictionCoeff);*/
    //------------------------------------------------------------------------------Data Output--------------------------------------------------------------------------------------------
    problem.setSaveCount(problem.numberOfCounts);//50
    //problem.dataFile.setFileType(FileType::NO_FILE);javascript:void(0)
    //problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.setParticlesWriteVTK(true);
    //problem.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    //
    problem.setTimeStep(problem.species->getMaxTimeStep());
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    problem.solve();
    return 0;
}
