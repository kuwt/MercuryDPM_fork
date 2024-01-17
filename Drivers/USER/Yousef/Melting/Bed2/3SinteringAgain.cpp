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


class SinteringAgain: public Mercury3D
{
public:
    //
    //Mdouble scale = 5e-3; //1e-1; //copy to speed and simulation boundaries
    //
    //Mdouble startHeating = 0.035;//0.2;
    //Mdouble stopHeating = 0.0075;//0.25;
    //
    SinteringAgain ()
    {
        setName("SinteringAddingParticles");
        readRestartFile();
        //setRestarted(false);
        //
        ThermalParticle p1;
        p1.setLocalHeatingTime(0.03 + 1.5 + 0.04);
        p1.setLocalHeatInput(1e-3);
        //
        setName("SinteringAgain");
        species = dynamic_cast<MeltableFrictionSpecies *>(speciesHandler.getObject(0));
    };
    //
    void setupInitialConditions() {

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
    SinteringAgain problem;
    //
    Mdouble gravityValue = -9.81;
    Mdouble XMaxDomain = 25e-5;//0.00025;//problem.scale*45e-1+(2.0*problem.scale);
    Mdouble YMaxDomain = 25e-5;//0.00025;//problem.scale*30e-1+(2.0*problem.scale);
    Mdouble ZMaxDomain = 0.01;//problem.scale*40e-1*(4.0);//*problem.scale); Mdouble ZMaxDomain = problem.scale*20e-1+(4.0*problem.scale);
    //
    Mdouble MaxSimTime = 0.03+1.5+14.0;
    //
    Mdouble MatDensity = 1400.0;
    Mdouble elasticModulus = 8e7;
    Mdouble possionsRation = 0.2;
    Mdouble latentHeat = 245e3;
    Mdouble heatCapacity = 2270;
    Mdouble thermalcond = 0.36;
    Mdouble thermalConvCoff = 15;
    Mdouble aTemp = 443.15;//413.15;
    Mdouble mTemp = 451.15;
    Mdouble damping = 0.1;
    Mdouble heatingTime = 0.0;//0.03+1.5+0.04;//0.03 ; 0.06;//0.075;//1.50;//0.0075;
    Mdouble heatSource = 0.0;//1e-3;//2.5e-4;//10.0;//2.5e-4;//problem.heat;//2.5e-4;
    Mdouble emmisivity = 0.9;
    //
    //
    //Mdouble restitutionCoeff = 0.1;
    //Mdouble tc = 90e-7;
    //Mdouble timeStep = 1e-7;//0.00005/50.0;//0.02*tc;
    //
    //-----------------------------------------------------------------------------------------------------------------------------------------------------
/*    Mdouble SlidingFrictionCoeff = 0.25;
    Mdouble RollingFrictionCoeff = 0.05;*/
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    //
    problem.setName("SinteringAgain");
    problem.setSystemDimensions(3);
    problem.setGravity(Vec3D(0.0,0.0,gravityValue));
    //problem.setXMin(0.0);
    //problem.setYMin(0.0);
    //problem.setZMin(0.0);
    problem.setXMax(XMaxDomain);
    problem.setYMax(YMaxDomain);
    problem.setZMax(ZMaxDomain);
    problem.setTimeMax(MaxSimTime);
    problem.setHGridMaxLevels(2);
    //
    //auto species = problem.speciesHandler.copyAndAddObject(MeltableFrictionSpecies());
    //auto species2 = problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    //auto wallSpecies = problem.speciesHandler.copyAndAddObject(species2);
    //auto wallParticleSpecies = problem.speciesHandler.getMixedObject(species2,wallSpecies);
    //
    problem.species->setDensity(MatDensity);
    problem.species->setElasticModulus(elasticModulus);
    problem.species->setPoissonRatio(possionsRation);
    problem.species->setLatentHeat(latentHeat);
    problem.species->setSolidHeatCapacity(heatCapacity);
    problem.species->setThermalConductivityCoefficient(thermalcond);
    problem.species->setThermalConvectionCoefficient(thermalConvCoff);
    problem.species->setDissipation(damping);
    problem.species->setAmbientTemperature(aTemp);
    problem.species->setMeltingTemperature(mTemp);
    problem.species->setMaterialEmissivity(emmisivity);
    //
    problem.species->setHeatingTime(heatingTime);
    problem.species->setHeatInput(heatSource);
    //
    problem.species->setMaxTimeStep(40e-6,MatDensity,elasticModulus);
    //
    //Mdouble mass = species->getMassFromRadius(13.183e-4/2);
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
    problem.setSaveCount(10000);//50
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
