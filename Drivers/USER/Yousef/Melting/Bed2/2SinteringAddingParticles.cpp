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


class SinteringAddingParticles: public Mercury3D
{
public:
    //
    SinteringAddingParticles ()
    {
        setName("SinteringUniform");
        readRestartFile();
        //setRestarted(false);
        //
        Mdouble ToolLength = 25e-5;//0.00025; //45e-1 * scale * 2.0; // unit cm
        Mdouble ToolWidth = 25e-5;//0.00025;//30e-1 * scale * 2.0;
        Mdouble ToolHight = 0.01;//40e-1 * scale * 4.0; //20e-1 * scale * 4.0;
        Mdouble ToolThickness = 25e-6;//0.000025;//5e-1 * scale;
        Mdouble Gap = 0.1e-3;
        //Front Dam
        Mdouble damThickness = ToolLength + ToolThickness;//(5.0e-1 * scale);
        //particles:
        Mdouble minRadius = 40e-6;
        Mdouble maxRadius = 50e-6;
        Mdouble meanRadius1 = 0.5*(minRadius+maxRadius);
        //
        //the volume to be added:
        Mdouble addVolume = Gap*ToolLength*ToolWidth;//(20e-1 * scale)*ToolLength*ToolWidth;//(ToolHight/4.0)*ToolLength*ToolWidth; //0.0018 cm3
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
        ThermalParticle p1;
        p1.setSpecies(speciesHandler.getObject(0));
        p1.setRadius(meanRadius1);
        p1.setTemperature(443.15);
        if(p1.getKineticEnergy() < 1e-5)
        {
            p1.setLocalHeatingTime(0.03 + 1.5 + 0.04);
            p1.setLocalHeatInput(1e-3);
        }
        Mdouble fillHeight = 0.0;
        while (addVolume > 0) {
            Mdouble x = random.getRandomNumber(ToolThickness, ToolLength);
            Mdouble y = random.getRandomNumber(ToolThickness, (ToolWidth-ToolThickness));
            Mdouble z = random.getRandomNumber(25e-5, fillHeight+25e-5);
            p1.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p1)) {
                particleHandler.copyAndAddObject(p1);
                addVolume -= p1.getVolume();
                do {
                    p1.setRadius(d(gen));
                } while (p1.getRadius() < minRadius || p1.getRadius() > maxRadius);
                if (particleHandler.getNumberOfObjects() % 100 == 0) std::cout << '.' << std::flush;
                if (particleHandler.getNumberOfObjects() % 1000 == 0) std::cout << ' ';
                if (particleHandler.getNumberOfObjects() % 10000 == 0) std::cout << addVolume << '\n';
            } else {
                fillHeight += 0.01*meanRadius1;//0.01 * scale * meanRadius1; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        logger(INFO, " Inserted % particles", particleHandler.getNumberOfObjects());
        //
        //
        setName("SinteringAddingParticles");
        //species = speciesHandler.copyAndAddObject(MeltableFrictionSpecies());
        species = dynamic_cast<MeltableFrictionSpecies *>(speciesHandler.getObject(0));
    };

    //
/*    void setupInitialConditions() {

        //TOOL:
        Mdouble ToolLength = 25e-5;//0.00025; //45e-1 * scale * 2.0; // unit cm
        Mdouble ToolWidth = 25e-5;//0.00025;//30e-1 * scale * 2.0;
        Mdouble ToolHight = 0.01;//40e-1 * scale * 4.0; //20e-1 * scale * 4.0;
        Mdouble ToolThickness = 25e-6;//0.000025;//5e-1 * scale;
        Mdouble Gap = 0.1e-3;
        //Front Dam
        Mdouble damThickness = ToolLength + ToolThickness;//(5.0e-1 * scale);
        //particles:
        Mdouble minRadius = 40e-6;
        Mdouble maxRadius = 50e-6;
        Mdouble meanRadius1 = 0.5*(minRadius+maxRadius);
        //
        //the volume to be added:
        Mdouble addVolume = Gap*ToolLength*ToolWidth;//(20e-1 * scale)*ToolLength*ToolWidth;//(ToolHight/4.0)*ToolLength*ToolWidth; //0.0018 cm3
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
        ThermalParticle p1;
        p1.setSpecies(speciesHandler.getObject(0));
        p1.setRadius(meanRadius1);
        p1.setTemperature(443.15);
        Mdouble fillHeight = 0.0;
        while (addVolume > 0) {
            Mdouble x = random.getRandomNumber(ToolThickness, ToolLength);
            Mdouble y = random.getRandomNumber(ToolThickness, (ToolWidth-ToolThickness));
            Mdouble z = random.getRandomNumber(35e-5, fillHeight);
            p1.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p1)) {
                particleHandler.copyAndAddObject(p1);
                addVolume -= p1.getVolume();
                do {
                    p1.setRadius(d(gen));
                } while (p1.getRadius() < minRadius || p1.getRadius() > maxRadius);
                if (particleHandler.getNumberOfObjects() % 100 == 0) std::cout << '.' << std::flush;
                if (particleHandler.getNumberOfObjects() % 1000 == 0) std::cout << ' ';
                if (particleHandler.getNumberOfObjects() % 10000 == 0) std::cout << addVolume << '\n';
            } else {
                fillHeight += 0.01*meanRadius1;//0.01 * scale * meanRadius1; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        logger(INFO, " Inserted % particles", particleHandler.getNumberOfObjects());
    }*/

    MeltableFrictionSpecies *species;
};
//
int main(int argc UNUSED, char *argv[] UNUSED)
{
    //logger(INFO,"Simple box for creating particles");
    SinteringAddingParticles problem;
    //
    Mdouble gravityValue = -9.81;
    Mdouble XMaxDomain = 25e-5;//0.00025;//problem.scale*45e-1+(2.0*problem.scale);
    Mdouble YMaxDomain = 25e-5;//0.00025;//problem.scale*30e-1+(2.0*problem.scale);
    Mdouble ZMaxDomain = 0.01;//problem.scale*40e-1*(4.0);//*problem.scale);
    //
    Mdouble MaxSimTime = 0.03+1.5+14.0;
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
/*    Mdouble SlidingFrictionCoeff = 0.25;//0.25;
    Mdouble RollingFrictionCoeff = 0.05;//0.05;*/
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    //
    problem.setName("SinteringAddingParticles");
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
    //auto species = problem.speciesHandler.copyAndAddObject(MeltableFrictionSpecies());
    //auto species2 = deposition_problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    //auto wallSpecies = deposition_problem.speciesHandler.copyAndAddObject(species2);
    //auto wallParticleSpecies = deposition_problem.speciesHandler.getMixedObject(species2,wallSpecies);
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
    //Mdouble mass = species->getMassFromRadius(13.183e-4/2);
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
    problem.setSaveCount(10000);//50
    //problem.restartFile.setSaveCount(10000);
    //deposition_problem.dataFile.setFileType(FileType::NO_FILE);javascript:void(0)
    //deposition_problem.fStatFile.setFileType(FileType::NO_FILE);
    problem.setParticlesWriteVTK(true);
    problem.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    //
    problem.setTimeStep(problem.species->getMaxTimeStep());
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    problem.solve();
    return 0;
}
