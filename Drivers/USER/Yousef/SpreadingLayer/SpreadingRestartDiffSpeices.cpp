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
#include "Boundaries/PeriodicBoundary.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include <random>
#include "MercuryTime.h"
#include "Boundaries/DeletionBoundary.h"
#include "Species/LinearViscoelasticFrictionJKRAdhesiveSpecies.h"
//#include "Species/LinearViscoelasticFrictionSpecies.h"
//#include "Species/LinearViscoelasticSpecies.h"
//#include "Boundaries/CubeInsertionBoundary.h"
//#include "Chute.h"
//#include <ChuteBottom.h>
//#include "stdio.h"

class SpreadingProcessInsertion: public Mercury3D
{
public:
    // -------------------------------------------------- initial parameters -------------------------------------------------------
    Mdouble TimeRemoveWall = 0.01;
    Mdouble TimeMoveTool = 0.05;
    //
    Mdouble tool_Speed = 0.05;//0.1; //first run -> 0.05 - 50 mm/s;//0.01;// m/s
    //
    //TOOL - initial configuration setup - unit: m:
    Mdouble SpreadingLength = 0.01; // 10 mm -> 1 cm -> 0.01 m;
    Mdouble SpreadingWidth  = 0.001; // 1 mm -> 0.1 cm -> 0.001 m
    Mdouble ToolHight  = 0.003; // 3 mm, 2 mm -> 0.2 cm -> 0.002 m
    Mdouble WallHight = ToolHight;//500e-6;
    Mdouble ToolThickness = 500e-6;// 0.5 mm 250e-6;// m
    Mdouble Gap = 100e-6;//100e-6; //100 microns -> 100e-6 m ;
    //Mdouble TrackLength = 900e-1*(scale*scale);
    //Mdouble Angel = constants::pi/180.0*30;
    //Mdouble AngelDim = scale*4e-1;
//----------------------------------------------------------------------------------------------------
    void setupInitialConditions() override {
        //TOOL
        // ----Balde:
        Vec3D minPoint_tool = Vec3D(0.0,0.0,Gap);
        Vec3D maxPoint_tool = Vec3D(ToolThickness,SpreadingWidth,ToolHight);
        // ----Roller:
        Mdouble rollerRadius = ToolThickness*2.0; // 0.001; // 1 mm
        Vec3D rollerPosition = Vec3D(-ToolThickness,0.0,rollerRadius+Gap);
        //round blade:
        Vec3D minPoint_toolround = Vec3D(0.0,0.0,Gap+ToolThickness);
        Vec3D maxPoint_toolround = Vec3D(ToolThickness,SpreadingWidth,ToolHight);
        //Back Dam
        Vec3D minPoint_back_dam = Vec3D(0.0,0.0,0.0);
        Vec3D maxPoint_back_dam = Vec3D(ToolThickness,SpreadingWidth,Gap);
        //-------------------------------------------------------------------------------particles / PSD-------------------------------------------------------------------
        //------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //Bottom Wall
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0,0.0,0.0));
        wallHandler.copyAndAddObject(w0);
        //TOOL:
        //------Blade:
        IntersectionOfWalls w1;
        w1.setSpecies(speciesHandler.getObject(1)); //using different species for the Tool: w0.setSpecies(speciesHandler.getObject(1));
        w1.addObject(Vec3D(1.0,0.0,0.0),minPoint_tool);
        w1.addObject(Vec3D(0.0,1.0,0.0),minPoint_tool);
        w1.addObject(Vec3D(0.0,0.0,1.0),minPoint_tool);
        w1.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_tool);
        w1.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_tool);
        w1.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_tool); //w1.addObject(Vec3D(cos(Angel),0.0,sin(Angel)),Vec3D(AngelDim,0.0,Gap)); //in next line
        wallHandler.copyAndAddObject(w1);
        //------Roller: w1
/*        AxisymmetricIntersectionOfWalls roller;
        roller.setSpecies(speciesHandler.getObject(0));
        roller.setPosition(rollerPosition);
        roller.setAxis({0, 1, 0});//setAxis
        roller.addObject({-1, 0, 0}, {rollerRadius, 0, 0});
        roller.addObject({0, 0, -1}, {0, 0, ToolThickness});
        roller.addObject({0, 0, 1}, {0, 0, ToolThickness});
        roller.setAngularVelocity({0,-tool_Speed/rollerRadius,0});
        wallHandler.copyAndAddObject(roller);*/
        //-----------Blade with round edge:
/*        IntersectionOfWalls w1;
        w1.setSpecies(speciesHandler.getObject(0)); //using different species for the Tool: w0.setSpecies(speciesHandler.getObject(1));
        w1.addObject(Vec3D(1.0,0.0,0.0),minPoint_toolround);
        w1.addObject(Vec3D(0.0,1.0,0.0),minPoint_toolround);
        w1.addObject(Vec3D(0.0,0.0,1.0),minPoint_toolround);
        w1.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_toolround);
        w1.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_toolround);
        w1.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_toolround);
        wallHandler.copyAndAddObject(w1);
        AxisymmetricIntersectionOfWalls roller;
        roller.setSpecies(speciesHandler.getObject(0));
        roller.setPosition(Vec3D(0.0,0.0,ToolThickness+Gap));
        roller.setAxis({0, 1, 0});//setAxis
        roller.addObject({-1, 0, 0}, {ToolThickness, 0, 0});
        roller.addObject({0, 0, -1}, {0, 0, ToolThickness});
        //roller.addObject({-1, 0, 0}, {ToolThickness, 0, 0});
        wallHandler.copyAndAddObject(roller);*/
        //Back Dam
        IntersectionOfWalls w2;
        w2.setSpecies(speciesHandler.getObject(0));
        w2.addObject(Vec3D(1.0,0.0,0.0),minPoint_back_dam);
        w2.addObject(Vec3D(0.0,1.0,0.0),minPoint_back_dam);
        w2.addObject(Vec3D(0.0,0.0,1.0),minPoint_back_dam);
        w2.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_back_dam);
        w2.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_back_dam);
        w2.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_back_dam);
        wallHandler.copyAndAddObject(w2);
        //
        PeriodicBoundary b0;
        b0.set(Vec3D(0, 1, 0), 0.0, SpreadingWidth);
        boundaryHandler.copyAndAddObject(b0);
        //----------------------------------------------------------------------- Particles insertion ------------------------------------------------------------------
        //autoNumber();
        // read data files
        auto s = std::to_string(DPMBase::readRunNumberFromFile()-1);
        readDataFile("SpreadingProcessInsertion."+s+".data",14);
        //----------------------------------------------------------------------------------------------------------------------------------------------
        // delete particles outside domain of interest:
        DeletionBoundary deletionBoundary1;
        deletionBoundary1.set(Vec3D(1, 0, 0), (SpreadingLength-0.002)+ToolThickness);//(SpreadingLength-0.005)+ToolThickness);
        boundaryHandler.copyAndAddObject(deletionBoundary1);
        DeletionBoundary deletionBoundary2;
        deletionBoundary2.set(Vec3D(-1, 0, 0), ToolThickness);//(SpreadingLength-0.005)+ToolThickness);
        boundaryHandler.copyAndAddObject(deletionBoundary2);
        //
    }
    //
    //----------------------------------------------------------------------------------------------------------------------------------------------
    void actionsAfterTimeStep() override {
        // Start Spreading
       static bool moveTool = false;
        if (!moveTool && getTime()>= TimeMoveTool) {
            //logger(INFO,"Tool is Moving");
            //blade:
            wallHandler.getObject(1)->setVelocity(Vec3D(tool_Speed,0,0));
            //roller:
            //wallHandler.getLastObject()->setVelocity(Vec3D(tool_Speed,0,0));
            //wallHandler.getObject(1)->setAngularVelocity({0,-rollerVelocity_/rollerRadius_,0});
            moveTool = true;
        }
        // fix static particles after spreading to increase computational efficiency
/*        if (moveTool) {
            //blade:
            Vec3D toolPosition = wallHandler.getObject(1)->getPosition();
            //roller:
            //Vec3D toolPosition = wallHandler.getLastObject()->getPosition();
            for (int i = 0; i < particleHandler.getNumberOfObjects(); ++i) {
                Vec3D particlePosition = particleHandler.getObject(i)->getPosition();
                Mdouble kinEnergy = particleHandler.getObject(i)->getKineticEnergy();
                Mdouble elaEnergy = particleHandler.getObject(i)->getE
                if (particlePosition.X < toolPosition.X &&  < 2e-10) {
                    particleHandler.getObject(i)->fixParticle();
                    //logger(INFO, "particles fixed behind tool");
                }
            }
        }*/
    }
};
// --------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc UNUSED, char *argv[] UNUSED)
{
    Time time;
    time.tic();

    //logger(INFO,"Simple box for creating particles");
    SpreadingProcessInsertion deposition_problem;
    //
    //------------------------------------------------------------------------------ Domain ------------------------------------------------------------------------------
    Mdouble gravityValue = -9.81; // g=9.81 m/s2 - 981 cm/s2
    Mdouble XMaxDomain = deposition_problem.SpreadingLength+deposition_problem.ToolThickness;
    //(deposition_problem.SpreadingLength+deposition_problem.ToolThickness)*2.0+deposition_problem.ToolThickness;
    Mdouble YMaxDomain = deposition_problem.SpreadingWidth;
    Mdouble ZMinDomain = -deposition_problem.ToolThickness;
    Mdouble ZMaxDomain = deposition_problem.ToolHight;
    //------------------------------------------------------------------------------ Parameters ------------------------------------------------------------------------------
    const char *simName = "SpreadingToolSpecies";
    //
    Mdouble MaxSimTime = 0.3;//0.2; // at speed = 50mm/s -> 0.5, but still too long make it 0.3;
    //
    Mdouble MatDensity = 4430;// density of Ti = 4.430g/cm3 - 4430kg/m3
    //
    Mdouble tc = 90e-7; //calculated as tc = 0.005-0.01 * tg : tg=sqrt(d50/g) => 0.005*sqrt(35e-6/9.81)
    //
    Mdouble restitutionCoeff = 0.4;//0.1;
    //
    Mdouble timeStep = 0.02*tc;
    //
    Mdouble minRaduis = 12.0e-6;//15.0e-6/2.0;//13.183e-6/2.0;
    //
    unsigned int saveCount = 10000; // 5000
    //
    Mdouble AdhStiffnessFactor = 0.5;
    Mdouble surfaceEnergy = 0.1e-3; //J/m2
    //------------------------------------------------------------------------------Parametric friction:------------------------------------------------------------------------------
    deposition_problem.autoNumber();
    std::vector<int> studyNum=deposition_problem.get2DParametersFromRunNumber(8,8);
    //
    std::vector<double> SlidingFrictionCoeffVector = {10.0,2.0,1.0,0.5,0.4,0.25,0.1,0.05};
    std::vector<double> RollingFrictionCoeffVector = {10.0,0.4,0.2,0.15,0.1,0.075,0.05,0.005};
    //
    Mdouble SlidingFrictionCoeff = SlidingFrictionCoeffVector[studyNum[1]-1];
    Mdouble RollingFrictionCoeff = RollingFrictionCoeffVector[studyNum[2]-1];
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    //
    deposition_problem.setName(simName);
    deposition_problem.setSystemDimensions(3);
    deposition_problem.setGravity(Vec3D(0.0,0.0,gravityValue));
    //deposition_problem.setXMin(0.0);
    //deposition_problem.setYMin(0.0);
    deposition_problem.setZMin(ZMinDomain);
    deposition_problem.setXMax(XMaxDomain);
    deposition_problem.setYMax(YMaxDomain);
    deposition_problem.setZMax(ZMaxDomain);
    deposition_problem.setTimeMax(MaxSimTime);
    //deposition_problem.setHGridMaxLevels(2);
    //
    //auto species = deposition_problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    auto species = deposition_problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionJKRAdhesiveSpecies());
    //
    species->setDensity(MatDensity);
    //
    // Here we assign the mass of the smallest particle to calculate the stiffiness and disppation as we want to get ------
    Mdouble mass = species->getMassFromRadius(minRaduis);
    // ----- set stiffness and dissipation from tc,COR,mass
    species->setCollisionTimeAndRestitutionCoefficient(tc,restitutionCoeff,mass);
    //------------------------------------------------------------------------------ Cohesion Parameters ------------------------------------
    Mdouble adhesionStiffness = AdhStiffnessFactor*species->getStiffness();
    species->setAdhesionStiffness(adhesionStiffness);
    species->setSurfaceEnergy(surfaceEnergy);
    //------------------------------------------------------------------------------Sliding and rolling stiffness and dissipation -------------------------------------
    species->setSlidingStiffness(2./7.*species->getStiffness());
    species->setSlidingDissipation(2./7.*species->getDissipation());
    species->setSlidingFrictionCoefficient(SlidingFrictionCoeff);
    //species->setSlidingFrictionCoefficientStatic(0.6);
    //
    species->setRollingStiffness(2./5.*species->getStiffness());
    species->setRollingDissipation(2./5.*species->getDissipation());
    species->setRollingFrictionCoefficient(RollingFrictionCoeff);
    //------------------------------------------------------------------------------ P-Tool species -------------------------------------------------------------
/*    auto toolSpecies = deposition_problem.speciesHandler.copyAndAddObject(species);
    auto toolParticleSpecies = deposition_problem.speciesHandler.getMixedObject(species,toolSpecies);
    toolParticleSpecies->setSlidingFrictionCoefficient(0.5*species->getSlidingFrictionCoefficient());
    toolParticleSpecies->setRollingFrictionCoefficient(0.5*species->getRollingFrictionCoefficient());
    //toolParticleSpecies->setAdhesionStiffness(0.0);
    toolParticleSpecies->setSurfaceEnergy(0.0);*/
    //------------------------------------------------------------------------------ Tool species -------------------------------------------------------------
    auto toolSpecies = deposition_problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionJKRAdhesiveSpecies());
    toolSpecies->setDensity(MatDensity);
    Mdouble massPT = toolSpecies->getMassFromRadius(minRaduis);
    // ----- set stiffness and dissipation from tc,COR,mass
    toolSpecies->setCollisionTimeAndRestitutionCoefficient(tc,restitutionCoeff,massPT);
    //
    toolSpecies->setSlidingStiffness(2./7.*species->getStiffness());
    toolSpecies->setSlidingDissipation(2./7.*species->getDissipation());
    toolSpecies->setSlidingFrictionCoefficient(species->getSlidingFrictionCoefficient()/4.0);
    toolSpecies->setRollingStiffness(2./5.*species->getStiffness());
    toolSpecies->setRollingDissipation(2./5.*species->getDissipation());
    toolSpecies->setRollingFrictionCoefficient(species->getRollingFrictionCoefficient()/4.0);
    //------------------------------------------------------------------------------ P-Tool species -------------------------------------------------------------
    auto toolParticleSpecies = deposition_problem.speciesHandler.getMixedObject(species,toolSpecies);
    //
    toolParticleSpecies->setStiffness(species->getStiffness());
    toolParticleSpecies->setDissipation(species->getDissipation());
    //
    toolParticleSpecies->setSlidingStiffness(2./7.*species->getStiffness());
    toolParticleSpecies->setSlidingDissipation(2./7.*species->getDissipation());
    toolParticleSpecies->setSlidingFrictionCoefficient(species->getSlidingFrictionCoefficient()/4.0);
    toolParticleSpecies->setRollingStiffness(2./5.*species->getStiffness());
    toolParticleSpecies->setRollingDissipation(2./5.*species->getDissipation());
    toolParticleSpecies->setRollingFrictionCoefficient(species->getRollingFrictionCoefficient()/4.0);
    //
    //--------------------------------------------------------------------------------Parallel-------------------------------------------------------------------------
    //deposition_problem.setNumberOfDomains({1,2,1});
    //------------------------------------------------------------------------------Data Output--------------------------------------------------------------------------------------------
    deposition_problem.setSaveCount(saveCount);//50
    deposition_problem.dataFile.setFileType(FileType::ONE_FILE);
    //deposition_problem.fStatFile.setFileType(FileType::NO_FILE);
    deposition_problem.setParticlesWriteVTK(true);
    deposition_problem.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    //
    deposition_problem.setTimeStep(timeStep);
    //
    deposition_problem.solve();//(argc, argv);
    //
    // record time:
    std::ofstream outTime ("simTime"+deposition_problem.getName(), std::ios::out);
    outTime << time.toc() << std::endl;
    outTime.close();
    //logger(INFO, "Total time to run this simulation: % s", time.toc());
    //
    return 0;
}