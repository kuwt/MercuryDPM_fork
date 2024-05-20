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

class SpreadingProcessRestart: public Mercury3D
{
public:
    // -------------------------------------------------- initial parameters -------------------------------------------------------
    Mdouble TimeRemoveWall = 0.01;
    Mdouble TimeMoveTool = 0.05;
    //
    Mdouble tool_Speed = 0.01;//0.1 - 100 mm/s,  0.05 - 50 mm/s , 0.01 - 10 mm/s
    //
    //TOOL - initial configuration setup - unit: m:
    Mdouble SpreadingLength = 0.01; // 10 mm -> 1 cm -> 0.01 m;
    Mdouble SpreadingWidth  = 0.001; // 1 mm -> 0.1 cm -> 0.001 m
    Mdouble ToolHight  = 0.003; // 3 mm, 2 mm -> 0.2 cm -> 0.002 m
    Mdouble WallHight = ToolHight;//500e-6;
    Mdouble ToolThickness = 500e-6;// 0.5 mm 250e-6;// m
    Mdouble Gap = 100e-6; //100 microns -> 100e-6 m, 180 microns -> 180e-6 m ;
    bool changeGapHeight = false;
    //
    //Mdouble Angel = constants::pi/180.0*30;
    //Mdouble AngelDim = scale*4e-1;
    //
    //----------------------------------------------------------------------------------------------------
    SpreadingProcessRestart(){
        autoNumber();
        // restart all files
        auto s = std::to_string(DPMBase::readRunNumberFromFile()-1);
        // restart from specific files
        //std::vector<int> set1 = {19,20,21,22,23,24,27,28,29,30,31,32,35,36,37,38,39,40,43,44,45,46,47,48,51,52,53,54,55,56,59,60,61,62,63,64};
        //std::vector<int> set2 = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,25,26,33,34,41,42,49,50,57,58};
        //std::vector<int> set3 = {11,12,13,14,15,16};
        //
        //std::vector<int> setSlow = {19,20,21,22,23,24,59,60,61,62,63,64};
        //
        //auto fileNumber = std::to_string(set3[DPMBase::readRunNumberFromFile()-2]);
        // start from restart
        setName("SpreadingProcessInsertion."+s);//fileNumber);
        readRestartFile();
        //setRestarted(false);
        setName("SpreadingProcessRestart."+s);//fileNumber);
        //species = dynamic_cast<LinearViscoelasticFrictionJKRAdhesiveSpecies *>(speciesHandler.getObject(0));
        //----------------------------------------------------------------------------------------------------
        if (changeGapHeight) {
            wallHandler.removeObject(2);
            wallHandler.removeObject(1);
            //----------------------------------------------------------------------------------------------------
            // ----Balde:
            Vec3D minPoint_tool = Vec3D(0.0, 0.0, Gap);
            Vec3D maxPoint_tool = Vec3D(ToolThickness, SpreadingWidth, ToolHight);
            //Back Dam 2
            Vec3D minPoint_back_dam = Vec3D(0.0, 0.0, 0.0);
            Vec3D maxPoint_back_dam = Vec3D(ToolThickness, SpreadingWidth, Gap);
            //Bottom Wall
/*            InfiniteWall w0;
            w0.setSpecies(speciesHandler.getObject(0));
            w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, 0.0));
            wallHandler.copyAndAddObject(w0);*/
            //------Blade:
            IntersectionOfWalls w1;
            w1.setSpecies(speciesHandler.getObject(0)); //using different species for the Tool: w0.setSpecies(speciesHandler.getObject(1));
            w1.addObject(Vec3D(1.0, 0.0, 0.0), minPoint_tool);
            w1.addObject(Vec3D(0.0, 1.0, 0.0), minPoint_tool);
            w1.addObject(Vec3D(0.0, 0.0, 1.0), minPoint_tool);
            w1.addObject(Vec3D(-1.0, 0.0, 0.0), maxPoint_tool);
            w1.addObject(Vec3D(0.0, -1.0, 0.0), maxPoint_tool);
            w1.addObject(Vec3D(0.0, 0.0, -1.0),maxPoint_tool); //w1.addObject(Vec3D(cos(Angel),0.0,sin(Angel)),Vec3D(AngelDim,0.0,Gap)); //in next line
            w1.setVelocity(Vec3D(tool_Speed, 0, 0));
            wallHandler.copyAndAddObject(w1);
            //Back Dam
            IntersectionOfWalls w2;
            w2.setSpecies(speciesHandler.getObject(0));
            w2.addObject(Vec3D(1.0, 0.0, 0.0), minPoint_back_dam);
            w2.addObject(Vec3D(0.0, 1.0, 0.0), minPoint_back_dam);
            w2.addObject(Vec3D(0.0, 0.0, 1.0), minPoint_back_dam);
            w2.addObject(Vec3D(-1.0, 0.0, 0.0), maxPoint_back_dam);
            w2.addObject(Vec3D(0.0, -1.0, 0.0), maxPoint_back_dam);
            w2.addObject(Vec3D(0.0, 0.0, -1.0), maxPoint_back_dam);
            wallHandler.copyAndAddObject(w2);
        }
    };
    //----------------------------------------------------------------------------------------------------
    void setupInitialConditions() override {
        //
        //
        PeriodicBoundary b0;
        b0.set(Vec3D(0, 1, 0), 0.0, SpreadingWidth);
        boundaryHandler.copyAndAddObject(b0);
        //----------------------------------------------------------------------- Particles insertion ------------------------------------------------------------------
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
    // Remove dams
    void actionsAfterTimeStep() override {

        // Start Spreading
        static bool moveTool = false;
        if (!changeGapHeight && !moveTool && getTime()>= TimeMoveTool) {
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
                if (particlePosition.X < toolPosition.X &&  kinEnergy < 6e-12) {
                    particleHandler.getObject(i)->fixParticle();
                    //logger(INFO, "particles fixed behind tool");
                }
            }
        }*/
    }
    //
    /*void printTime() const override
    {
        logger(INFO,"t=%, tMax=%, N=%", getTime(),getTimeMax(), particleHandler.getSize());
    }*/
};
// --------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc UNUSED, char *argv[] UNUSED)
{
    Time time;
    time.tic();

    //logger(INFO,"Simple box for creating particles");
    SpreadingProcessRestart deposition_problem;
    //
    //------------------------------------------------------------------------------ Domain ------------------------------------------------------------------------------
    Mdouble gravityValue = -9.81; // g=9.81 m/s2 - 981 cm/s2
    Mdouble XMaxDomain = deposition_problem.SpreadingLength+deposition_problem.ToolThickness;
    //(deposition_problem.SpreadingLength+deposition_problem.ToolThickness)*2.0+deposition_problem.ToolThickness;
    Mdouble YMaxDomain = deposition_problem.SpreadingWidth;
    Mdouble ZMinDomain = -deposition_problem.ToolThickness;
    Mdouble ZMaxDomain = deposition_problem.ToolHight;
    //------------------------------------------------------------------------------ Parameters ------------------------------------------------------------------------------
    //const char *simName = "SpreadingProcessRestart";
    //
    Mdouble MaxSimTime = 0.3; // at speed = 100 mm/s -> 0.2 , = 50mm/s -> 0.3, = 10 mm/s -> 0.9
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
    //Mdouble AdhStiffnessFactor = 0.5;
    //Mdouble surfaceEnergy = 0.1e-3; //J/m2
    //static bool activeParallel = false;
    //------------------------------------------------------------------------------Parametric friction:------------------------------------------------------------------------------
    //deposition_problem.autoNumber();
    std::vector<int> studyNum=deposition_problem.get2DParametersFromRunNumber(8,8);
    //
    std::vector<double> SlidingFrictionCoeffVector = {10.0,2.0,1.0,0.5,0.4,0.25,0.1,0.05};
    std::vector<double> RollingFrictionCoeffVector = {10.0,0.4,0.2,0.15,0.1,0.075,0.05,0.005};
    //
/*    std::vector<double> SlidingFrictionCoeffVector = {1.0,0.5,0.4,0.25,0.1,0.05};
    std::vector<double> RollingFrictionCoeffVector = {0.2,0.15,0.1,0.075,0.05,0.005};*/
    //
/*    std::vector<double> SlidingFrictionCoeffVector = {1.0,0.5,0.4,0.25,0.1,0.05};
    std::vector<double> RollingFrictionCoeffVector = {0.4};*/
    //
    Mdouble SlidingFrictionCoeff = SlidingFrictionCoeffVector[studyNum[1]-1];
    Mdouble RollingFrictionCoeff = RollingFrictionCoeffVector[studyNum[2]-1];
    //TOOL -> P-Tool
    //Mdouble WSlidingFriCoeff = 0.4;
    //Mdouble WRollingFriCoeff = 0.04;
    //
    //std::cout << "studyNum2=" << studyNum[1] << std::endl;
    //std::cout << "Sliding Friction" << SlidingFrictionCoeff << std::endl;
    //std::cout << "Rolling Friction" << RollingFrictionCoeff << std::endl;
    //std::cout << "Wall Sliding Friction" << WSlidingFriCoeff << std::endl;
    //std::cout << "Wall Rolling Friction" << WRollingFriCoeff << std::endl;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    //
    //deposition_problem.setName(simName);
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
    //species->setStiffness(10000);
    //species->setDissipation(100);
    //
    // Here we assign the mass of the smallest particle to calculate the stiffiness and disppation as we want to get ------
    Mdouble mass = species->getMassFromRadius(minRaduis);
    // ----- set stiffness and dissipation from tc,COR,mass
    species->setCollisionTimeAndRestitutionCoefficient(tc,restitutionCoeff,mass);
    //------------------------------------------------------------------------------ Cohesion Parameters ------------------------------------
    //Mdouble adhesionStiffness = AdhStiffnessFactor*species->getStiffness();
    //species->setAdhesionStiffness(adhesionStiffness);
    //species->setSurfaceEnergy(surfaceEnergy);
    //------------------------------------------------------------------------------Sliding and rolling stiffness and dissipation -------------------------------------
    species->setSlidingStiffness(2./7.*species->getStiffness());
    species->setSlidingDissipation(2./7.*species->getDissipation());
    species->setSlidingFrictionCoefficient(SlidingFrictionCoeff);
    //species->setSlidingFrictionCoefficientStatic(0.6);
    //
    species->setRollingStiffness(2./5.*species->getStiffness());
    species->setRollingDissipation(2./5.*species->getDissipation());
    species->setRollingFrictionCoefficient(RollingFrictionCoeff);
    //------------------------------------------------------------------------------ Tool species -------------------------------------------------------------
    //auto wallSpecies = deposition_problem.speciesHandler.copyAndAddObject(species);
    //auto wallParticleSpecies = deposition_problem.speciesHandler.getMixedObject(species,wallSpecies);
    //
    //wallParticleSpecies->setSlidingFrictionCoefficient(WSlidingFriCoeff);
    //wallParticleSpecies->setRollingFrictionCoefficient(WRollingFriCoeff);
    //
    //deposition_problem.setFixedParticleRadius(0.0019);
    //deposition_problem.setRoughBottomType(FLAT);
    //--------------------------------------------------------------------------------Parallel-------------------------------------------------------------------------
/*    if (activeParallel)
    {}*/

    //deposition_problem.setNumberOfDomains({6,1,1});
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