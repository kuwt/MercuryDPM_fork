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

class SpreadingOverRB: public Mercury3D
{
public:
    //
    bool deleteBoundary = true;
    bool insertParticlesOnly = false;
    bool RB = true; // true -> RB, false -> previous Sim
    bool fixBL = false;
    //
    int s = DPMBase::readRunNumberFromFile();
    const char *previousSimName = "SpreadingProcessRestart.";
    //char previousSimName = "SpreadingProcessRestart." + std::to_string(s) +".data";
    //
    const char *simName = "SpreadingProcessRB";//"SpreadingProcessInsertionRB";
    // ---------------------------------------PROCESS PARAMETERS ----------------------------------------------------------------
    // only insert RB: 0.1 s // at speed = 10 mm/s -> 1.1 s, 50mm/s -> 0.4 s , 100mm/s -> 0.3 s ;
    // from Previous Sim: at 10 speed = 2.0 s / 50 Speed = 0.8 s / 100 speed = 0.6 s
    Mdouble MaxSimTime = 1.1;
    //
    Mdouble TimeRemoveWall = 0.08;// RB -> 0.08, previous sim 1.0
    Mdouble TimeMoveTool = 0.1;// RB -> 0.1, previous sim 1.05
    //
    bool blade = true;//true; // false -> roller , true -> blade
    Mdouble tool_Speed = 0.01; //0.1 m/s - 100 mm/s / 0.05 m/s - 50 mm/s / 0.01 m/s - 10 mm/s
    //
    Mdouble Gap = 100e-6; // if from RB = 100e-6 , if from previous Sim = 200e-6
    Mdouble startPointBackDam = -Gap; // if from RB = -Gap (-layerThickness = 100e-6), if from previous sim = 0.0
    Mdouble startPointFrontDam = 100e-6; // in case of RB = 0.0xxx or (100e-6), if from previous Sim = layerThickness = 100e-6 (or Gap)
    //
    Mdouble particlesInsertionZ0 = 0.0;// if RB = 0.0, if previous Sim = 100e-6;
    // ---------------------------------------Material PARAMETERS - cohesion ------------------------------------------------------
    bool cohesiveON = false;
    Mdouble AdhStiffnessFactor = 0.5;
    Mdouble surfaceEnergy = 0.1e-3; //J/m2
    //
    //------------------------------------------------------------------------------------------------------------------------------
    //
    int particleGenerationSEED = 1;//23;//1; // first-> 1 / second -> 23
    //
    int numOfStudies1 = 7;//5;//2;//5;//2;//5;
    int numOfStudies2 = 6;//1;//6;//2;//1;//6;//5;
    //
    std::vector<double> SlidingFrictionCoeffVector = {0.5,0.4,0.3,0.25,0.2,0.1,0.05};
    std::vector<double> RollingFrictionCoeffVector = {0.4,0.3,0.2,0.1,0.05,0.005};
    //
    //std::vector<double> SlidingFrictionCoeffVector = {0.5,0.4,0.25,0.1,0.05};//{0.3,0.2};
    //std::vector<double> RollingFrictionCoeffVector = {0.3};//{0.4,0.3,0.2,0.1,0.05,0.005};
    // 4 cases -> diff seeds // 5,12 cases // intial cases
    //{0.4,0.25};//{0.5,0.4,0.25,0.1,0.05};//{0.2,0.3};//{0.5,0.4,0.25,0.1,0.05};
    //{0.2,0.1};//{0.3};//{0.4,0.3,0.2,0.1,0.05,0.005};//{0.4,0.2,0.1,0.05,0.005};
    //
    // new Sets:
    //std::vector<double> SlidingFrictionCoeffVector = {0.5,0.4,0.3,0.25,0.2,0.1,0.05};
    //std::vector<double> RollingFrictionCoeffVector = {0.4,0.3,0.2,0.1,0.05,0.005};
    //
    /*
    std::vector<double> SlidingFrictionCoeffVector = {10.0,2.0,1.0,0.5,0.4,0.25,0.1,0.05};
    std::vector<double> RollingFrictionCoeffVector = {10.0,0.4,0.2,0.15,0.1,0.075,0.05,0.005};*/
    //
    // ---------------------------------------Material PARAMETERS ---------------------------------------------------------------
    //
    Mdouble MatDensity = 4430;// density of Ti = 4.430g/cm3 - 4430kg/m3
    //
    Mdouble tc = 90e-7; //calculated as tc = 0.005-0.01 * tg : tg=sqrt(d50/g) => 0.005*sqrt(35e-6/9.81)=
    Mdouble timeStep = 0.02*tc;
    Mdouble restitutionCoeff = 0.4;//0.1;
    Mdouble minRaduis = 12.0e-6;///2.0;//MinRadius;//ALL PREVIOUS CASES THIS WAS USED -> 12.0e-6;   //15.0e-6/2.0;//13.183e-6/2.0;
    //------------------------------------------------------------------------------------------------------------------------------
    Mdouble gravityValue = -9.81; // g=9.81 m/s2 - 981 cm/s2
    //
    unsigned int saveCount = 10000; // 5000
    //
    //unsigned int numberCoresX = 12;
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------particles / PSD-------------------------------------------------------------------
    Mdouble meanRadius1 = 19.0e-6;//1.8577e-05;//0.0019; // in cm
    Mdouble MinRadius = 12.0e-6/2.0;//(13.183e-6/2.0); //in cm
    Mdouble MaxRadius = 79.0e-6/2.0;//(79.433e-6/2.0); //in cm
    // --------------------------------log-normal distribution mean and std, calculated from data points via lognfit matlab function-----------------------------------
    Mdouble LogmeanRadius = -10.8936;//-6.2884; //log(meanRadius1);
    Mdouble LogstdRadius = 0.3384; //log(stdRadius1);
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // -------------------------------------------------- initial parameters --------------------------------------------------------------------------------------------------------------
    //
    Mdouble scale = 1;//4e-1;//1e-2;
    //TOOL - initial configuration setup - unit: m:
    Mdouble SpreadingLength = 0.01*scale; // 10 mm -> 1 cm -> 0.01 m;
    Mdouble SpreadingWidth  = 0.001*scale; // 1 mm -> 0.1 cm -> 0.001 m
    Mdouble layerThickness = 100e-6;
    Mdouble ToolHight  = 0.003*scale; // 3 mm, 2 mm -> 0.2 cm -> 0.002 m
    Mdouble WallHight = ToolHight; //500e-6;
    Mdouble ToolThickness = 500e-6*scale; // 0.5 mm 250e-6;// m
    //Mdouble Gap = 100e-6; // 100 microns -> 100e-6 m ;
    //Mdouble TrackLength = 900e-1*(scale*scale);
    //Mdouble Angel = constants::pi/180.0*30;
    //Mdouble AngelDim = scale*4e-1;
    Mdouble rollerRadius = ToolThickness*2.0; // 0.001; // 1 mm
    //------------------------------------------------------------------------------ Domain ------------------------------------------------------------------------------
    Mdouble XMaxDomain = SpreadingLength+ToolThickness;//SpreadingLength/2.0+ToolThickness;//<- insertion domain+remove wall / -> full spreading domain SpreadingLength+ToolThickness;
    //(deposition_problem.SpreadingLength+deposition_problem.ToolThickness)*2.0+deposition_problem.ToolThickness;
    Mdouble YMinDomain = -ToolThickness;
    Mdouble YMaxDomain = SpreadingWidth;
    Mdouble ZMinDomain = -ToolThickness;
    Mdouble ZMaxDomain = ToolHight;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //TOOL
    // ----Balde:
    Vec3D minPoint_tool = Vec3D(0.0,0.0,Gap);
    Vec3D maxPoint_tool = Vec3D(ToolThickness,SpreadingWidth,ToolHight);
    // ----Roller:
    Vec3D rollerPosition = Vec3D(-ToolThickness,0.0,rollerRadius+Gap);
    //round blade:
    //Vec3D minPoint_toolround = Vec3D(0.0,0.0,Gap+ToolThickness);
    //Vec3D maxPoint_toolround = Vec3D(ToolThickness,SpreadingWidth,ToolHight);
    //Back Dam, small wall under tool
    Vec3D minPoint_back_dam = Vec3D(-ToolThickness,0.0,startPointBackDam);
    Vec3D maxPoint_back_dam = Vec3D(ToolThickness,SpreadingWidth,Gap);
    //front dam, present while inserting particles -> removed during spreading
    Mdouble damThickness = SpreadingLength/5.0+ToolThickness;// -> 2mm insertion domain
    Vec3D minPoint_dam = Vec3D(damThickness,0.0,startPointFrontDam);
    Vec3D maxPoint_dam = Vec3D(damThickness+ToolThickness,SpreadingWidth,WallHight);
    //
    //
    Vec3D minPoint_endWall = Vec3D(SpreadingLength+ToolThickness,0.0,startPointBackDam);
    Vec3D maxPoint_endWall = Vec3D(SpreadingLength+(2.0*ToolThickness),SpreadingWidth,Gap);
    //
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Mdouble minXparticles = ToolThickness;
    Mdouble maxXparticles = SpreadingLength/5.0+ToolThickness;
    Mdouble maxYparticles = SpreadingWidth;
    //
    //---------------------------------------------------------------- the volume to be added: ------------------------------------------------------------------------
    // correspond to spreading volume of desired layer 5 x 1 x 0.1 mm + initial insertion volume = 7 x 1 x 0.1 mm
    Mdouble addVolume0 = (SpreadingLength-(SpreadingLength/2.0))*SpreadingWidth*layerThickness + (SpreadingLength/5.0)*SpreadingWidth*layerThickness;
    //---------------------------------------------------------------------------------------------------
    /*Mdouble random_number(Mdouble min, Mdouble max)
    {
        // use thread_local to make this function thread safe
        thread_local static std::mt19937 mt{std::random_device{}()};
        thread_local static std::lognormal_distribution<Mdouble> dist;
        using pick = std::lognormal_distribution<Mdouble>::param_type;

        return dist(mt, pick(min, max));
    }*/
    //-------------------------------------------------- --------------------------------------------------
    // restart simulation:
    /*SpreadingLayerInsertAndSpread() {
        autoNumber();
        // restart all files
        auto s = std::to_string(DPMBase::readRunNumberFromFile() - 1);
        // restart from specific files
        //std::vector<int> set1 = {19,20,21,22,23,24,27,28,29,30,31,32,35,36,37,38,39,40,43,44,45,46,47,48,51,52,53,54,55,56,59,60,61,62,63,64};
        //std::vector<int> set2 = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,25,26,33,34,41,42,49,50,57,58};
        //std::vector<int> set3 = {11,12,13,14,15,16};
        //
        //std::vector<int> setSlow = {19,20,21,22,23,24,59,60,61,62,63,64};
        //
        //auto fileNumber = std::to_string(set3[DPMBase::readRunNumberFromFile()-2]);
        // start from restart
        setName("SpreadingProcessInsertion." + s);//fileNumber);
        readRestartFile();
        //setRestarted(false);
        setName("SpreadingProcessRestart." + s);//fileNumber);
        //species = dynamic_cast<LinearViscoelasticFrictionJKRAdhesiveSpecies *>(speciesHandler.getObject(0));
    };*/
    //----------------------------------------------------------------------------------------------------
    unsigned nParticlesLayerRB;
    Mdouble getInfo(const BaseParticle& P) const override {
        if (P.getIndex()<nParticlesLayerRB)
            return 0;
        else
            return 1;
    }
   //----------------------------------------------------------------------------------------------------
    void setupInitialConditions() override {
        //Bottom Wall
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0,0.0,startPointBackDam));
        wallHandler.copyAndAddObject(w0);
        //TOOL:
        //------Blade:
        IntersectionOfWalls w1;
        w1.setSpecies(speciesHandler.getObject(0)); //using different species for the Tool: w0.setSpecies(speciesHandler.getObject(1));
        w1.addObject(Vec3D(1.0,0.0,0.0),minPoint_tool);
        w1.addObject(Vec3D(0.0,1.0,0.0),minPoint_tool);
        w1.addObject(Vec3D(0.0,0.0,1.0),minPoint_tool);
        w1.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_tool);
        w1.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_tool);
        w1.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_tool); //w1.addObject(Vec3D(cos(Angel),0.0,sin(Angel)),Vec3D(AngelDim,0.0,Gap)); //in next line
        wallHandler.copyAndAddObject(w1);
        if(!blade)
        {
            //------Roller: w1
            AxisymmetricIntersectionOfWalls roller;
            roller.setSpecies(speciesHandler.getObject(0));
            roller.setPosition(rollerPosition);
            roller.setAxis({0, 1, 0});//setAxis
            roller.addObject({-1, 0, 0}, {rollerRadius, 0, 0});
            roller.addObject({0, -1, 0}, {0, 0, 0});
            //roller.setAngularVelocity({0,-tool_Speed/rollerRadius,0});
            wallHandler.copyAndAddObject(roller);
        }
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
        //Front wall
        IntersectionOfWalls w3;
        w3.setSpecies(speciesHandler.getObject(0));
        w3.addObject(Vec3D(1.0,0.0,0.0),minPoint_dam);
        w3.addObject(Vec3D(0.0,1.0,0.0),minPoint_dam);
        w3.addObject(Vec3D(0.0,0.0,1.0),minPoint_dam);
        w3.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_dam);
        w3.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_dam);
        w3.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_dam);
        wallHandler.copyAndAddObject(w3);
        //
        // Dam at end of RB
       IntersectionOfWalls endWall;
       endWall.setSpecies(speciesHandler.getObject(0));
       endWall.addObject(Vec3D(1.0,0.0,0.0),minPoint_endWall);
       endWall.addObject(Vec3D(0.0,1.0,0.0),minPoint_endWall);
       endWall.addObject(Vec3D(0.0,0.0,1.0),minPoint_endWall);
       endWall.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_endWall);
       endWall.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_endWall);
       endWall.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_endWall);
       wallHandler.copyAndAddObject(endWall);
        //
        PeriodicBoundary b0;
        b0.set(Vec3D(0, 1, 0), 0.0, SpreadingWidth);
        boundaryHandler.copyAndAddObject(b0);
        //----------------------------------------------------------------------- Particles insertion ------------------------------------------------------------------
       // read data files
       if (RB)
       {
           readDataFile("roughBaseFullPSD.data", 14);
           if(fixBL)
           {
               for (int i = 0; i < particleHandler.getNumberOfObjects(); ++i) {
                   Vec3D pP = particleHandler.getObject(i)->getPosition();
                   if (pP.Z <= 0.0)
                       particleHandler.getObject(i)->fixParticle();
               }
               logger(INFO,"Base particles fixed", particleHandler.getSize());
           } else{
               logger(INFO,"Base Particles NOT fixed", particleHandler.getSize());
           }
           nParticlesLayerRB = particleHandler.getSize();
       }
       else
       {
           readDataFile(previousSimName + std::to_string(s) + ".data", 14);
           nParticlesLayerRB = particleHandler.getSize();
           logger(INFO,"Read from previous Sim", particleHandler.getSize());
       }

       //
        //set up random number generator
        //std::random_device rd;
        //std::mt19937 gen(rd());
        std::mt19937 gen;
        gen.seed(particleGenerationSEED);
        std::lognormal_distribution<> d(LogmeanRadius, LogstdRadius);
        //
        //add particles until the volume to be added is zero
        //logger(INFO,"Adding particles ...");
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(meanRadius1);
        Mdouble fillHeight0 = particlesInsertionZ0;
        while (addVolume0>0) {
            Mdouble x = random.getRandomNumber(minXparticles, maxXparticles);
            Mdouble y = random.getRandomNumber(0.0, maxYparticles); //(ToolThickness, (ToolWidth-ToolThickness));
            Mdouble z = random.getRandomNumber(particlesInsertionZ0, fillHeight0);
            p0.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p0)) {
                particleHandler.copyAndAddObject(p0);
                addVolume0 -= p0.getVolume();
                do {
                    p0.setRadius(d(gen));
                } while (p0.getRadius()< MinRadius || p0.getRadius() > MaxRadius);
            } else {
                fillHeight0 += 1e-4*meanRadius1; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        //
        logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());
        //
        // delete particles outside domain of interest:
        if(deleteBoundary)
        {
            DeletionBoundary deletionBoundary1;
            deletionBoundary1.set(Vec3D(1, 0, 0), SpreadingLength); // - (SpreadingLength / 5.0)) + 4.0*ToolThickness);//(SpreadingLength-0.005)+ToolThickness);
            boundaryHandler.copyAndAddObject(deletionBoundary1);
            DeletionBoundary deletionBoundary2;
            deletionBoundary2.set(Vec3D(-1, 0, 0), ToolThickness);//(SpreadingLength-0.005)+ToolThickness);
            boundaryHandler.copyAndAddObject(deletionBoundary2);
        }
        //
    }
    //
    //----------------------------------------------------------------------------------------------------------------------------------------------
    void actionsAfterTimeStep() override {
        //
        // Remove dams
        static bool wallRemoved = false;
        if (wallRemoved == false && getTime() >= TimeRemoveWall) {
            //logger(INFO,"walls removed");
            if (blade) {
                wallHandler.removeObject(3);
            } else {
                wallHandler.removeObject(4); // front dam -> 3 / when using a roller wallHandler.removeObject(4);
                wallHandler.removeObject(1); // tool,blade -> remove when using a roller
            }
            wallRemoved = true;
        }
        // Start Spreading
        if (!insertParticlesOnly)
        {
            static bool moveTool = false;
            if (!moveTool && getTime() >= TimeMoveTool)
            {
                //logger(INFO,"Tool is Moving");
                if (blade)
                {
                    //blade:
                    wallHandler.getObject(1)->setVelocity(Vec3D(tool_Speed, 0, 0));
                }
                else
                {
                    //roller:
                    wallHandler.getLastObject()->setVelocity(Vec3D(tool_Speed, 0, 0));
                    wallHandler.getLastObject()->setAngularVelocity({0, -tool_Speed / rollerRadius,0});//getObject(1)->setAngularVelocity({0,-tool_Speed/rollerRadius,0});
                }
                moveTool = true;
            }
        }
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
    SpreadingOverRB deposition_problem;
    //
    //------------------------------------------------------------------------------Parametric friction:------------------------------------------------------------------------------
    deposition_problem.autoNumber();
    //
    std::vector<int> studyNum=deposition_problem.get2DParametersFromRunNumber(deposition_problem.numOfStudies1,deposition_problem.numOfStudies2);
    //
    Mdouble SlidingFrictionCoeff = deposition_problem.SlidingFrictionCoeffVector[studyNum[1]-1];
    Mdouble RollingFrictionCoeff = deposition_problem.RollingFrictionCoeffVector[studyNum[2]-1];
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
    deposition_problem.setName(deposition_problem.simName);
    deposition_problem.setSystemDimensions(3);
    deposition_problem.setGravity(Vec3D(0.0,0.0,deposition_problem.gravityValue));
    //deposition_problem.setXMin(0.0);
    deposition_problem.setYMin(deposition_problem.YMinDomain);
    deposition_problem.setZMin(deposition_problem.ZMinDomain);
    deposition_problem.setXMax(deposition_problem.XMaxDomain);
    deposition_problem.setYMax(deposition_problem.YMaxDomain);
    deposition_problem.setZMax(deposition_problem.ZMaxDomain);
    deposition_problem.setTimeMax(deposition_problem.MaxSimTime);
    //deposition_problem.setHGridMaxLevels(2);
    //
    //auto species = deposition_problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    auto species = deposition_problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionJKRAdhesiveSpecies());
    //
    species->setDensity(deposition_problem.MatDensity);
    //
    // Here we assign the mass of the smallest particle to calculate the stiffiness and disppation as we want to get ------
    Mdouble mass = species->getMassFromRadius(deposition_problem.minRaduis);
    // ----- set stiffness and dissipation from tc,COR,mass
    species->setCollisionTimeAndRestitutionCoefficient(deposition_problem.tc,deposition_problem.restitutionCoeff,mass);
    //------------------------------------------------------------------------------ Cohesion Parameters ------------------------------------
    if (deposition_problem.cohesiveON)
    {
        Mdouble adhesionStiffness = deposition_problem.AdhStiffnessFactor * species->getStiffness();
        species->setAdhesionStiffness(adhesionStiffness);
        species->setSurfaceEnergy(deposition_problem.surfaceEnergy);
    }
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
/*    auto toolSpecies = deposition_problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionJKRAdhesiveSpecies());
    toolSpecies->setDensity(deposition_problem.MatDensity);
    Mdouble massPT = toolSpecies->getMassFromRadius(deposition_problem.minRaduis);
    // ----- set stiffness and dissipation from tc,COR,mass
    toolSpecies->setCollisionTimeAndRestitutionCoefficient(deposition_problem.tc,deposition_problem.restitutionCoeff,massPT);
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
    //*/
    //--------------------------------------------------------------------------------Parallel-------------------------------------------------------------------------
    //deposition_problem.setNumberOfDomains({deposition_problem.numberCoresX,1,1});
    //------------------------------------------------------------------------------Data Output--------------------------------------------------------------------------------------------
    deposition_problem.setSaveCount(deposition_problem.saveCount);//50
    deposition_problem.dataFile.setFileType(FileType::ONE_FILE);
    //deposition_problem.fStatFile.setFileType(FileType::NO_FILE);
    deposition_problem.setParticlesWriteVTK(true);
    deposition_problem.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    //
    deposition_problem.setTimeStep(deposition_problem.timeStep);
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