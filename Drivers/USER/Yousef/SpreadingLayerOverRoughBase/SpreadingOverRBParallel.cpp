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

class SpreadingOverRBParallel: public Mercury3D
{
public:
    // -------------------------------------------------- initial parameters -------------------------------------------------------
    Mdouble TimeRemoveWall = 0.08;
    Mdouble TimeMoveTool = 0.1;
    //
    Mdouble tool_Speed = 0.05;//0.1; //first run -> 0.05 - 50 mm/s;//0.01;// m/s
    //
    Mdouble scale = 1.0;//3e-1;
    //TOOL - initial configuration setup - unit: m:
    Mdouble SpreadingLength = 0.01*scale; // 10 mm -> 1 cm -> 0.01 m;
    Mdouble SpreadingWidth  = 0.001*scale; // 1 mm -> 0.1 cm -> 0.001 m
    Mdouble ToolHight  = 0.003*scale; // 3 mm, 2 mm -> 0.2 cm -> 0.002 m
    Mdouble WallHight = ToolHight;//500e-6;
    Mdouble ToolThickness = 500e-6*scale;// 0.5 mm 250e-6;// m
    Mdouble Gap = 100e-6;//100e-6; //100 microns -> 100e-6 m ;
    //Mdouble TrackLength = 900e-1*(scale*scale);
    //Mdouble Angel = constants::pi/180.0*30;
    //Mdouble AngelDim = scale*4e-1;
    //------------------------------------------------------------------------------ Domain ------------------------------------------------------------------------------
    Mdouble XMaxDomain = SpreadingLength+ToolThickness;
    //(deposition_problem.SpreadingLength+deposition_problem.ToolThickness)*2.0+deposition_problem.ToolThickness;
    Mdouble YMaxDomain = SpreadingWidth;
    Mdouble ZMinDomain = -ToolThickness;
    Mdouble ZMaxDomain = ToolHight;
    //----------------------------------------------------------------------------------------------------
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
    Vec3D minPoint_back_dam = Vec3D(0.0,0.0,-200e-6);
    Vec3D maxPoint_back_dam = Vec3D(ToolThickness,SpreadingWidth,Gap);
    //
    Mdouble damThickness = SpreadingLength/5.0+ToolThickness;// -> 2mm insertion domain
    Vec3D minPoint_dam = Vec3D(damThickness,0.0,50e-6);
    Vec3D maxPoint_dam = Vec3D(damThickness+ToolThickness,SpreadingWidth,WallHight);
    //-------------------------------------------------------------------------------particles / PSD-------------------------------------------------------------------
    int particleGenerationSEED = 1;
    Mdouble meanRadius1 = 19.0e-6;//1.8577e-05;//0.0019; // in cm
    Mdouble MinRadius = 12.0e-6/2.0;//(13.183e-6/2.0); //in cm
    Mdouble MaxRadius = 79.0e-6/2.0;//(79.433e-6/2.0); //in cm
    //Mdouble mode_log_normal = 34.674e-4/2.0;
    //Mdouble stdRadius1 = 1.4027; //0.2 * meanRadius1;
    // --------------------------------log-normal distribution mean and std, calculated from data points via lognfit matlab function-----------------------------------
    Mdouble LogmeanRadius = -10.8936;//-6.2884; //log(meanRadius1);
    Mdouble LogstdRadius = 0.3384; //log(stdRadius1);
    //---------------------------------------------------------------- the volume to be added: ------------------------------------------------------------------------
    Mdouble addVolume0 = (SpreadingLength-(SpreadingLength/2.0))*SpreadingWidth*Gap + (SpreadingLength/5.0)*SpreadingWidth*Gap;// correspond to spreading volume of desired layer 5 x 1 x 0.1 mm
    //
    Mdouble minXparticles = ToolThickness;
    Mdouble maxXparticles = SpreadingLength/5.0+ToolThickness;
    Mdouble maxYparticles = SpreadingWidth;
    //------------------------------------------------------------------------------ Parameters ------------------------------------------------------------------------------
    Mdouble gravityValue = -9.81; // g=9.81 m/s2 - 981 cm/s2
    //
    const char *simName = "SpreadingOverRB";//Insertion
    //
    Mdouble MaxSimTime = 0.4;// at 100 mm/s -> 0.3; // at  50mm/s -> 0.4 // at 10 mm/s -> 1.0;
    //
    Mdouble MatDensity = 4430;// density of Ti = 4.430g/cm3 - 4430kg/m3
    //
    Mdouble tc = 90e-7; //calculated as tc = 0.005-0.01 * tg : tg=sqrt(d50/g) => 0.005*sqrt(35e-6/9.81)
    //
    Mdouble restitutionCoeff = 0.4;//0.1;
    //
    Mdouble timeStep = 0.02*tc;
    //
    Mdouble minRaduis = 12.0e-6;//MinRadius
    //
    unsigned int saveCount = 10000; // 5000
    //
    Mdouble AdhStiffnessFactor = 0.5;
    Mdouble surfaceEnergy = 0.1e-3; //J/m2
    //----------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------- Particles insertion ------------------------------------------------------------------
    int numOfStudies = 5;
    std::vector<double> SlidingFrictionCoeffVector = {0.5,0.4,0.25,0.1,0.05};
    std::vector<double> RollingFrictionCoeffVector = {0.4,0.2,0.1,0.05,0.005};
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // To make two species fixed and free particles
/*    unsigned nParticlesLayerRB;

    Mdouble getInfo(const BaseParticle& P) const override {
        if (P.getIndex()<nParticlesLayerRB)
            return 0;
        else
            return 1;
    }*/
//----------------------------------------------------------------------------------------------------
    void setupInitialConditions() override {

        //-------------------------------------------------------------------------------particles / PSD-------------------------------------------------------------------
        //------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //Bottom Wall
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0, 0.0, -200e-6));
        wallHandler.copyAndAddObject(w0);
        //TOOL:
        //------Blade:
        IntersectionOfWalls w1;
        w1.setSpecies(speciesHandler.getObject(
                1)); //using different species for the Tool: w0.setSpecies(speciesHandler.getObject(1));
        w1.addObject(Vec3D(1.0, 0.0, 0.0), minPoint_tool);
        w1.addObject(Vec3D(0.0, 1.0, 0.0), minPoint_tool);
        w1.addObject(Vec3D(0.0, 0.0, 1.0), minPoint_tool);
        w1.addObject(Vec3D(-1.0, 0.0, 0.0), maxPoint_tool);
        w1.addObject(Vec3D(0.0, -1.0, 0.0), maxPoint_tool);
        w1.addObject(Vec3D(0.0, 0.0, -1.0),
                     maxPoint_tool); //w1.addObject(Vec3D(cos(Angel),0.0,sin(Angel)),Vec3D(AngelDim,0.0,Gap)); //in next line
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
        w2.addObject(Vec3D(1.0, 0.0, 0.0), minPoint_back_dam);
        w2.addObject(Vec3D(0.0, 1.0, 0.0), minPoint_back_dam);
        w2.addObject(Vec3D(0.0, 0.0, 1.0), minPoint_back_dam);
        w2.addObject(Vec3D(-1.0, 0.0, 0.0), maxPoint_back_dam);
        w2.addObject(Vec3D(0.0, -1.0, 0.0), maxPoint_back_dam);
        w2.addObject(Vec3D(0.0, 0.0, -1.0), maxPoint_back_dam);
        wallHandler.copyAndAddObject(w2);
        //
        IntersectionOfWalls w3;
        w3.setSpecies(speciesHandler.getObject(0));
        w3.addObject(Vec3D(1.0, 0.0, 0.0), minPoint_dam);
        w3.addObject(Vec3D(0.0, 1.0, 0.0), minPoint_dam);
        w3.addObject(Vec3D(0.0, 0.0, 1.0), minPoint_dam);
        w3.addObject(Vec3D(-1.0, 0.0, 0.0), maxPoint_dam);
        w3.addObject(Vec3D(0.0, -1.0, 0.0), maxPoint_dam);
        w3.addObject(Vec3D(0.0, 0.0, -1.0), maxPoint_dam);
        wallHandler.copyAndAddObject(w3);
        //
        PeriodicBoundary b0;
        b0.set(Vec3D(0, 1, 0), 0.0, SpreadingWidth);
        boundaryHandler.copyAndAddObject(b0);
        //------------------------------------------------------------------------
        //autoNumber();
        // read data files
        readDataFile("roughBaseFullPSD.data", 14);
        //if (readDataFile("roughBase.data", 14))
        //{
            for (int i = 0; i < particleHandler.getNumberOfRealObjectsLocal(); ++i)
            {
                Vec3D pP = particleHandler.getObject(i)->getPosition();
                if (pP.Z < 0.0)
                    particleHandler.getObject(i)->fixParticle();
            }
        //}
        //
/*        readDataFile("roughBase2.data",14);
        if (readDataFile("roughBase2.data",14))
        {
            for (int j = 0; j < particleHandler.getNumberOfRealObjects(); ++j)
            {
                Vec3D pP = particleHandler.getObject(j)->getPosition();
                //Mdouble pRadius = particleHandler.getObject(i)->getRadius();
                particleHandler.getObject(j)->setPosition(Vec3D(pP.X + 0.0105, pP.Y, pP.Z));
                if (pP.Z < 0.0)
                    particleHandler.getObject(j)->fixParticle();
            }
        }*/
        //nParticlesLayerRB = particleHandler.getSize();
        //------------------------------------------------------------------------
        std::mt19937 gen;
        gen.seed(1);
        std::lognormal_distribution<> d(LogmeanRadius, LogstdRadius);
        //
        //add particles until the volume to be added is zero
        //logger(INFO,"Adding particles ...");
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(meanRadius1);
        Mdouble fillHeight0 = 0.0;
        while (addVolume0>0) {
            Mdouble x = random.getRandomNumber(minXparticles, maxXparticles);
            Mdouble y = random.getRandomNumber(0.0, maxYparticles); //(ToolThickness, (ToolWidth-ToolThickness));
            Mdouble z = random.getRandomNumber(0.0, fillHeight0);
            p0.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p0)) {
                particleHandler.copyAndAddObject(p0);
                addVolume0 -= p0.getVolume();
                do {
                    //p.setRadius(random_number(LogmeanRadius,LogstdRadius));
                    p0.setRadius(d(gen));
                } while (p0.getRadius()< MinRadius || p0.getRadius() > MaxRadius);
            } else {
                fillHeight0 += 1e-4*meanRadius1; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        //----------------------------------------------------------------------------------------------------------------------------------------------
        // delete particles outside domain of interest:
        // 1. delete particles that falls outside the doamin at the end of the layer
/*        DeletionBoundary deletionBoundary1;
        deletionBoundary1.set(Vec3D(0, 0, -1), SpreadingLength+ToolThickness);//(SpreadingLength-0.002)+ToolThickness);
        boundaryHandler.copyAndAddObject(deletionBoundary1);
        // 2. delete particles that go outside domain before spreading
        DeletionBoundary deletionBoundary2;
        deletionBoundary2.set(Vec3D(-1, 0, 0), ToolThickness);//(SpreadingLength-0.005)+ToolThickness);
        boundaryHandler.copyAndAddObject(deletionBoundary2);*/
        //
    }
    //
    //----------------------------------------------------------------------------------------------------------------------------------------------
    void actionsAfterTimeStep() override {
        //
        static bool wallRemoved = false;
        if (wallRemoved==false && getTime()>=TimeRemoveWall) {
            wallHandler.removeObject(3); // front dam -> 3
            wallRemoved = true;
        }
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
    SpreadingOverRBParallel deposition_problem;
    //
    //------------------------------------------------------------------------------Parametric friction:---------------------------------------------------------------------
    deposition_problem.autoNumber();
    std::vector<int> studyNum=deposition_problem.get2DParametersFromRunNumber(deposition_problem.numOfStudies,deposition_problem.numOfStudies);
    //
    Mdouble SlidingFrictionCoeff = deposition_problem.SlidingFrictionCoeffVector[studyNum[1]-1];
    Mdouble RollingFrictionCoeff = deposition_problem.RollingFrictionCoeffVector[studyNum[2]-1];
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    //
    deposition_problem.setName(deposition_problem.simName);
    deposition_problem.setSystemDimensions(3);
    deposition_problem.setGravity(Vec3D(0.0,0.0,deposition_problem.gravityValue));
    //deposition_problem.setXMin(0.0);
    //deposition_problem.setYMin(0.0);
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
    Mdouble adhesionStiffness = deposition_problem.AdhStiffnessFactor*species->getStiffness();
    species->setAdhesionStiffness(adhesionStiffness);
    species->setSurfaceEnergy(deposition_problem.surfaceEnergy);
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
    auto toolSpecies = deposition_problem.speciesHandler.copyAndAddObject(species);
    auto toolParticleSpecies = deposition_problem.speciesHandler.getMixedObject(species,toolSpecies);
    //
/*    auto toolSpecies = deposition_problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionJKRAdhesiveSpecies());
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
    toolParticleSpecies->setRollingFrictionCoefficient(species->getRollingFrictionCoefficient()/4.0);*/
    //
    //--------------------------------------------------------------------------------Parallel-------------------------------------------------------------------------
    deposition_problem.setNumberOfDomains({12,1,1});
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