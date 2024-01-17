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

class makeRoughBase: public Mercury3D
{
public:
    // -------------------------------------------------- initial parameters -------------------------------------------------------
    Mdouble TimeRemoveWall = 0.05;
    Mdouble TimeFixParticles = 0.048;
    //
    //TOOL - initial configuration setup - unit: m:
    Mdouble SpreadingLength = 0.01; // 10 mm -> 1 cm -> 0.01 m;
    Mdouble SpreadingWidth  = 0.001; // 1 mm -> 0.1 cm -> 0.001 m
    Mdouble ToolHight  = 0.003; // 3 mm, 2 mm -> 0.2 cm -> 0.002 m
    Mdouble WallHight = ToolHight;//500e-6;
    Mdouble ToolThickness = 500e-6;// 0.5 mm 250e-6;// m
    Mdouble Gap = 100e-6;//100e-6; //100 microns -> 100e-6 m ;
    //
    Mdouble particleLayerThickness = 100e-6;
    //-------------------------------------------------- --------------------------------------------------
    //-------------------------------------------------- --------------------------------------------------
/*    unsigned nParticlesLayerRB;
    unsigned nParticlesLayer1;

    Mdouble getInfo(const BaseParticle& P) const override {
        if (P.getIndex()<nParticlesLayerRB)
            return 10;
    }*/
    //-------------------------------------------------- --------------------------------------------------
    void setupInitialConditions() override {
        //TOOL
        // ----Balde:
        Vec3D minPoint_tool = Vec3D(0.0,0.0,Gap);
        Vec3D maxPoint_tool = Vec3D(ToolThickness,SpreadingWidth,ToolHight);
        //Back Dam
        Vec3D minPoint_back_dam = Vec3D(0.0,0.0,-particleLayerThickness);
        Vec3D maxPoint_back_dam = Vec3D(ToolThickness,SpreadingWidth,Gap);
        //
        //front Dam1
        Mdouble damThicknessF = (SpreadingLength + ToolThickness);//(5.0e-1*scale);
        Vec3D minPoint_damF = Vec3D(damThicknessF, 0.0, -particleLayerThickness);//ToolLength*2.0
        Vec3D maxPoint_damF = Vec3D(damThicknessF + ToolThickness, SpreadingWidth, 0.0);
        //front Dam2
        Vec3D minPoint_damF2 = Vec3D(damThicknessF,0.0,0.0);
        Vec3D maxPoint_damF2 = Vec3D(damThicknessF+ToolThickness,SpreadingWidth,ToolHight);
        //-------------------------------------------------------------------------------particles / PSD-------------------------------------------------------------------
        Mdouble meanRadius1 = 19.0e-6;//1.8577e-05;//0.0019; // in cm
        Mdouble MinRadius = 12.0e-6/2.0;;// first one -> 22.0e-6/2.0; //in cm
        Mdouble MaxRadius = 79.0e-6/2.0;//first one -> 55.0e-6/2.0;////in cm
        //Mdouble mode_log_normal = 34.674e-4/2.0;
        //Mdouble stdRadius1 = 1.4027; //0.2 * meanRadius1;
        // --------------------------------log-normal distribution mean and std, calculated from data points via lognfit matlab function-----------------------------------
        Mdouble LogmeanRadius = -10.8936;//-6.2884; //log(meanRadius1);
        Mdouble LogstdRadius = 0.3384; //log(stdRadius1);
        //---------------------------------------------------------------- the volume to be added: ------------------------------------------------------------------------
        Mdouble addVolumeRB = SpreadingLength * SpreadingWidth * particleLayerThickness;//(Gap*2.0);
        //------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //Bottom Wall
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0,0.0,-particleLayerThickness));
        wallHandler.copyAndAddObject(w0);
        //Back Dam
        IntersectionOfWalls w1;
        w1.setSpecies(speciesHandler.getObject(0));
        w1.addObject(Vec3D(1.0,0.0,0.0),minPoint_back_dam);
        w1.addObject(Vec3D(0.0,1.0,0.0),minPoint_back_dam);
        w1.addObject(Vec3D(0.0,0.0,1.0),minPoint_back_dam);
        w1.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_back_dam);
        w1.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_back_dam);
        w1.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_back_dam);
        wallHandler.copyAndAddObject(w1);
        //Front Dam1
        IntersectionOfWalls w2;
        w2.setSpecies(speciesHandler.getObject(0));
        w2.addObject(Vec3D(1.0, 0.0, 0.0), minPoint_damF);
        w2.addObject(Vec3D(0.0, 1.0, 0.0), minPoint_damF);
        w2.addObject(Vec3D(0.0, 0.0, 1.0), minPoint_damF);
        w2.addObject(Vec3D(-1.0, 0.0, 0.0), maxPoint_damF);
        w2.addObject(Vec3D(0.0, -1.0, 0.0), maxPoint_damF);
        w2.addObject(Vec3D(0.0, 0.0, -1.0), maxPoint_damF);
        wallHandler.copyAndAddObject(w2);
        // Front - Tool
        IntersectionOfWalls w3;
        w3.setSpecies(speciesHandler.getObject(0));
        w3.addObject(Vec3D(1.0,0.0,0.0),minPoint_tool);
        w3.addObject(Vec3D(0.0,1.0,0.0),minPoint_tool);
        w3.addObject(Vec3D(0.0,0.0,1.0),minPoint_tool);
        w3.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_tool);
        w3.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_tool);
        w3.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_tool);
        wallHandler.copyAndAddObject(w3);
        //front Dam 2
        IntersectionOfWalls w4;
        w4.setSpecies(speciesHandler.getObject(0));
        w4.addObject(Vec3D(1.0,0.0,0.0),minPoint_damF2);
        w4.addObject(Vec3D(0.0,1.0,0.0),minPoint_damF2);
        w4.addObject(Vec3D(0.0,0.0,1.0),minPoint_damF2);
        w4.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_damF2);
        w4.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_damF2);
        w4.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_damF2);
        wallHandler.copyAndAddObject(w4);
        //
        PeriodicBoundary b0;
        b0.set(Vec3D(0, 1, 0), 0.0, SpreadingWidth);
        boundaryHandler.copyAndAddObject(b0);
        //----------------------------------------------------------------------- Particles insertion ------------------------------------------------------------------
        //set up random number generator
        //std::random_device rd;
        //std::mt19937 gen(rd());
        std::mt19937 gen;
        gen.seed(1);
        std::lognormal_distribution<> d(LogmeanRadius, LogstdRadius);
        //
        SphericalParticle pRB;
        pRB.setSpecies(speciesHandler.getObject(0)); // RBSpecies -> no friction
        pRB.setRadius(meanRadius1);
        Mdouble fillHeightRB = 0.0;
        while (addVolumeRB > 0) {
            Mdouble x = random.getRandomNumber(ToolThickness, SpreadingLength+ToolThickness);
            Mdouble y = random.getRandomNumber(0.0, SpreadingWidth); //(ToolThickness, (ToolWidth-ToolThickness));
            Mdouble z = random.getRandomNumber(-particleLayerThickness, fillHeightRB);
            pRB.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(pRB)) {
                particleHandler.copyAndAddObject(pRB);
                addVolumeRB -= pRB.getVolume();
                do {
                    //p.setRadius(random_number(LogmeanRadius,LogstdRadius));
                    pRB.setRadius(d(gen));
                } while (pRB.getRadius() < MinRadius || pRB.getRadius() > MaxRadius);
                if (particleHandler.getNumberOfObjects()%100==0) std::cout << '.' << std::flush;
                if (particleHandler.getNumberOfObjects()%1000==0) std::cout << ' ';
                if (particleHandler.getNumberOfObjects()%10000==0) std::cout << addVolumeRB << '\n';
            } else {
                fillHeightRB += 0.8e-4 * meanRadius1; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());
        //nParticlesLayerRB = particleHandler.getSize();
        //
    }
    //
    //----------------------------------------------------------------------------------------------------------------------------------------------
    // Remove dams
    void actionsAfterTimeStep() override {
        //
        static bool wallRemoved = false;
        if (wallRemoved==false && getTime()>=TimeRemoveWall) {
            wallHandler.removeObject(4);
            wallHandler.removeObject(3); //
            wallRemoved = true;
        }
        // Make rough base:
        if (getTime()>= TimeFixParticles) {
            for (int i = 0; i < particleHandler.getNumberOfObjects(); ++i) {
                Vec3D particlePosition = particleHandler.getObject(i)->getPosition();
                Mdouble kinEnergy = particleHandler.getObject(i)->getKineticEnergy();
                if (particlePosition.Z > 0.0)
                {
                    particleHandler.removeObject(i);
                }
                else if (particlePosition.Z < 0.0 && kinEnergy < 2e-11) {
                    particleHandler.getObject(i)->fixParticle();
                    //logger(INFO, "particles fixed behind tool");
                }
            }
        }
    }
};
// --------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc UNUSED, char *argv[] UNUSED)
{
    Time time;
    time.tic();

    //logger(INFO,"Simple box for creating particles");
    makeRoughBase deposition_problem;
    //
    //------------------------------------------------------------------------------ Domain ------------------------------------------------------------------------------
    Mdouble gravityValue = -9.81; // g=9.81 m/s2 - 981 cm/s2
    Mdouble XMaxDomain = deposition_problem.SpreadingLength+2*deposition_problem.ToolThickness;
    //(deposition_problem.SpreadingLength+deposition_problem.ToolThickness)*2.0+deposition_problem.ToolThickness;
    Mdouble YMaxDomain = deposition_problem.SpreadingWidth;
    Mdouble ZMinDomain = -deposition_problem.particleLayerThickness;
    Mdouble ZMaxDomain = deposition_problem.ToolHight;
    //------------------------------------------------------------------------------ Parameters ------------------------------------------------------------------------------
    const char *simName = "makeRoughBaseFullPSD";
    //
    Mdouble MaxSimTime = 0.05;//0.2; // at speed = 50mm/s -> 0.5, but still too long make it 0.3;
    //
    Mdouble MatDensity = 4430;// density of Ti = 4.430g/cm3 - 4430kg/m3
    //
    Mdouble tc = 90e-7; //calculated as tc = 0.005-0.01 * tg : tg=sqrt(d50/g) => 0.005*sqrt(35e-6/9.81)
    //
    Mdouble restitutionCoeff = 0.4;//0.1;
    //
    Mdouble timeStep = 0.02*tc;
    //
    Mdouble minRaduis = 12.0e-6/2.0;//15.0e-6/2.0;//13.183e-6/2.0;
    //
    unsigned int saveCount = 10000; // 5000
    //
    Mdouble AdhStiffnessFactor = 0.5;
    Mdouble surfaceEnergy = 0.1e-3; //J/m2
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
    // RB properties:
    // RB species
    auto speciesRB = deposition_problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionJKRAdhesiveSpecies());
    speciesRB->setDensity(MatDensity);
    Mdouble mass = speciesRB->getMassFromRadius(minRaduis);
    speciesRB->setCollisionTimeAndRestitutionCoefficient(tc,restitutionCoeff,mass);
/*    speciesRB->setSlidingStiffness(2./7.*speciesRB->getStiffness());
    speciesRB->setSlidingDissipation(2./7.*speciesRB->getDissipation());
    speciesRB->setSlidingFrictionCoefficient(0.0);
    speciesRB->setRollingStiffness(2./5.*speciesRB->getStiffness());
    speciesRB->setRollingDissipation(2./5.*speciesRB->getDissipation());
    speciesRB->setRollingFrictionCoefficient(0.0);*/
    //
    //------------------------------------------------------------------------------ Cohesion Parameters ------------------------------------
/*    Mdouble adhesionStiffness = AdhStiffnessFactor*species->getStiffness();
    species->setAdhesionStiffness(adhesionStiffness);
    species->setSurfaceEnergy(surfaceEnergy);*/
    //------------------------------------------------------------------------------Sliding and rolling stiffness and dissipation -------------------------------------
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