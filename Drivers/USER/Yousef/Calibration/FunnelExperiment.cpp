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
#include "Walls/AxisymmetricIntersectionOfWalls.h"
//#include "Species/LinearViscoelasticFrictionSpecies.h"
//#include "Species/LinearViscoelasticSpecies.h"
//#include "Boundaries/CubeInsertionBoundary.h"
//#include "Chute.h"
//#include <ChuteBottom.h>
//#include "stdio.h"

class FunnelExperiment: public Mercury3D
{
public:
    //
    //bool deleteBoundary = false;
    //bool insertParticlesOnly = true;
    // ---------------------------------------Material PARAMETERS ---------------------------------------------------------------
    Mdouble MaxSimTime = 1.0;
    //
    Mdouble MatDensity = 4430;// density of Ti = 4.430g/cm3 - 4430kg/m3
    //
    Mdouble tc = 90e-7; //calculated as tc = 0.005-0.01 * tg : tg=sqrt(d50/g) => 0.005*sqrt(35e-6/9.81)=
    Mdouble timeStep = 0.02*tc;
    Mdouble restitutionCoeff = 0.4;//0.1;
    Mdouble minRaduis = 12.0e-6;///2.0;//MinRadius;//ALL PREVIOUS CASES THIS WAS USED -> 12.0e-6;   //15.0e-6/2.0;//13.183e-6/2.0;
    // ---------------------------------------Material PARAMETERS - cohesion ------------------------------------------------------
    bool cohesiveON = true;
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
    //------------------------------------------------------------------------------------------------------------------------------
    Mdouble gravityValue = -9.81; // g=9.81 m/s2 - 981 cm/s2
    //
    const char *simName = "FunnelExperiment";//"SpreadingProcess";
    //
    unsigned int saveCount = 10000; // 5000
    //
    //-------------------------------------------------------------------------------particles / PSD-------------------------------------------------------------------
    Mdouble meanRadius1 = 19.0e-6;//1.8577e-05;//0.0019; // in cm
    Mdouble MinRadius = 12.0e-6/2.0;//(13.183e-6/2.0); //in cm
    Mdouble MaxRadius = 79.0e-6/2.0;//(79.433e-6/2.0); //in cm
    // --------------------------------log-normal distribution mean and std, calculated from data points via lognfit matlab function-----------------------------------
    Mdouble LogmeanRadius = -10.8936;//-6.2884; //log(meanRadius1);
    Mdouble LogstdRadius = 0.3384; //log(stdRadius1);
    //
    //------------------------------------------------------------------------------ Domain ------------------------------------------------------------------------------
    Mdouble XMinDomain = -3e-3;
    Mdouble XMaxDomain = 3e-3;
    Mdouble YMinDomain = -3e-3;
    Mdouble YMaxDomain = 3e-3;
    Mdouble ZMinDomain = -100e-3;
    Mdouble ZMaxDomain = 5e-3;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //
    Mdouble FunnelMinRadius = 0.2e-3;
    Mdouble FunnelMaxRadius = 2e-3;//FunnelMinRadius*10;
    Mdouble FunnelHeight = FunnelMaxRadius/mathsFunc::tan(60. * constants::pi / 180.);//(FunnelMaxRadius - FunnelMinRadius) / mathsFunc::sin(60. * constants::pi / 180.); // /funnel angle
    Vec3D FunnelPointOnAxis = Vec3D(0.5 * (getXMax() + getXMin()), 0.5 * (getYMax() + getYMin()), 0);
    //
    Mdouble baseSideLength = 2.5e-3/4;
    Mdouble baseBottomePoint = -2e-3;
    Mdouble baseTopPoint = -0.1e-3;
    Mdouble lideSideLength = 2.75e-3/4;
    Mdouble lidBottomePoint = baseTopPoint;//-0.05e-3;
    //BottomWAll:
    Vec3D minPointBWall1 = Vec3D(-baseSideLength,-baseSideLength,baseBottomePoint);//Vec3D(-2.5e-3/2,-2.5e-3/2,-2e-3);
    Vec3D maxPointBWall1 = Vec3D(baseSideLength,baseSideLength,baseTopPoint);//-0.1e-3);//Vec3D(2.5e-3/2,2.5e-3/2,-0.1e-3);
    //
    // lids:
    Vec3D minPointLidLeft = Vec3D(-lideSideLength,-lideSideLength,lidBottomePoint);
    Vec3D maxPointLidLeft = Vec3D(-baseSideLength,baseSideLength,0);//-lideSideLength,lideSideLength,0);
    //
    Vec3D minPointLidFront = Vec3D(-baseSideLength,-lideSideLength,lidBottomePoint);
    Vec3D maxPointLidFront = Vec3D(baseSideLength,-baseSideLength,0);
    //
    Vec3D minPointLidRight = Vec3D(baseSideLength,-baseSideLength,lidBottomePoint);
    Vec3D maxPointLidRight = Vec3D(lideSideLength,lideSideLength,0);
    //
    Vec3D minPointLidBack = Vec3D(-baseSideLength,baseSideLength,lidBottomePoint);
    Vec3D maxPointLidBack = Vec3D(lideSideLength,lideSideLength,0);
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Mdouble minXparticles = -1e-3;
    Mdouble maxXparticles = 1e-3;
    Mdouble minYparticles = -1e-3;
    Mdouble maxYparticles = 1e-3;
    //---------------------------------------------------------------- the volume to be added: ------------------------------------------------------------------------
    Mdouble addVolume0 = 3e-10;//3*3e-10;
    //----------------------------------------------------------------------------------------------------
    void setupInitialConditions() override {
        //---------------------- Species ----------------------------------------
        //------------------- walls ----------------------------------
        //
        //
        IntersectionOfWalls base;
        base.setSpecies(speciesHandler.getObject(0));
        base.addObject(Vec3D(1.0,0.0,0.0),minPointBWall1);
        base.addObject(Vec3D(0.0,1.0,0.0),minPointBWall1);
        base.addObject(Vec3D(0.0,0.0,1.0),minPointBWall1);
        base.addObject(Vec3D(-1.0,0.0,0.0),maxPointBWall1);
        base.addObject(Vec3D(0.0,-1.0,0.0),maxPointBWall1);
        base.addObject(Vec3D(0.0,0.0,-1.0),maxPointBWall1);
        wallHandler.copyAndAddObject(base);
        //
        //for a prism wall, define points ABC (or more points)
        //as the corners of the prism base (ordered in clockwise
        //direction) and D as the 'periodic direction' of the prism
        //(upwards: Vec3D::Cross(A-B,D) has to point into the wall)
        Vec3D A(FunnelMinRadius, 0, 0);//2e-3 - FunnelHeight);
        Vec3D B(FunnelMaxRadius, 0, FunnelHeight);//2e-3);
        Vec3D C(FunnelMaxRadius, 0, 0);//2e-3 - FunnelHeight);
        Vec3D D(0, 1, 0); //Periodic direction of the prism
        //set walls
        AxisymmetricIntersectionOfWalls funnel;
        funnel.setSpecies(speciesHandler.getObject(1));
        funnel.setPosition(FunnelPointOnAxis);
        funnel.setAxis(Vec3D(0, 0, 1));
        funnel.addObject(Vec3D::cross(A - B, D), A);
        funnel.addObject(Vec3D::cross(B - C, D), B);
        funnel.addObject(Vec3D::cross(C - A, D), C);
        wallHandler.copyAndAddObject(funnel);
        //
        // LID
        IntersectionOfWalls lidL;
        lidL.setSpecies(speciesHandler.getObject(0));
        lidL.addObject(Vec3D(1.0,0.0,0.0),minPointLidLeft);
        lidL.addObject(Vec3D(0.0,1.0,0.0),minPointLidLeft);
        lidL.addObject(Vec3D(0.0,0.0,1.0),minPointLidLeft);
        lidL.addObject(Vec3D(-1.0,0.0,0.0),maxPointLidLeft);
        lidL.addObject(Vec3D(0.0,-1.0,0.0),maxPointLidLeft);
        lidL.addObject(Vec3D(0.0,0.0,-1.0),maxPointLidLeft);
        wallHandler.copyAndAddObject(lidL);
        //
        IntersectionOfWalls lidF;
        lidF.setSpecies(speciesHandler.getObject(0));
        lidF.addObject(Vec3D(1.0,0.0,0.0),minPointLidFront);
        lidF.addObject(Vec3D(0.0,1.0,0.0),minPointLidFront);
        lidF.addObject(Vec3D(0.0,0.0,1.0),minPointLidFront);
        lidF.addObject(Vec3D(-1.0,0.0,0.0),maxPointLidFront);
        lidF.addObject(Vec3D(0.0,-1.0,0.0),maxPointLidFront);
        lidF.addObject(Vec3D(0.0,0.0,-1.0),maxPointLidFront);
        wallHandler.copyAndAddObject(lidF);
        //
        IntersectionOfWalls lidR;
        lidR.setSpecies(speciesHandler.getObject(0));
        lidR.addObject(Vec3D(1.0,0.0,0.0),minPointLidRight);
        lidR.addObject(Vec3D(0.0,1.0,0.0),minPointLidRight);
        lidR.addObject(Vec3D(0.0,0.0,1.0),minPointLidRight);
        lidR.addObject(Vec3D(-1.0,0.0,0.0),maxPointLidRight);
        lidR.addObject(Vec3D(0.0,-1.0,0.0),maxPointLidRight);
        lidR.addObject(Vec3D(0.0,0.0,-1.0),maxPointLidRight);
        wallHandler.copyAndAddObject(lidR);
        //
        IntersectionOfWalls lidB;
        lidB.setSpecies(speciesHandler.getObject(0));
        lidB.addObject(Vec3D(1.0,0.0,0.0),minPointLidBack);
        lidB.addObject(Vec3D(0.0,1.0,0.0),minPointLidBack);
        lidB.addObject(Vec3D(0.0,0.0,1.0),minPointLidBack);
        lidB.addObject(Vec3D(-1.0,0.0,0.0),maxPointLidBack);
        lidB.addObject(Vec3D(0.0,-1.0,0.0),maxPointLidBack);
        lidB.addObject(Vec3D(0.0,0.0,-1.0),maxPointLidBack);
        wallHandler.copyAndAddObject(lidB);
        //
        //
// Container:
        /*Vec3D minPointLeftSide = Vec3D(-2e-3,-2e-3,2e-3);
        Vec3D maxPointLeftSide = Vec3D(-2e-3+1e-5,2e-3,4e-3);
        //
        Vec3D minPointFront = Vec3D(-2e-3 + 1e-5,-2e-3,2e-3);
        Vec3D maxPointFront = Vec3D(2e-3,-2e-3 + 1e-5,4e-3);
        AxisymmetricIntersectionOfWalls w;
        w.setSpecies(speciesHandler.getObject(0));
        w.setPosition(Vec3D(0,0,4e-3));
        w.setAxis(Vec3D(0,0,1));
        w.addObject(Vec3D(1,0,0), Vec3D(2e-3,0,0));  //Cylindric wall
        wallHandler.copyAndAddObject(w);
        //
        IntersectionOfWalls sideLeft ;//w2;
        sideLeft.setSpecies(speciesHandler.getObject(0));
        sideLeft.addObject(Vec3D(1.0,0.0,0.0),minPointLeftSide);
        sideLeft.addObject(Vec3D(0.0,1.0,0.0),minPointLeftSide);
        sideLeft.addObject(Vec3D(0.0,0.0,1.0),minPointLeftSide);
        sideLeft.addObject(Vec3D(-1.0,0.0,0.0),maxPointLeftSide);
        sideLeft.addObject(Vec3D(0.0,-1.0,0.0),maxPointLeftSide);
        sideLeft.addObject(Vec3D(0.0,0.0,-1.0),maxPointLeftSide);
        wallHandler.copyAndAddObject(sideLeft);
        //
        IntersectionOfWalls front ;//w3;
        front.setSpecies(speciesHandler.getObject(0));
        front.addObject(Vec3D(1.0,0.0,0.0),minPointFront);
        front.addObject(Vec3D(0.0,1.0,0.0),minPointFront);
        front.addObject(Vec3D(0.0,0.0,1.0),minPointFront);
        front.addObject(Vec3D(-1.0,0.0,0.0),maxPointFront);
        front.addObject(Vec3D(0.0,-1.0,0.0),maxPointFront);
        front.addObject(Vec3D(0.0,0.0,-1.0),maxPointFront);
        wallHandler.copyAndAddObject(front);*/
        //----------------------------------------------------------------------- Particles insertion ------------------------------------------------------------------
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
        Mdouble fillHeight0 = 0.0;
        while (addVolume0>0) {
            Mdouble x = random.getRandomNumber(minXparticles, maxXparticles);
            Mdouble y = random.getRandomNumber(minYparticles, maxYparticles); //(ToolThickness, (ToolWidth-ToolThickness));
            Mdouble z = random.getRandomNumber(1e-4, 1e-4+fillHeight0);
            p0.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p0)) {
                particleHandler.copyAndAddObject(p0);
                addVolume0 -= p0.getVolume();
                do {
                    p0.setRadius(d(gen));
                } while (p0.getRadius()< MinRadius || p0.getRadius() > MaxRadius);
                if (particleHandler.getNumberOfObjects() % 100 == 0) std::cout << '.' << std::flush;
                if (particleHandler.getNumberOfObjects() % 1000 == 0) std::cout << ' ';
                if (particleHandler.getNumberOfObjects() % 10000 == 0) std::cout << addVolume0 << '\n';
            } else {
                fillHeight0 += 1e-3*meanRadius1; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        //
        logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());
        //
        DeletionBoundary deletionBoundaryZ;
        deletionBoundaryZ.set(Vec3D(0, 0, -1), 2e-3);
        boundaryHandler.copyAndAddObject(deletionBoundaryZ);
        //
    }
    //
    //----------------------------------------------------------------------------------------------------------------------------------------------
    void actionsAfterTimeStep() override {
        Mdouble baseSpeed = 2e-3;//10e-3;
        // Start Exp
            static bool moveTool = false;
            //if (!moveTool && getTime() >= TimeMoveTool)
            if(!moveTool && getTime() >= 0.1) //getKineticEnergy()<1e-5*getElasticEnergy())
            {
                logger(INFO,"Base is Moving");
                wallHandler.getObject(0)->setVelocity(Vec3D(0, 0, -baseSpeed));
                //
                wallHandler.getObject(2)->setVelocity(Vec3D(0, 0, -baseSpeed));
                wallHandler.getObject(3)->setVelocity(Vec3D(0, 0, -baseSpeed));
                wallHandler.getObject(4)->setVelocity(Vec3D(0, 0, -baseSpeed));
                wallHandler.getObject(5)->setVelocity(Vec3D(0, 0, -baseSpeed));
                moveTool = true;
            }
            //
            if(getTime()>=0.6)
            {
                wallHandler.getObject(0)->setVelocity(Vec3D(0, 0, 0));
                //
                // bug?? if I stop the walls then try to remove them -> Error, why?
/*                wallHandler.getObject(2)->setVelocity(Vec3D(0, 0, 0));
                wallHandler.getObject(3)->setVelocity(Vec3D(0, 0, 0));
                wallHandler.getObject(4)->setVelocity(Vec3D(0, 0, 0));
                wallHandler.getObject(5)->setVelocity(Vec3D(0, 0, 0));*/
                //
            }
            //
        static bool lidRemoved = false;
            if(!lidRemoved && getTime()>0.65)
            {
                wallHandler.removeObject(5);
                wallHandler.removeObject(4);
                wallHandler.removeObject(3);
                wallHandler.removeObject(2);
                lidRemoved = true;
            }
/*            Mdouble timeStarted = getTime();
            //
            if(moveTool && getTime() > timeStarted && getKineticEnergy()<1e-5*getElasticEnergy())
                logger(INFO,"time %, stop here", getTime());*/
            }
    //
};
// --------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc UNUSED, char *argv[] UNUSED)
{
    //logger(INFO,"Simple box for creating particles");
    FunnelExperiment FunnelExp;
    //
    //------------------------------------------------------------------------------Parametric friction:------------------------------------------------------------------------------
    FunnelExp.autoNumber();
    //
    std::vector<int> studyNum=FunnelExp.get2DParametersFromRunNumber(FunnelExp.numOfStudies1,FunnelExp.numOfStudies2);
    //
    Mdouble SlidingFrictionCoeff = FunnelExp.SlidingFrictionCoeffVector[studyNum[1]-1];
    Mdouble RollingFrictionCoeff = FunnelExp.RollingFrictionCoeffVector[studyNum[2]-1];
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    //
    FunnelExp.setName(FunnelExp.simName);
    FunnelExp.setSystemDimensions(3);
    FunnelExp.setGravity(Vec3D(0.0,0.0,FunnelExp.gravityValue));
    FunnelExp.setXMin(FunnelExp.XMinDomain);
    FunnelExp.setYMin(FunnelExp.YMinDomain);
    FunnelExp.setZMin(FunnelExp.ZMinDomain);
    FunnelExp.setXMax(FunnelExp.XMaxDomain);
    FunnelExp.setYMax(FunnelExp.YMaxDomain);
    FunnelExp.setZMax(FunnelExp.ZMaxDomain);
    FunnelExp.setTimeMax(FunnelExp.MaxSimTime);
    //FunnelExp.setHGridMaxLevels(2);
    //
    //auto species = FunnelExp.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    auto species = FunnelExp.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionJKRAdhesiveSpecies());
    //
    species->setDensity(FunnelExp.MatDensity);
    //
    // Here we assign the mass of the smallest particle to calculate the stiffiness and disppation as we want to get ------
    Mdouble mass = species->getMassFromRadius(FunnelExp.minRaduis);
    // ----- set stiffness and dissipation from tc,COR,mass
    species->setCollisionTimeAndRestitutionCoefficient(FunnelExp.tc,FunnelExp.restitutionCoeff,mass);
    //------------------------------------------------------------------------------ Cohesion Parameters ------------------------------------
    if (FunnelExp.cohesiveON)
    {
        Mdouble adhesionStiffness = FunnelExp.AdhStiffnessFactor * species->getStiffness();
        species->setAdhesionStiffness(adhesionStiffness);
        species->setSurfaceEnergy(FunnelExp.surfaceEnergy);
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
    //------------------------------------------------------------------------------ P-Funnel species -------------------------------------------------------------
    auto FunnelSpecies = FunnelExp.speciesHandler.copyAndAddObject(species);
    FunnelSpecies->setSurfaceEnergy(0.0);
    // --------- added later - 27 Sept 2020 --------
    //FunnelSpecies->setSlidingFrictionCoefficient(0.0);
    //FunnelSpecies->setRollingFrictionCoefficient(0.0);
    //----------------------------------------------
    auto FunnelParticleSpecies = FunnelExp.speciesHandler.getMixedObject(species,FunnelSpecies);
    //
    //FunnelParticleSpecies->setSlidingStiffness(2./7.*species->getStiffness());
    //FunnelParticleSpecies->setSlidingDissipation(2./7.*species->getDissipation());
    FunnelParticleSpecies->setSlidingFrictionCoefficient(SlidingFrictionCoeff);
    //
    //FunnelParticleSpecies->setRollingStiffness(2./5.*species->getStiffness());
    //FunnelParticleSpecies->setRollingDissipation(2./5.*species->getDissipation());
    FunnelParticleSpecies->setRollingFrictionCoefficient(RollingFrictionCoeff);
    //
    FunnelParticleSpecies->setSurfaceEnergy(0.0);
    //------------------------------------------------------------------------------Data Output--------------------------------------------------------------------------------------------
    FunnelExp.setSaveCount(FunnelExp.saveCount);//50
    FunnelExp.dataFile.setFileType(FileType::ONE_FILE);
    //FunnelExp.fStatFile.setFileType(FileType::NO_FILE);
    FunnelExp.setParticlesWriteVTK(true);
    FunnelExp.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    //
    FunnelExp.setTimeStep(FunnelExp.timeStep);
    //
    FunnelExp.solve();//(argc, argv);
    //
    return 0;
}