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
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include "Walls/AxisymmetricIntersectionOfWalls.h"
#include <random>
#include "MercuryTime.h"
#include "Boundaries/DeletionBoundary.h"

#include "Species/LinearViscoelasticFrictionSpecies.h"
//#include "Species/LinearViscoelasticFrictionJKRAdhesiveSpecies.h"


class SpreadingLayerInsertAndSpread: public Mercury3D
{
private:

    // simulation files name
    const char *simName = "SpreadingProcessInsertion";

    // save data after each number of time steps:
    unsigned int saveCount = 10000;

    // to delete particles out of the domain along the spreading direction
    bool deleteBoundary = false;
    // only insert particles ot insert and spread (false)
    bool insertParticlesOnly = false;

    //  PROCESS PARAMETERS
    bool roller = false; // true -> roller , false -> blade
    Mdouble toolSpeed = 0.01; // m/s
    // Gap between tool and bottom wall
    Mdouble Gap = 100e-6; //100 microns -> 100e-6 m ;

    //  Material PARAMETERS
    // set time taking into account the time needed for the system to relax
    Mdouble maxSimTime = 0.05;
    // time to remove frony wall after particles settle down
    Mdouble TimeRemoveWall = 0.01;
    // time to move the spreading tool after the particles have settled down and the system is relaxed
    Mdouble TimeMoveTool = 0.05;

    Mdouble materialDensity = 4430;
    Mdouble tc = 90e-7; //calculated as tc = 0.005-0.01 * tg : tg=sqrt(d50/g) => 0.005*sqrt(35e-6/9.81)=
    Mdouble timeStep = 0.02*tc;
    Mdouble restitutionCoeff = 0.4;//0.1;
    Mdouble minRaduis = 12.0e-6;

    //Material PARAMETERS - cohesion
    //
    //bool cohesiveON = false;
    //Mdouble AdhStiffnessFactor = 0.5;
    //Mdouble surfaceEnergy = 0.1e-3; //J/m2

    // Parametric study parameters
    // seed number for particle generation
    int particleGenerationSEED = 1;

    // number of studies for each parameter
    int numOfStudies1 = 7;
    int numOfStudies2 = 6;
    // parametric study of friction coeffecients
    std::vector<double> SlidingFrictionCoeffVector = {0.5,0.4,0.3,0.25,0.2,0.1,0.05};
    std::vector<double> RollingFrictionCoeffVector = {0.4,0.3,0.2,0.1,0.05,0.005};


    Mdouble gravityValue = -9.81;


    // particles sizes / PSD
    Mdouble meanRadius1 = 19.0e-6;
    Mdouble MinRadius = 12.0e-6/2.0;
    Mdouble MaxRadius = 79.0e-6/2.0;
    // log-normal distribution mean and std, calculated from data points via lognfit matlab function
    Mdouble LogmeanRadius = -10.8936;
    Mdouble LogstdRadius = 0.3384;


    //
    //  spatial and domain  parameters
    Mdouble scale = 1;
    //TOOL - initial configuration setup - unit: m:
    Mdouble SpreadingLength = 0.01*scale;
    Mdouble SpreadingWidth  = 0.001*scale;
    Mdouble ToolHight  = 0.003*scale;
    Mdouble WallHight = ToolHight;
    Mdouble ToolThickness = 500e-6*scale;
    Mdouble rollerRadius = ToolThickness*2.0;

    // Domain
    Mdouble XMaxDomain = SpreadingLength+ToolThickness;
    Mdouble YMaxDomain = SpreadingWidth;
    Mdouble ZMinDomain = -ToolThickness;
    Mdouble ZMaxDomain = ToolHight;


    //TOOL
    // ---- tool Balde:
    Vec3D minPoint_tool = Vec3D(0.0,0.0,Gap);
    Vec3D maxPoint_tool = Vec3D(ToolThickness,SpreadingWidth,ToolHight);
    // ---- tool Roller:
    Vec3D rollerPosition = Vec3D(-ToolThickness,0.0,rollerRadius+Gap);

    //Back Dam, small wall under tool
    Vec3D minPoint_back_dam = Vec3D(-ToolThickness,0.0,0.0);
    Vec3D maxPoint_back_dam = Vec3D(ToolThickness,SpreadingWidth,Gap);
    //front dam, present while inserting particles -> removed during spreading
    Mdouble damThickness = SpreadingLength/5.0+ToolThickness;// -> 2mm insertion domain
    Vec3D minPoint_dam = Vec3D(damThickness,0.0,0.0);
    Vec3D maxPoint_dam = Vec3D(damThickness+ToolThickness,SpreadingWidth,WallHight);

    //
    Mdouble minXparticles = ToolThickness;
    Mdouble maxXparticles = SpreadingLength/5.0+ToolThickness;
    Mdouble maxYparticles = SpreadingWidth;

    // the volume of particles to be added:
    // correspond to spreading volume of desired layer 5 x 1 x 0.1 mm + initial insertion volume = 7 x 1 x 0.1 mm
    Mdouble addVolume0 = (SpreadingLength-(SpreadingLength/2.0))*SpreadingWidth*Gap + (SpreadingLength/5.0)*SpreadingWidth*Gap;

public:
    void setupInitialConditions() override
    {
        // Parametric friction:
        autoNumber();
        std::vector<int> studyNum=get2DParametersFromRunNumber(numOfStudies1,numOfStudies2);
        Mdouble SlidingFrictionCoeff = SlidingFrictionCoeffVector[studyNum[1]-1];
        Mdouble RollingFrictionCoeff = RollingFrictionCoeffVector[studyNum[2]-1];

        //
        setName(simName);
        setSystemDimensions(3);
        setGravity(Vec3D(0.0,0.0,gravityValue));
        setZMin(ZMinDomain);
        setXMax(XMaxDomain);
        setYMax(YMaxDomain);
        setZMax(ZMaxDomain);
        setTimeMax(maxSimTime);

        //
        auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
        // this species uses adhesion in terms of JKR pull-off force
        //auto species = speciesHandler.copyAndAddObject(LinearViscoelasticFrictionJKRAdhesiveSpecies());


        species->setDensity(materialDensity);
        // Here we assign the mass of the smallest particle to calculate the stiffiness and disppation
        Mdouble mass = species->getMassFromRadius(minRaduis);
        // ----- set stiffness and dissipation from tc,COR,mass
        species->setCollisionTimeAndRestitutionCoefficient(tc,restitutionCoeff,mass);


        // Cohesion Parameters
        /*if (cohesiveON)
        {
            Mdouble adhesionStiffness = AdhStiffnessFactor * species->getStiffness();
            species->setAdhesionStiffness(adhesionStiffness);
            species->setSurfaceEnergy(surfaceEnergy);
        }*/
        //

        // Sliding and rolling stiffness and dissipation
        species->setSlidingStiffness(2./7.*species->getStiffness());
        species->setSlidingDissipation(2./7.*species->getDissipation());
        species->setSlidingFrictionCoefficient(SlidingFrictionCoeff);
        //
        species->setRollingStiffness(2./5.*species->getStiffness());
        species->setRollingDissipation(2./5.*species->getDissipation());
        species->setRollingFrictionCoefficient(RollingFrictionCoeff);

        // Data Output
        setSaveCount(saveCount);//50
        dataFile.setFileType(FileType::ONE_FILE);
        setParticlesWriteVTK(true);
        setWallsWriteVTK(FileType::MULTIPLE_FILES);
        setTimeStep(timeStep);


        //Bottom Wall
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0,0.0,0.0));
        wallHandler.copyAndAddObject(w0);

        //TOOLs:
        // Blade:
        IntersectionOfWalls w1;
        w1.setSpecies(speciesHandler.getObject(0)); //using different species for the Tool: w0.setSpecies(speciesHandler.getObject(1));
        w1.addObject(Vec3D(1.0,0.0,0.0),minPoint_tool);
        w1.addObject(Vec3D(0.0,1.0,0.0),minPoint_tool);
        w1.addObject(Vec3D(0.0,0.0,1.0),minPoint_tool);
        w1.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_tool);
        w1.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_tool);
        w1.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_tool); //w1.addObject(Vec3D(cos(Angel),0.0,sin(Angel)),Vec3D(AngelDim,0.0,Gap)); //in next line
        wallHandler.copyAndAddObject(w1);
        if(roller)
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
        //mid Dam 1
        IntersectionOfWalls w3;
        w3.setSpecies(speciesHandler.getObject(0));
        w3.addObject(Vec3D(1.0,0.0,0.0),minPoint_dam);
        w3.addObject(Vec3D(0.0,1.0,0.0),minPoint_dam);
        w3.addObject(Vec3D(0.0,0.0,1.0),minPoint_dam);
        w3.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_dam);
        w3.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_dam);
        w3.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_dam);
        wallHandler.copyAndAddObject(w3);

        // PBC in Y-dir
        PeriodicBoundary b0;
        b0.set(Vec3D(0, 1, 0), 0.0, SpreadingWidth);
        boundaryHandler.copyAndAddObject(b0);

        // Particles insertion
        std::mt19937 gen;
        // set seed to avoid randon packing for paramteric study i.e. every sim will have same initial packing
        gen.seed(particleGenerationSEED);
        std::lognormal_distribution<> d(LogmeanRadius, LogstdRadius);
        //add particles until the volume to be added is zero
        SphericalParticle p0;
        p0.setSpecies(speciesHandler.getObject(0));
        p0.setRadius(meanRadius1);
        Mdouble fillHeight0 = 0.0;
        while (addVolume0>0) {
            Mdouble x = random.getRandomNumber(minXparticles, maxXparticles);
            Mdouble y = random.getRandomNumber(0.0, maxYparticles);
            Mdouble z = random.getRandomNumber(0.0, fillHeight0);
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
        logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());


        // delete particles outside domain of interest:
        if(deleteBoundary)
        {
            DeletionBoundary deletionBoundary1;
            deletionBoundary1.set(Vec3D(1, 0, 0), (SpreadingLength - (SpreadingLength / 5.0)) + 4.0*ToolThickness);//(SpreadingLength-0.005)+ToolThickness);
            boundaryHandler.copyAndAddObject(deletionBoundary1);
            DeletionBoundary deletionBoundary2;
            deletionBoundary2.set(Vec3D(-1, 0, 0), ToolThickness);//(SpreadingLength-0.005)+ToolThickness);
            boundaryHandler.copyAndAddObject(deletionBoundary2);
        }
    }


    void actionsAfterTimeStep() override
    {

        // Remove dams
        static bool wallRemoved = false;
        if (wallRemoved == false && getTime() >= TimeRemoveWall) {
            logger(INFO,"walls removed");
            if (!roller) {
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
                logger(INFO,"Tool is Moving");
                if (!roller)
                {
                    //blade:
                    wallHandler.getObject(1)->setVelocity(Vec3D(toolSpeed, 0, 0));
                }
                else
                {
                    //roller:
                    wallHandler.getLastObject()->setVelocity(Vec3D(toolSpeed, 0, 0));
                    wallHandler.getLastObject()->setAngularVelocity({0, -toolSpeed / rollerRadius,0});//getObject(1)->setAngularVelocity({0,-tool_Speed/rollerRadius,0});
                }
                moveTool = true;
            }
        }
    }

    // to print data on terminal
    /*void printTime() const override
    {
        logger(INFO,"t=%, tMax=%, N=%", getTime(),getTimeMax(), particleHandler.getSize());
    }*/
};


int main(int argc UNUSED, char *argv[] UNUSED)
{
    Time time;
    time.tic();

    SpreadingLayerInsertAndSpread powderSpreading;
    powderSpreading.solve();

    // record cpu time:
    std::ofstream outTime ("simTime"+powderSpreading.getName(), std::ios::out);
    outTime << time.toc() << std::endl;
    outTime.close();

    return 0;
}