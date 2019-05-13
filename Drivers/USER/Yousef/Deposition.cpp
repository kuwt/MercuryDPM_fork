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

//
// Created by yousef on 8/24/17.
//
#include "Mercury3D.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
////#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"



class Deposition: public Mercury3D
{
public:

    void setupInitialConditions() {

        //
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0,0.0,0.0));
        wallHandler.copyAndAddObject(w0);
        //TOOL:
        Mdouble scale = 1e-2; //1e-1; //copy to speed and simulation boundaries
        //
        Mdouble ToolLength = 45e-1; // unit cm
        Mdouble ToolWidth = 30e-1;
        Mdouble ToolHight = 20e-1;
        Mdouble ToolThickness = 5e-1;
        Mdouble Gap = 0.1e-1;
        Mdouble Angel = constants::pi/180.0*30;
        Mdouble AngelDim = scale*4e-1;
        //
        //side wall
        Vec3D minPoint1_tool = Vec3D(0.0,0.0,0.0);
        Vec3D maxPoint1_tool = Vec3D(scale*ToolLength,scale*ToolThickness,scale*ToolHight); //Vec3D(getXMax(),0.05*getYMax(),getZMax())
        //back wall
        Vec3D minPoint2_tool = Vec3D(0.0,0.0,Gap);
        Vec3D maxPoint2_tool = Vec3D(scale*ToolThickness,scale*(ToolWidth-ToolThickness),scale*ToolHight); //Vec3D(0.2*getXMax(),getYMax(),getZMax())
        //side wall
        Vec3D minPoint3_tool = Vec3D(0.0,scale*(ToolWidth-ToolThickness),0.0);
        Vec3D maxPoint3_tool = Vec3D(scale*ToolLength,scale*ToolWidth,scale*ToolHight); // Vec3D(getXMax(),getYMax(),getZMax())
        //
        //in z-x plane - tool
        IntersectionOfWalls w1;
        w1.setSpecies(speciesHandler.getObject(0));
        w1.addObject(Vec3D(1.0,0.0,0.0),minPoint1_tool);
        w1.addObject(Vec3D(0.0,1.0,0.0),minPoint1_tool);
        w1.addObject(Vec3D(0.0,0.0,1.0),minPoint1_tool);
        w1.addObject(Vec3D(-1.0,0.0,0.0),maxPoint1_tool);
        w1.addObject(Vec3D(0.0,-1.0,0.0),maxPoint1_tool);
        w1.addObject(Vec3D(0.0,0.0,-1.0),maxPoint1_tool);
        wallHandler.copyAndAddObject(w1);
        //
        IntersectionOfWalls w2;
        w2.setSpecies(speciesHandler.getObject(0));
        w2.addObject(Vec3D(1.0,0.0,0.0),minPoint2_tool);
        w2.addObject(Vec3D(0.0,1.0,0.0),minPoint2_tool);
        w2.addObject(Vec3D(0.0,0.0,1.0),minPoint2_tool);
        w2.addObject(Vec3D(-1.0,0.0,0.0),maxPoint2_tool);
        w2.addObject(Vec3D(0.0,-1.0,0.0),maxPoint2_tool);
        w2.addObject(Vec3D(0.0,0.0,-1.0),maxPoint2_tool);
        w2.addObject(Vec3D(cos(Angel),0.0,sin(Angel)),Vec3D(AngelDim,0.0,Gap));
        wallHandler.copyAndAddObject(w2);
        //
        IntersectionOfWalls w3;
        w3.setSpecies(speciesHandler.getObject(0));
        w3.addObject(Vec3D(1.0,0.0,0.0),minPoint3_tool);
        w3.addObject(Vec3D(0.0,1.0,0.0),minPoint3_tool);
        w3.addObject(Vec3D(0.0,0.0,1.0),minPoint3_tool);
        w3.addObject(Vec3D(-1.0,0.0,0.0),maxPoint3_tool);
        w3.addObject(Vec3D(0.0,-1.0,0.0),maxPoint3_tool);
        w3.addObject(Vec3D(0.0,0.0,-1.0),maxPoint3_tool);
        wallHandler.copyAndAddObject(w3);
        //
        //particles container - deleted
        //
        //Front Dam
        Mdouble damThickness = ToolLength+5.0e-1;
        //
        Vec3D minPoint_dam = Vec3D(scale*ToolLength,0.0,0.0);
        Vec3D maxPoint_dam = Vec3D(scale*damThickness,scale*ToolWidth,scale*ToolHight);
        //
        //
        IntersectionOfWalls w4;
        w4.setSpecies(speciesHandler.getObject(0));
        w4.addObject(Vec3D(1.0,0.0,0.0),minPoint_dam);
        w4.addObject(Vec3D(0.0,1.0,0.0),minPoint_dam);
        w4.addObject(Vec3D(0.0,0.0,1.0),minPoint_dam);
        w4.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_dam);
        w4.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_dam);
        w4.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_dam);
        wallHandler.copyAndAddObject(w4);
        //
        //particles:
        //
        Mdouble MinRadius = 10e-4; //diameter = 20//10microM = 10e-3 mm = 10e-4 cm
        Mdouble MaxRadius = 27.5e-4; // diameter = 55 // 27.5microM = 27.5e-3 mm = 27.5e-4cm
        Vec3D minPoint_particles = Vec3D(scale*ToolThickness,scale*ToolThickness,0.0);
        Vec3D maxPoint_particles = Vec3D(scale*ToolLength,scale*(ToolWidth-ToolThickness),scale*ToolHight);
        //
        SphericalParticle insertionBoundaryParticle;
        insertionBoundaryParticle.setSpecies(speciesHandler.getObject(0));
        //CubeInsertionBoundary::set(BaseParticle* particleToCopy, int maxFailed, Vec3D posMin, Vec3D posMax, Vec3D velMin, Vec3D velMax, double radMin, double radMax)
        CubeInsertionBoundary insertionBoundary;
        insertionBoundary.set(&insertionBoundaryParticle,0, minPoint_particles, maxPoint_particles, Vec3D(0, 0, 0), Vec3D(0, 0, 0), MinRadius, MaxRadius);
        int count = 0;
        while (particleHandler.getNumberOfObjects()<700)
        {
            if (++count%10==0) std::cout << particleHandler.getNumberOfObjects() << std::endl;
            insertionBoundary.checkBoundaryBeforeTimeStep(this);
        }
        //boundaryHandler.copyAndAddObject(insertionBoundary);
        //particleHandler.getVolume();
        //logger(INFO,"particles number= ",insertionBoundary.getNumberOfParticlesInserted());
    }
    //
    //

    void actionsAfterTimeStep() override {
        //
        static bool wallRemoved = false;
        if (wallRemoved==false && getTime()>=0.18) {
            //wallHandler.getObject(7)->setPosition(Vec3D(0,0,4));
            logger(INFO,"wall removed");
            wallHandler.removeObject(4);
            wallRemoved = true;
        }
        //
        Mdouble tool_Speed = 1e-2*200e-1; //cm/s
        static bool moveTool = false;
        if (!moveTool && getTime()>= 0.185) {
            logger(INFO,"Tool is Moving");
            wallHandler.getObject(1)->setVelocity(Vec3D(tool_Speed,0,0));
            wallHandler.getObject(2)->setVelocity(Vec3D(tool_Speed,0,0));
            wallHandler.getObject(3)->setVelocity(Vec3D(tool_Speed,0,0));
            moveTool = true;
        }
    }
    //
    void printTime() const override
    {
        logger(INFO,"t=%, tMax=%, N=%", getTime(),getTimeMax(), particleHandler.getSize());
    }

   /* bool readNextArgument(int& i, int argc, char* argv[]) override
    {
        // The argument argv[i] identifies the label of the flag, and subsequent arguments (usually 1)
        // contain the content.
        //
        // For example...
        // Checks if the "-name" flag has been passed
        // The strcmp returns 0 if "argv[i]" is "-name" (i.e. !strcmp(argv[i], "-name") --> 1)
        // In this case, the setName function is run with the relevant input (i.e. the value which
        // immediately follows the "-name" flag
        if (!strcmp(argv[i], "-width"))
        {
            Mdouble width = atof(argv[i+1]);
            logger(INFO,"Setting the width to %",width);
            setYMax(width);
            setName(getName()+"W"+argv[i+1]);
            logger(INFO,"Name set to %",getName());
        }
        else
        {
            DPMBase::readNextArgument(i,argc,argv);
        }
    }*/
};

int main(int argc UNUSED, char *argv[] UNUSED)
{
    logger(INFO,"Simple box for creating particles");

    Deposition deposition_problem;

    Mdouble scaleMain = 1e-2;


    deposition_problem.setName("Deposition");
    deposition_problem.setSystemDimensions(3);
    deposition_problem.setGravity(Vec3D(0.0,0.0,-981));
    //deposition_problem.setXMin(0.0);
    //deposition_problem.setYMin(0.0);
    //deposition_problem.setZMin(0.0);
    deposition_problem.setXMax(scaleMain*(45.0e-1+5e-1)); // cm
    deposition_problem.setYMax(scaleMain*30.0e-1);
    deposition_problem.setZMax(scaleMain*100.0e-1);
    deposition_problem.setTimeMax(0.4); //0.4
    deposition_problem.setHGridMaxLevels(2);

    auto species = deposition_problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    species->setDensity(4.430); //g/cm3 - 4430kg/m3 //2000
    //species->setStiffness(10000); //Where the value came from:
    //species->setDissipation(100);
    //
    Mdouble mass = species->getMassFromRadius(10e-4);
    Mdouble tc = 2e-4; //0.01
    Mdouble restitutionCoeff = 0.1;
    species->setCollisionTimeAndRestitutionCoefficient(tc,restitutionCoeff,mass);
    //
    species->setSlidingStiffness(2./7.*species->getStiffness());
    species->setSlidingDissipation(2./7.*species->getDissipation());
    species->setSlidingFrictionCoefficient(0.5);
    //species->setSlidingFrictionCoefficientStatic(0.6);
    //
    species->setRollingStiffness(2./5.*species->getStiffness());
    species->setRollingDissipation(2./5.*species->getDissipation());
    species->setRollingFrictionCoefficient(0.1);
    //
    deposition_problem.setSaveCount(500);//50
    deposition_problem.setParticlesWriteVTK(true);
    deposition_problem.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    //
    deposition_problem.setTimeStep(0.02*tc); //0.02*tc
    //deposition_problem.setTimeStep(1e-4);

    //deposition_problem.autoNumber();

    deposition_problem.solve(argc, argv);
    return 0;
}
