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
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
//#include "Species/LinearViscoelasticSpecies.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
#include <random>

class Deposition: public Mercury3D
{
public:
    //
    Mdouble scale = 0.5e-1; //1e-1; //copy to speed and simulation boundaries
    //
    //
    //Mdouble TimeRemoveBackWall = 0.2;
    Mdouble TimeRemoveWall = 0.032;//0.2;
    Mdouble TimeMoveTool = 0.035;//0.25;
    Mdouble tool_Speed = 200e-1*scale;// = (1e-1)*200e-1; //cm/s
    //
    /*Mdouble random_number(Mdouble min, Mdouble max)
    {
        // use thread_local to make this function thread safe
        thread_local static std::mt19937 mt{std::random_device{}()};
        thread_local static std::lognormal_distribution<Mdouble> dist;
        using pick = std::lognormal_distribution<Mdouble>::param_type;

        return dist(mt, pick(min, max));
    }*/
    void setupInitialConditions() {

        //TOOL:
        Mdouble ToolLength = 42e-1*scale; // unit cm
        Mdouble ToolWidth = 30e-1*scale;
        Mdouble ToolHight = 25e-1*scale;
        Mdouble ToolThickness = 5e-1*scale;
        //For front wall:
        Mdouble FrontWallThickness = 18e-1*scale;
        Mdouble FrontWallDistance = ToolLength-FrontWallThickness;
        Mdouble Angel = constants::pi/180.0*30;
        Mdouble AngelDim = FrontWallDistance;//ToolLength-ToolThickness;//scale*5e-1;
        // For Scraper:
        Mdouble Gap = 0.1e-1; //0.1e-1; // keeping the gap as the actual dim requires scaling down the track volume by e-2 instead of (0.5e-1) which corresponds to 9000 particles
        //Mdouble TrackLength = 900e-1*(scale*scale);
        Mdouble ScraperDistance = 6e-1*scale;
        Mdouble ScraperWidth = 15e-1*scale;
        Mdouble ScraperActiveLength = 19.54e-1*scale; //(25-5.13030214426)*1e-1*scale;
        Mdouble ScraperAngle = constants::pi/180.0*70;
        //
        //
        //Front Dam
        Mdouble damThickness = ToolLength+(5.0e-1*scale);
        //
        //Perform parametric study:
        //
        /*if (DPMBase::getRunNumber() == 2)
         {
             ToolWidth = 40e-1*scale; // corresponds to 30e-1 for the cavity
         }
         //
         if(DPMBase::getRunNumber() == 3)
         {
             ToolWidth = 40e-1*scale; // corresponds to 30e-1 for the cavity
             tool_Speed = 200e-1;
         }
         //
         if(DPMBase::getRunNumber() == 4)
         {
             ToolWidth = 20e-1*scale; // corresponds to 10e-1 for the cavity
             tool_Speed = 200e-1*scale;
         }
         //
         if(DPMBase::getRunNumber() == 5)
         {
             ToolWidth = 20e-1*scale; // corresponds to 10e-1 for the cavity
             tool_Speed = 200e-1;
         }
         //
         if(DPMBase::getRunNumber() == 6)
         {
             scale = (0.5)*1e-1;
             tool_Speed = 200e-1*scale;
         }
         //
         if(DPMBase::getRunNumber() == 7)
         {
             scale = (0.5)*1e-1;
             tool_Speed = 200e-1;
         }*/
        //
        //TOOL
        //side wall
        Vec3D minPoint1_tool = Vec3D(0.0,0.0,0.0);
        Vec3D maxPoint1_tool = Vec3D(ToolLength,ToolThickness,ToolHight);
        //back wall
        Vec3D minPoint2_tool = Vec3D(FrontWallDistance,ToolThickness,0.0);//Vec3D(FrontWallDistance,0.0,0.0); //
        Vec3D maxPoint2_tool = Vec3D(ToolLength,(ToolWidth-ToolThickness),ToolHight);//Vec3D(ToolLength,(ToolWidth-(2*ToolThickness)),ToolHight); //
        //side wall
        Vec3D minPoint3_tool = Vec3D(0.0,(ToolWidth-ToolThickness),0.0);
        Vec3D maxPoint3_tool = Vec3D(ToolLength,ToolWidth,ToolHight);
        //
        //Scraper:
        Vec3D minPoint_scraper = Vec3D(ScraperDistance,ToolThickness,Gap);//Vec3D(ScraperDistance,0.0,Gap); //
        Vec3D maxPoint_scraper = Vec3D(ScraperDistance+ScraperWidth,ToolWidth-ToolThickness,ToolHight+Gap);//Vec3D(ScraperDistance+ScraperWidth,ToolWidth-(2*ToolThickness),ToolHight+Gap); //
        //Back Dam
        Vec3D minPoint_back_dam = Vec3D(0.0,ToolThickness,0.0);
        Vec3D maxPoint_back_dam = Vec3D(ToolThickness,(ToolWidth-ToolThickness),Gap);
        //
        //Bottom Wall
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(0)); //w0.setSpecies(speciesHandler.getObject(1));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0,0.0,0.0));
        wallHandler.copyAndAddObject(w0);
        //TOOL:
        //in z-x plane - tool / side wall

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
        /*
        InfiniteWall w1;
        w1.setSpecies(speciesHandler.getObject(1));
        w1.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0,0.0,0.0));
        wallHandler.copyAndAddObject(w1);
        */
        //Front wall:
        IntersectionOfWalls w2;
        w2.setSpecies(speciesHandler.getObject(0));
        w2.addObject(Vec3D(1.0,0.0,0.0),minPoint2_tool);
        w2.addObject(Vec3D(0.0,1.0,0.0),minPoint2_tool);
        w2.addObject(Vec3D(0.0,0.0,1.0),minPoint2_tool);
        w2.addObject(Vec3D(-1.0,0.0,0.0),maxPoint2_tool);
        w2.addObject(Vec3D(0.0,-1.0,0.0),maxPoint2_tool);
        w2.addObject(Vec3D(0.0,0.0,-1.0),maxPoint2_tool);
        w2.addObject(Vec3D(cos(Angel),0.0,-sin(Angel)),Vec3D(AngelDim,0.0,0.0));
        wallHandler.copyAndAddObject(w2);
        //side wall

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
        /*
        InfiniteWall w3;
        w3.setSpecies(speciesHandler.getObject(1));
        w3.set(Vec3D(0.0, -1.0, 0.0), Vec3D(0.0,ToolWidth-(2*ToolThickness),0.0));
        wallHandler.copyAndAddObject(w3);
        */
        //
        //scraper
        //
        IntersectionOfWalls w4;
        w4.setSpecies(speciesHandler.getObject(0));
        w4.addObject(Vec3D(1.0,0.0,0.0),minPoint_scraper);
        w4.addObject(Vec3D(0.0,1.0,0.0),minPoint_scraper);
        w4.addObject(Vec3D(0.0,0.0,1.0),minPoint_scraper);
        w4.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_scraper);
        w4.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_scraper);
        w4.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_scraper);
        w4.addObject(Vec3D(cos(ScraperAngle),0.0,sin(ScraperAngle)),Vec3D(ScraperDistance,0.0,(ToolHight+Gap)-ScraperActiveLength));
        wallHandler.copyAndAddObject(w4);
        /*
        //Back Dam
        IntersectionOfWalls w5;
        w5.setSpecies(speciesHandler.getObject(0));
        w5.addObject(Vec3D(1.0,0.0,0.0),minPoint_back_dam);
        w5.addObject(Vec3D(0.0,1.0,0.0),minPoint_back_dam);
        w5.addObject(Vec3D(0.0,0.0,1.0),minPoint_back_dam);
        w5.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_back_dam);
        w5.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_back_dam);
        w5.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_back_dam);
        wallHandler.copyAndAddObject(w5);*/
        //
        //particles:
        Mdouble meanRadius1 = 0.0019; // V2_1 -> 18.8075e-4; //17.3370e-4; //17.572e-4;  //35.14e-4/2;    //17.59e-4;//35.18e-4/2;
        Mdouble mode_log_normal = 34.674e-4/2.0;
        //Mdouble stdRadius1 = 1.4027; // 1.7824; //1.8549; //0.2 * meanRadius1;
        // log-normal distribution mean and std, calculated from data points via lognfit matlab function
        Mdouble LogmeanRadius = -6.2884; //-6.4266;//-6.3575;//log(meanRadius1);
        Mdouble LogstdRadius = 0.3384; //0.6178; //log(stdRadius1);
        Mdouble MinRadius = (13.183e-4/2.0); //10e-4; //diameter = 20//10microM = 10e-3 mm = 10e-4 cm
        Mdouble MaxRadius = (79.433e-4/2.0); //27.5e-4; // diameter = 55 // 27.5microM = 27.5e-3 mm = 27.5e-4cm
        //
        //the volume to be added:
        Mdouble addVolume = (90*scale*scale)*(((ToolWidth/scale)-10e-1)*scale)*(Gap); //0.0018 cm3
        // //set up random number generator
        //std::random_device rd;
        //std::mt19937 gen(rd());
        std::mt19937 gen;
        gen.seed(1);
        std::lognormal_distribution<> d(LogmeanRadius, LogstdRadius);
        //
        //add particles until the volume to be added is zero
        //logger(INFO,"Adding particles ...");
        SphericalParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(meanRadius1);
        Mdouble fillHeight = 0.0;
        while (addVolume>0) {
            Mdouble x = random.getRandomNumber(ScraperDistance+ScraperWidth, ToolLength-ToolThickness);
            Mdouble y = random.getRandomNumber(ToolThickness, (ToolWidth-ToolThickness));//(0.0, ToolWidth-(2.0*ToolThickness));//
            Mdouble z = random.getRandomNumber(0.0, fillHeight);
            p.setPosition({x, y, z});
            // check if particle can be inserted
            if (checkParticleForInteraction(p)) {
                particleHandler.copyAndAddObject(p);
                addVolume -= p.getVolume();
                do {
                    //p.setRadius(random_number(LogmeanRadius,LogstdRadius));
                    p.setRadius(d(gen));
                } while (p.getRadius()<MinRadius || p.getRadius()>MaxRadius);
                //(p.getRadius()<0.3505*meanRadius1 || p.getRadius()>2.112*meanRadius1);
                //(p.getRadius()<0.2*meanRadius || p.getRadius()>1.8*meanRadius); //reject too small or large radii
                //if (particleHandler.getNumberOfObjects()%100==0) std::cout << '.' << std::flush;
                //if (particleHandler.getNumberOfObjects()%1000==0) std::cout << ' ';
                //if (particleHandler.getNumberOfObjects()%10000==0) std::cout << addVolume << '\n';
            } else {
                fillHeight += 0.01*scale*scale*meanRadius1; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        //
        //
        //
        //logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());


    }
    //
    //Start tool motion
    //
    void actionsAfterTimeStep() override {
        //
        /*static bool wallRemoved = false;
        if (wallRemoved==false && getTime()>=TimeRemoveWall) {
            //wallHandler.getObject(7)->setPosition(Vec3D(0,0,4));
            //logger(INFO,"walls removed");
            //wallHandler.removeObject(5);
            //wallHandler.removeObject(4);
            wallRemoved = true;
        }*/
        //
        static bool moveTool = false;
        if (!moveTool && getTime()>= TimeMoveTool) {
            //logger(INFO,"Tool is Moving");
            wallHandler.getObject(1)->setVelocity(Vec3D(tool_Speed,0,0));
            wallHandler.getObject(2)->setVelocity(Vec3D(tool_Speed,0,0));
            wallHandler.getObject(3)->setVelocity(Vec3D(tool_Speed,0,0));
            wallHandler.getObject(4)->setVelocity(Vec3D(tool_Speed,0,0));
            moveTool = true;
        }
    }
    //
    /*void printTime() const override
    {
        logger(INFO,"t=%, tMax=%, N=%", getTime(),getTimeMax(), particleHandler.getSize());
    }*/

    //Another approach to perform the parametric study
    /*bool readNextArgument(int& i, int argc, char* argv[]) override
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
//
int main(int argc UNUSED, char *argv[] UNUSED)
{
    //logger(INFO,"Simple box for creating particles");
    Deposition deposition_problem;
    //
    Mdouble gravityValue = -981; // g=9.81 m/s2 - 981 cm/s2
    Mdouble XMaxDomain = (deposition_problem.scale)*(deposition_problem.scale)*(4*900e-1); // 1stScale = 0.9 cm - 2nd Scale = 0.225 cm
    Mdouble YMaxDomain = (deposition_problem.scale)*30e-1; // 1stScale = 0.4 cm - 2ndScale = 0.2 cm
    Mdouble ZMaxDomain = (deposition_problem.scale)*30e-1;//25e-1; // 1st Scale = 0.2 cm - 2ndScale = 0.1 cm
    //
    Mdouble MaxSimTime = 1.2;//3.0;
    //
    Mdouble MatDensity = 4.430;// density of Ti = 4.430g/cm3 - 4430kg/m3
    //
    Mdouble tc = 90e-7;//2e-4; //calculated as tc = 0.05-0.1 * tg : tg=sqrt(d50/g)
    Mdouble restitutionCoeff = 0.1;
    Mdouble timeStep = 0.02*tc;
    //
    //------------------------------------------------------------------------------Parametric friction:------------------------------------------------------------------------------
    deposition_problem.autoNumber();
    std::vector<int> studyNum=deposition_problem.get2DParametersFromRunNumber(5,5);
    //
    Mdouble SlidingFrictionCoeff = 1.0/studyNum[1]; //0.5;
    Mdouble RollingFrictionCoeff = 2.0/(studyNum[2]*10.0); //0.1;
    //Mdouble WSlidingFriCoeff = 1.0;
    //Mdouble WRollingFriCoeff = 0.2;
    //
    //std::cout << "studyNum2=" << studyNum[1] << std::endl;
    //std::cout << "Sliding Friction" << SlidingFrictionCoeff << std::endl;
    //std::cout << "Rolling Friction" << RollingFrictionCoeff << std::endl;
    //std::cout << "Wall Sliding Friction" << WSlidingFriCoeff << std::endl;
    //std::cout << "Wall Rolling Friction" << WRollingFriCoeff << std::endl;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    //
    deposition_problem.setName("SpreadingTiOldToolCase1");
    deposition_problem.setSystemDimensions(3);
    deposition_problem.setGravity(Vec3D(0.0,0.0,gravityValue));
    //deposition_problem.setXMin(0.0);
    //deposition_problem.setYMin(0.0);
    //deposition_problem.setZMin(0.0);
    deposition_problem.setXMax(XMaxDomain);
    deposition_problem.setYMax(YMaxDomain);
    deposition_problem.setZMax(ZMaxDomain);
    deposition_problem.setTimeMax(MaxSimTime);
    deposition_problem.setHGridMaxLevels(2);
    //
    auto species = deposition_problem.speciesHandler.copyAndAddObject(LinearViscoelasticFrictionSpecies());
    //
    species->setDensity(MatDensity);
    //species->setStiffness(10000);
    //species->setDissipation(100);
    //
    Mdouble mass = species->getMassFromRadius(13.183e-4/2); //13.183e-4/2 , Here we assign the mass of the smallest particle to calculate the stiffiness and disppation as we want to get ------
    species->setCollisionTimeAndRestitutionCoefficient(tc,restitutionCoeff,mass);
    //
    //------------------------------------------------------------------------------Friction coeffiecients Values---------------------------------------------------------------------------
    species->setSlidingStiffness(2./7.*species->getStiffness());
    species->setSlidingDissipation(2./7.*species->getDissipation());
    species->setSlidingFrictionCoefficient(SlidingFrictionCoeff);
    //species->setSlidingFrictionCoefficientStatic(0.6);
    //
    species->setRollingStiffness(2./5.*species->getStiffness());
    species->setRollingDissipation(2./5.*species->getDissipation());
    species->setRollingFrictionCoefficient(RollingFrictionCoeff);
    //
    //auto wallSpecies = deposition_problem.speciesHandler.copyAndAddObject(species);
    //auto wallParticleSpecies = deposition_problem.speciesHandler.getMixedObject(species,wallSpecies);
    //
    //wallParticleSpecies->setSlidingFrictionCoefficient(WSlidingFriCoeff);
    //wallParticleSpecies->setRollingFrictionCoefficient(WRollingFriCoeff);
    //------------------------------------------------------------------------------Data Output--------------------------------------------------------------------------------------------
    deposition_problem.setSaveCount(1000);//50
    deposition_problem.setParticlesWriteVTK(true);
    deposition_problem.setWallsWriteVTK(FileType::MULTIPLE_FILES);
    //
    deposition_problem.setTimeStep(timeStep);
    //deposition_problem.setTimeStep(1e-4);
    //------------------------------------------------------------------------------Parametric studies--------------------------------------------------------------------------------------------
    //deposition_problem.autoNumber();
    ////----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    deposition_problem.solve();//(argc, argv);
    return 0;
}


