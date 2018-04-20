//
// Created by yousef on 4/4/18.
//


#include "Mercury3D.h"
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
#include "Particles/BaseParticle.h"
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
    Mdouble TimeRemoveWall = 0.04;//0.2;
    Mdouble TimeMoveTool = 0.05;//0.25;
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
        Mdouble ToolLength = 45e-1*scale; // unit cm
        Mdouble ToolWidth = 30e-1*scale;
        Mdouble ToolHight = 20e-1*scale;
        Mdouble ToolThickness = 5e-1*scale;
        Mdouble Gap = 0.1e-1; //0.1e-1; // keeping the gap as the actual dim requires scaling down the track volume by e-2 instead of (0.5e-1) which corresponds to 9000 particles
        //Mdouble TrackLength = 900e-1*(scale*scale);
        Mdouble Angel = constants::pi/180.0*30;
        Mdouble AngelDim = scale*4e-1;
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
        Vec3D maxPoint1_tool = Vec3D(ToolLength,ToolThickness,ToolHight); //Vec3D(getXMax(),0.05*getYMax(),getZMax())
        //back wall
        Vec3D minPoint2_tool = Vec3D(0.0,0.0,Gap);
        Vec3D maxPoint2_tool = Vec3D(ToolThickness,(ToolWidth-ToolThickness),ToolHight); //Vec3D(0.2*getXMax(),getYMax(),getZMax())
        //side wall
        Vec3D minPoint3_tool = Vec3D(0.0,(ToolWidth-ToolThickness),0.0);
        Vec3D maxPoint3_tool = Vec3D(ToolLength,ToolWidth,ToolHight); // Vec3D(getXMax(),getYMax(),getZMax())
        //
        //Front Dam
        Vec3D minPoint_dam = Vec3D(ToolLength,0.0,0.0);
        Vec3D maxPoint_dam = Vec3D(damThickness,ToolWidth,ToolHight);
        //Back Dam
        Vec3D minPoint_back_dam = Vec3D(0.0,ToolThickness,0.0);
        Vec3D maxPoint_back_dam = Vec3D(ToolThickness,(ToolWidth-ToolThickness),Gap);
        //
        //Bottom Wall
        InfiniteWall w0;
        w0.setSpecies(speciesHandler.getObject(1));
        w0.set(Vec3D(0.0, 0.0, -1.0), Vec3D(0.0,0.0,0.0));
        wallHandler.copyAndAddObject(w0);
        //TOOL:
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
        //Back wall:
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
        //Back Dam
        IntersectionOfWalls w5;
        w5.setSpecies(speciesHandler.getObject(0));
        w5.addObject(Vec3D(1.0,0.0,0.0),minPoint_back_dam);
        w5.addObject(Vec3D(0.0,1.0,0.0),minPoint_back_dam);
        w5.addObject(Vec3D(0.0,0.0,1.0),minPoint_back_dam);
        w5.addObject(Vec3D(-1.0,0.0,0.0),maxPoint_back_dam);
        w5.addObject(Vec3D(0.0,-1.0,0.0),maxPoint_back_dam);
        w5.addObject(Vec3D(0.0,0.0,-1.0),maxPoint_back_dam);
        wallHandler.copyAndAddObject(w5);
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
        std::random_device rd;
        std::mt19937 gen(rd());
        std::lognormal_distribution<> d(LogmeanRadius, LogstdRadius);
        //
        //add particles until the volume to be added is zero
        //logger(INFO,"Adding particles ...");
        BaseParticle p;
        p.setSpecies(speciesHandler.getObject(0));
        p.setRadius(meanRadius1);
        Mdouble fillHeight = 0.0;
        while (addVolume>0) {
            Mdouble x = random.getRandomNumber(ToolThickness, ToolLength);
            Mdouble y = random.getRandomNumber(ToolThickness, (ToolWidth-ToolThickness));
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
                fillHeight += 0.01*meanRadius1; //increase fill height (slowly to insert particles as low as possible)
            }
        }
        //
        //
        //
        logger(INFO," Inserted % particles",particleHandler.getNumberOfObjects());


    }
    //
    //Start tool motion
    //
    void actionsAfterTimeStep() override {
        //
        static bool wallRemoved = false;
        if (wallRemoved==false && getTime()>=TimeRemoveWall) {
            //wallHandler.getObject(7)->setPosition(Vec3D(0,0,4));
            //logger(INFO,"walls removed");
            wallHandler.removeObject(5);
            wallHandler.removeObject(4);
            wallRemoved = true;
        }
        //
        static bool moveTool = false;
        if (!moveTool && getTime()>= TimeMoveTool) {
            //logger(INFO,"Tool is Moving");
            wallHandler.getObject(1)->setVelocity(Vec3D(tool_Speed,0,0));
            wallHandler.getObject(2)->setVelocity(Vec3D(tool_Speed,0,0));
            wallHandler.getObject(3)->setVelocity(Vec3D(tool_Speed,0,0));
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
    Mdouble XMaxDomain = (deposition_problem.scale)*(deposition_problem.scale)*900e-1; // 1stScale = 0.9 cm - 2nd Scale = 0.225 cm
    Mdouble YMaxDomain = (deposition_problem.scale)*40e-1; // 1stScale = 0.4 cm - 2ndScale = 0.2 cm
    Mdouble ZMaxDomain = (deposition_problem.scale)*20e-1; // 1st Scale = 0.2 cm - 2ndScale = 0.1 cm
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
    Mdouble SlidingFrictionCoeff = 1.0/studyNum[1];//0.5;
    Mdouble RollingFrictionCoeff = 2.0/(studyNum[2]*10.0);//0.1;
    Mdouble WSlidingFriCoeff = 1.0/studyNum[1];//0.5;
    Mdouble WRollingFriCoeff = 2.0/(studyNum[2]*10.0);//0.1;
    //
    //std::cout << "studyNum2=" << studyNum[1] << std::endl;
    //std::cout << "Sliding Friction" << SlidingFrictionCoeff << std::endl;
    //std::cout << "Rolling Friction" << RollingFrictionCoeff << std::endl;
    //std::cout << "Wall Sliding Friction" << WSlidingFriCoeff << std::endl;
    //std::cout << "Wall Rolling Friction" << WRollingFriCoeff << std::endl;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    //
    deposition_problem.setName("Deposition_Ti_Calibration_Reduced"); //Deposition_gaussian_parametric
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
    Mdouble mass = species->getMassFromRadius(13.183e-4/2); // 11.092 //10e-4 , 13.183e-4/2 // Here we assign the mass of the smallest particle to calculate the stiffiness and disppation as we want to get ------
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
    auto wallSpecies = deposition_problem.speciesHandler.copyAndAddObject(species);
    auto wallParticleSpecies = deposition_problem.speciesHandler.getMixedObject(species,wallSpecies);
    //
    wallParticleSpecies->setSlidingFrictionCoefficient(WSlidingFriCoeff);
    wallParticleSpecies->setRollingFrictionCoefficient(WRollingFriCoeff);
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

