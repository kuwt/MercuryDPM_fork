#include <iostream>
#include "Mercury3D.h"
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include "Walls/IntersectionOfWalls.h"
//#include "Boundaries/CubeInsertionBoundary.h"
#include "Boundaries/PeriodicBoundary.h"
//#define DEBUG_OUTPUT

class my_problem : public Mercury3D{

public:
    
    my_problem()
    {
        
        //This all the default files in case you forgot to set, something. \todo at some point remove the defaults and flag you have not set a parameter.
        plateLength_=0.2;
        plateAngle_=0.0;
        particleDensity_=2000.0;
        tc=1e-5;
        particleRadius_=2.5/1000;
        rWall=1.0;
        rParticle=1.0;
        shakerAmp_=0.0;
        shakerFreq_=0.0;
        
        plateAngle_=0.0;
        motorAngle_=0.0;
        
        feederHeight_= 0.05;
        feederWidth_ = 0.0001;
        feederWedgeAngle_ = 35.0;
        feederDepth_ = 0.03;
        
        
        numberOfParticlesInserted_=0;
        particlesPerHour_=0;
        
    }
	
	void setupInitialConditions()
	{
        
        // Set up species (forces laws) for particles and walls
        auto speciesWall = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
        auto speciesParticles = speciesHandler.copyAndAddObject(speciesWall);
        auto speciesParticleWall = speciesHandler.getMixedObject(speciesParticles, speciesWall);
        
        
        //Set the properties of particles and walls.
        double particle_mass=4.0/3.0*constants::pi*pow(particleRadius_,3) * particleDensity_;
        speciesParticles->setDensity(particleDensity_);
        speciesParticles->setCollisionTimeAndRestitutionCoefficient(tc,rParticle, particle_mass); 
        speciesParticleWall->setCollisionTimeAndRestitutionCoefficient(tc,rWall, 2*particle_mass);
   
        //adjust for the feeder channel
        setXMin(- feederWidth_ - feederDepth_);

        //Now add thw alls
        //First solid wall at the left
        //@dducks, moved this to the far end of the feeder now
        InfiniteWall w0;
        w0.set(Vec3D(-1.0,0.0,0.0),Vec3D(getXMin(),0,0));
        wallHandler.copyAndAddObject(w0);
        
        //Construct the wedge. feederWedgeAngle is in degrees.
        Vec3D tmpNormal = {-1.0*std::sin(feederWedgeAngle_ * constants::pi / 180.0),
                           -1.0*std::cos(feederWedgeAngle_ * constants::pi / 180.0),
                           0.0};
        w0.set(tmpNormal,{0,0,0});
        //Add the wedge
        wallHandler.copyAndAddObject(w0);
        
        //Construct the middle wall of the feeder.
        IntersectionOfWalls w1;
        w1.addObject({-1.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
        w1.addObject({0.0, 1.0, 0.0}, {0.0, feederHeight_, 0.0});
        w1.addObject({1.0, 0.0, 0.0}, {-feederWidth_, 0.0, 0.0});
        wallHandler.copyAndAddObject(w1);
        
        //Now add the moving plate.
        w0.set(Vec3D(9.8*sin(plateAngle_),-9.8*cos(plateAngle_),0.0),Vec3D(0,0,0));
        w0.setPrescribedPosition([this] (double time)
            {
                double t = getTime();
                if (t > 0.0)
                {
                    return Vec3D(            shakerAmp_ * std::sin(t * 2.0 * shakerFreq_ * constants::pi)*std::sin(motorAngle_),
                                 getYMin() + shakerAmp_ * std::sin(t * 2.0 * shakerFreq_ * constants::pi)*std::cos(motorAngle_),
                                 0.0);
                }
                else
                {
                    return Vec3D(0.0,getYMin(),0.0);
                }
            });
        wallHandler.copyAndAddObject(w0);
        
        //Add an extra slopping wall to roll the particles in
        //w0.set(Vec3D(-1.0,-1.0,0.0),Vec3D(getYMin()+5.0*particleRadius_,getYMin()+5.0*particleRadius_,0);
        //wallHandler.copyAndAddObject(w0);
    
        //One periodic boundary in z
        PeriodicBoundary b0;
        b0.set(Vec3D(0.0,0.0,1.0),getZMin(),getZMax());
        boundaryHandler.copyAndAddObject(b0);
       
        //Gravity is now fixed to point down with strength 9.81
        setGravity(Vec3D(0.0,-9.81,0.0));

        //Now add in particles in the insertion region.
        addParticles();

        //Set the correct arguments for plotting
        setXBallsAdditionalArguments(" -3dturn 1 -solidf -noborder 4");    
        
//        setXMin(-feederWidth_);

        
    }


    /**
     * Sets the width of the wall of the feeder, measured from start to end
     * of the middle segment
     *        | |   |
     *        |_|   /
     *  -----______/
     *        | |
               W
     */
    void setFeederWallWidth(double width) { feederWidth_ = width; }
    
    
    /**
     * Sets the depth of the feed channel, measured from first to last wall
     *        | |   |
     *        |_|   /
     *  -----______/
     *          |    |
                   D
     */
    void setFeederDepth(double depth) { feederDepth_ = depth; }
    
    
    /*
     *  Sets the angle of the wedge. 0 means completely horizontal,
     *  90 means completely vertical. Beware of degrees.
     *
     *        | |   |
     *        |_|   /
     *  -----______/__ a (degrees)
     *        
              
     */
    void setFeederWedgeAngle(double angle) { feederWedgeAngle_ = angle; }
    
    /*
     *  Sets feeder height: the height of the entrance onto the moving plate
     *  
     *        | |   |
     *        |_|   /____
     *  -----______/_____ H
     *        
              
     */
    void setFeederHeight(double height) { feederHeight_ = height; }

    
    
    // This function insert the particles in a cube region on the left side of the domain
    void addParticles()
    {
        /// \todo remove the code below as it is a repeat, this is a short hack.
        double particle_mass=4.0/3.0*constants::pi*pow(particleRadius_,3);
        
        auto speciesWall = speciesHandler.getObject(0);
        auto speciesParticles = speciesHandler.getObject(1);
        
        BaseParticle p0;
		

        //Current insertion is in (xmin+pr, ymin+pr,zmin-pr) to (xmin+7pm,ymax-pr,zmax+pr) : recall peridoic in z but have insertion check bug for peroidic
        //@dducks: Nope. It is now inserted into the feeder.    
        unsigned int targetNumber=std::ceil(getTime()/(60.0*60.0)*particlesPerHour_);
        // THis is hack for an xballs bug. You need one particle to plot.
        if (targetNumber==0) targetNumber=1;
        unsigned int fail=0;
        while ((numberOfParticlesInserted_<targetNumber) && (fail<10))
        {
            
            
/*            p0.setPosition(Vec3D(random.getRandomNumber(getXMin()+particleRadius_,getXMin()+7.0*particleRadius_),
                                 random.getRandomNumber(getYMin()+1.0*particleRadius_,getYMin()+30.0*particleRadius_),
                                 random.getRandomNumber(getZMin()+particleRadius_,getZMax()-particleRadius_)
                                 ));
*/
            //We now have to insert particles into the feeder.
            //The feeder exists from (-width-depth, getYMin(), getZMin())
            //                    to (0,            getYMax(), getZMax())
            p0.setPosition({ random.getRandomNumber( -feederWidth_ - feederDepth_, -feederWidth_),
                             random.getRandomNumber( getYMin()                   , getYMax()),
                             random.getRandomNumber( getZMin()                   , getZMax())});
            p0.setVelocity(Vec3D(random.getRandomNumber(0.0,0.0),random.getRandomNumber(0.0,0.0),random.getRandomNumber(0.0,0.0)));
            p0.setRadius(particleRadius_);
            p0.setSpecies(speciesParticles);
            //p0.setMass(0.163/1000);
            if (checkParticleForInteraction(p0))
            {
                particleHandler.copyAndAddObject(p0);
                numberOfParticlesInserted_++;
                //std::cout << "Inserting particle" << "Number inserted" <<numberOfParticlesInserted_<< " Target number" << targetNumber << std::endl;
                fail=0;
            }
            else
            {
                fail++;
            }
        }
        
        if (fail>10)
        {
            std::cout << "WARNING :: Request particle insert rate not possible" << std::endl;
        }
        
        
    
        
        // Insertion boundary version, will be used in the future but currently not working.
        //        BaseParticle* insertionBoundaryParticle =new BaseParticle;
        //        insertionBoundaryParticle->setSpecies(speciesParticles);
        //
        //        CubeInsertionBoundary* insertionBoundary;
        //        insertionBoundary = boundaryHandler.copyAndAddObject(CubeInsertionBoundary());
        //
        //        insertionBoundary->set(insertionBoundaryParticle,1,Vec3D(particleRadius_,particleRadius_+0.02,getZMin()+particleRadius_*1.01),Vec3D(getXMax()/10.0-particleRadius_,getYMax()-particleRadius_,getZMax()-particleRadius_),Vec3D(-0.1,-0.1,-0.1),Vec3D(0.1,0.1,0.1),particleRadius_,particleRadius_);
    }
    
    void deleteParticles()
    {
        
        
        // delete all outflowing particles
        for (unsigned int i = 0; i < particleHandler.getNumberOfObjects();i++)
        {
            // If the particle is more than 5 particle dimeters off the end of plate delete.
            if (particleHandler.getObject(i)->getPosition().X > (getXMax()+10.0*particleRadius_) ) //||particleHandler.getObject(i)->Position.Z+particleHandler.getObject(i)->Radius<zMin_)
                
            {

                particleHandler.removeObject(i);
            }
        }
        
    }
    
 
        
    
    void setParticleRadius(double pr){particleRadius_=pr;}

    void setWallCOR(double cor){rWall=cor;}

    void setParticleCOR(double cor){rParticle=cor;}

    void setParticleDensity(double density){particleDensity_=density;}

    void setCollisionTime(double tc_in){tc=tc_in;}

    void setParticlesPerHour(int np){particlesPerHour_=np;}
    
    void setFrequency(double f){shakerFreq_=f;}

    void setAmplitude(double a){shakerAmp_=a;}
    
    void setMotorAngle(double angle){motorAngle_=angle/180.0*constants::pi;}
    
    void setPlateAngle(double angle)
        {
            plateAngle_=angle/180.0*constants::pi;
            setXMax(plateLength_*cos(plateAngle_));
        }
    
    void setPlateLength(double length)
        {
            plateLength_=length;
            setXMax(plateLength_*cos(plateAngle_));
        }

protected:

   void actionsBeforeTimeStep()
   {
        addParticles();
   }


   void actionsAfterTimeStep()
   {
        static int count;
        count++;
        if (count < 0) count=0;
        if (count > 1000)
        {
    
            deleteParticles();
            count = 0;
            
        }
    }


private:

    double particleRadius_;
    double particleDensity_;
    double tc;

    double rWall;
    double rParticle;
   
    double shakerAmp_;
    double shakerFreq_;
    
    double motorAngle_;
    double plateAngle_;
        
    double plateLength_;
    
    unsigned int numberOfParticlesInserted_;
    double particlesPerHour_;

    double feederWidth_;
    double feederDepth_;
    double feederHeight_;
    double feederWedgeAngle_;
    
};

int main(int argc, char *argv[]) {
    //Set problem up
    my_problem problem;
    problem.setName("ModelSystem");
    problem.setTimeMax(20);


    //Set Container Geometry (note this should be removed shortly, as these should be set automatically from plate and insertion properties.
    problem.setYMax(200.0/1000);
    problem.setZMax(50.0/1000);


    //Set Particle  properties
    //Note, this set the numberOfParticlesInsetedPerHour; however, note 0, inserts one particles onces (this is a special case).
    problem.setParticlesPerHour(1000000);
    problem.setParticleRadius(2.5/1000);
    problem.setParticleDensity(1400);
    double tc=1e-5;
    problem.setCollisionTime(tc);
    problem.setParticleCOR(0.5);
    
    //Set wall properties
    problem.setWallCOR(0.5);
    problem.setPlateAngle(2.83);
    problem.setPlateLength(7.0);
    problem.setFrequency(10);
    problem.setAmplitude(0.004);
    problem.setMotorAngle(20.0);
    

    //Set the feeder properties
    //Set the first wall to be very thin
    problem.setFeederWallWidth(0.1);
    //Set the depth of the channel to be 0.04
    problem.setFeederDepth(0.2);
    //Set the height of the entrance onto the plate to be 0.05
    problem.setFeederHeight(0.1);
    //Set the angle of the wedge to be 20 degrees
    problem.setFeederWedgeAngle(20.0);
   
    
	//Now run the code and solve  - time is set here because I am using the autodetec
    problem.setTimeStep(tc/50);
    problem.setSaveCount(1000*21);

	problem.autoNumber();
    problem.solve(argc,argv);
    
    return 0;
}
