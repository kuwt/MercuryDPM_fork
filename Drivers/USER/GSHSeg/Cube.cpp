#include <iostream>
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <Walls/InfiniteWall.h>
#include <Mercury3D.h>
#include <cmath>
#include <Boundaries/PeriodicBoundary.h>
//Defining a class for a fully controlled cube which is energised in all spatial
//dimensions, walls can be heated and vibrated and gravity turned on as well.
// The idea of this code create a class that can used to test kinetic stress, gravity driven and vibrofluidised segregation models in issolcation.
class GSHCube: public Mercury3D {
//ask Ant what the point of making this public is! 
public:
  
 
  //declaring 3 integers to control the positioning
  //of particles within the system
  int xUp = 0, yUp = 0, zUp = 0;

  void setupInitialConditions (){

    ////////////////////////////////////////////////////////////////////////////
    //////////////////////DEFINING PARTICLE PROPERTIES//////////////////////////
    ////////////////////////////////////////////////////////////////////////////
	  
    //Computing the masses of particles such that they can then be used
    //to calculate relevant elastic and frictional properties
    double massLarge = (4 / 3) * constants::pi * pow(radiusL, 3.0) * rhoLarge;
    double massSmall = (4 / 3) * constants::pi * pow(radiusS, 3.0) * rhoSmall;
    //setting up an initial species of large particles and naming 
    //them appropriately!
    auto largeParticle = new LinearViscoelasticFrictionSpecies();
    largeParticle->setDensity(rhoLarge);
    largeParticle->setCollisionTimeAndRestitutionCoefficient(tc, largeCOR, massLarge);
      
    //largeParticle->setSlidingStiffness(2.0/7.0*largeParticle->getStiffness());
    //largeParticle->setSlidingDissipation(largeParticle->getDissipation());
    //largeParticle->setSlidingFrictionCoefficient(0.5);
      
    //DO LATER!!
    //largeParticle->setSlidingDissipation(speciesDrumWall->getDissipation());
    speciesHandler.addObject(largeParticle);
    //setting up a second species of smaller particle
    auto smallPartilce = new LinearViscoelasticFrictionSpecies();
    smallPartilce->setDensity(rhoSmall);
    smallPartilce->setCollisionTimeAndRestitutionCoefficient(tc, smallCOR, massSmall);
    //smallPartilce->setSlidingStiffness(2.0/7.0*smallPartilce->getStiffness());
    //smallPartilce->setSlidingDissipation(smallPartilce->getDissipation());
    //smallPartilce->setSlidingFrictionCoefficient(0.5);
    speciesHandler.addObject(smallPartilce);
    //DO LATER!!
    //largeParticle->setSlidingDissipation(speciesDrumWall->getDissipation());

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////DEFINING WALL PROPERTIES////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    //This is a wall with a coefficient of resitution of 1; and use
    auto insulatedWall = new LinearViscoelasticFrictionSpecies();
    //It will only be used for walls so teh density does not matter but it still a good idea to set one
    insulatedWall->setDensity(rhoLarge);
    insulatedWall->setCollisionTimeAndRestitutionCoefficient(tc, restWallCOR, massLarge);
    speciesHandler.addObject(insulatedWall);
    
    //Set the properties of the collisions between the particles and the walls
    auto largeParticleInsulatedWall = speciesHandler.getMixedObject(largeParticle, insulatedWall);
    auto smallParticleInsulatedWall = speciesHandler.getMixedObject(smallPartilce, insulatedWall);
    largeParticleInsulatedWall->setCollisionTimeAndRestitutionCoefficient(tc,restWallCOR,massLarge);
    smallParticleInsulatedWall->setCollisionTimeAndRestitutionCoefficient(tc,restWallCOR,massSmall);
    
      
    //This is a hot wall with a coefficient of resitution greater than 1.0
    auto hotWall = new LinearViscoelasticFrictionSpecies();
    hotWall->setDensity(rhoLarge);
    hotWall->setCollisionTimeAndRestitutionCoefficient(tc,leftWallCOR,massLarge);
    speciesHandler.addObject(hotWall);
      
    //Set the properties of the collisions between the particles and the walls
    auto largeParticleHotWall = speciesHandler.getMixedObject(largeParticle, hotWall);
    auto smallParticleHotWall = speciesHandler.getMixedObject(smallPartilce, hotWall);
    largeParticleHotWall->setCollisionTimeAndRestitutionCoefficient(tc,leftWallCOR,massLarge);
    smallParticleHotWall->setCollisionTimeAndRestitutionCoefficient(tc,leftWallCOR,massSmall);
      
      
    
    //This a cold wall with a coeffiecient of resitution less than 1.0
    auto coldWall = new LinearViscoelasticFrictionSpecies();
    coldWall->setDensity(rhoLarge);
    coldWall->setCollisionTimeAndRestitutionCoefficient(tc,rightWallCOR,massLarge);
    speciesHandler.addObject(coldWall);
      
    //Set the properties of the collisions between the particles and the walls
    auto largeParticleColdWall = speciesHandler.getMixedObject(largeParticle, coldWall);
    auto smallParticleColdWall = speciesHandler.getMixedObject(smallPartilce, coldWall);
    largeParticleColdWall->setCollisionTimeAndRestitutionCoefficient(tc,rightWallCOR,massLarge);
    smallParticleColdWall->setCollisionTimeAndRestitutionCoefficient(tc,rightWallCOR,massSmall);
 
      
      
    //Defining a moving wall bottom (base) wall
    //Declaring the wall as a 'new InfiniteWall'
    baseWall = new InfiniteWall;
      
    baseWall->setSpecies(insulatedWall);
      
    

    //setting the wall's orientation and position
    //set( (orientation) , (position) )
    baseWall->set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,getZMin()));
    //choosing the motion of the drum, i.e. prescribing its position at each
    //given point in time
    baseWall->setPrescribedPosition([this] (double time)
				    {
				      return Vec3D(0.0,0.0,getZMin()+(baseAmp * std::sin(time * 2 * constants::pi * baseFreq ) ) );
				    });
    wallHandler.addObject(baseWall);
      
    //similarly to the above, declaring the rest of the system's walls in such
    //a manner as they may all be prescribed differing motions
    //starting with the right (XMax) wall...
    rightWall = new InfiniteWall;
    rightWall->setSpecies(coldWall);
    rightWall->set(Vec3D(1.0,0.0,0.0),Vec3D(getXMax(),0,0));
    rightWall->setPrescribedPosition([this] (double time)
				    {
				      return Vec3D( (getXMax() + (rightAmp * std::sin(time * 2 * constants::pi * rightFreq ) ) ) ,0.0,0.0);
				    });
    wallHandler.addObject(rightWall);
    //...and left wall...
    leftWall = new InfiniteWall;
    leftWall->setSpecies(hotWall);
    leftWall->set(Vec3D(-1.0,0.0,0.0),Vec3D(getXMin(),0,0));
    leftWall->setPrescribedPosition([this] (double time)
    				    {
    				      return Vec3D( (getXMin() - (leftAmp * std::sin(time * 2 * constants::pi * leftFreq ) ) ) ,0.0,0.0);
    				    });
    wallHandler.addObject(leftWall);
    //...and the top wall, allowing a quasi-2D system to be 
    //symmetrically shaken.
    topWall = new InfiniteWall;
    topWall->setSpecies(insulatedWall);
    topWall->set(Vec3D(0.0,0.0,1.0),Vec3D(0.0,0.0,getZMax()));
    topWall->setPrescribedPosition([this] (double time)
    				    {
    				      return Vec3D(0.0,0.0,getZMax()-(topAmp * std::sin(time * 2 * constants::pi * topFreq ) ) );
    				    });
    wallHandler.addObject(topWall);
    //finally, making the front and back walls moveable to allow for 
    //a FULLY-3D symmetrically-shaken system!
    //back wall
    backWall = new InfiniteWall;
    backWall->setSpecies(insulatedWall);
    backWall->set(Vec3D(0.0,-1.0,0.0),Vec3D(0.0,getYMin(),0.0));
    backWall->setPrescribedPosition([this] (double time)
    				    {
    				      return Vec3D(0.0,getYMin()+(backAmp * std::sin(time * 2 * constants::pi * backFreq ) ),0.0 );
    				    });
    wallHandler.addObject(backWall);
    //front wall
    frontWall = new InfiniteWall;
    frontWall->setSpecies(insulatedWall);
    frontWall->set(Vec3D(0.0,1.0,0.0),Vec3D(0.0,getYMax(),0.0));
    frontWall->setPrescribedPosition([this] (double time)
    				    {
    				      return Vec3D(0.0,getYMax()-(frontAmp * std::sin(time * 2 * constants::pi * frontFreq ) ),0.0 );
    				    });
    wallHandler.addObject(frontWall);
      
    
    //*******************IF YOU WANT A PERIODIC BOUNDARY************************
    //PeriodicBoundary p1;
    //p1.set(Vec3D(1,0,0),getXMin(),getXMax());
    //boundaryHandler.copyAndAddObject(p1);
    
    BaseParticle p0;
    //declaring a 3D vector named 'position' to (unsurprisingly)
    //set up the initial positions of particles
    Vec3D position;
    /////////////////////////////////////////////////////////////////
    ///////////////ADDING PARTICLES//////////////////////////////////
    /////////////////////////////////////////////////////////////////
    //starting by adding the LARGER species to the system
    //this should (I hope!) make packing all particles in easier!
    for (int i=0; i < particleNumberLarge;i++)
    {
      p0.setRadius(radiusL);
      p0.setPosition(Vec3D(random.getRandomNumber(getXMin(),getXMax()),random.getRandomNumber(getYMin(),getYMax()),random.getRandomNumber(getZMin(),getZMax())));
      p0.setVelocity(Vec3D(random.getRandomNumber(-0.5,0.5),random.getRandomNumber(-0.5,0.5),random.getRandomNumber(-0.5,0.5)));

      p0.setVelocity(Vec3D(random.getRandomNumber(-0.5,0.5),random.getRandomNumber(-0.5,0.5),random.getRandomNumber(-0.5,0.5)));
      p0.setSpecies(largeParticle);
      particleHandler.copyAndAddObject(p0);
    }
    //once large particles have been added, we now start a new vertical row in which
    //to begin storing small particles. This is not the most efficient method, 
    //but a nice simple hack for now.
    //note that the choice to start large particles at the bottom is deliberate, to 
    //avoid any possible pre-segregation
    xUp = 0;
    yUp = 0;
    zUp++;
    zUp++;
    //inserting small particles above large!
    for (int i=0; i < particleNumberSmall;i++)
    {
        
      p0.setRadius(radiusS);
      p0.setPosition(Vec3D(random.getRandomNumber(getXMin(),getXMax()),random.getRandomNumber(getYMin(),getYMax()),random.getRandomNumber(getZMin(),getZMax())));
      p0.setVelocity(Vec3D(random.getRandomNumber(-0.5,0.5),random.getRandomNumber(-0.5,0.5),random.getRandomNumber(-0.5,0.5)));
      p0.setSpecies(smallPartilce);
      particleHandler.copyAndAddObject(p0);
    }
      
   
  }
  
  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////////SET FUNCTIONS/////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  
  //declaring and defining an access function to allow the radii of 
  //small particles...
  void setSmallParticleRadius(double radS) {
    radiusS = radS; 
  }
  //...and large particles! Note that I am using two separate functions
  //in order to make things nice and obvious for the user!
  void setLargeParticleRadius(double radL) {
    radiusL = radL; 
  }
  //similarly, declaring access functions to allow the user to choose a
  //number of small and large particles to insert into the system
  //note that these need only take an int as we cannot have fractions
  //of particles!
  void setSmallParticleNumber(int nS) {
    particleNumberSmall = nS;
  }
  void setLargeParticleNumber(int nL) {
    particleNumberLarge = nL;
  }
  //an additional pair of functions to allow user to easily set particles' 
  //restitution coefficients!
  void setCoefficientOfRestitutionLarge(double newCOR)
  {   
    largeCOR=newCOR;  
  }
  void setCoefficientOfRestitutionSmall(double newCOR)
  {        
    smallCOR=newCOR;
  }
  void setCoefficientOfRestitutionWall(double newCOR)
  {    
    wallCOR=newCOR;      
  }
  //yet another pair of functions to allow the **density** of the 
  //two particle species to be freely and easily changed
  void setDensityLarge(double density)
  {   
    rhoLarge = density;  
  }
  void setDensitySmall(double density)
  {   
    rhoSmall = density;  
  }
  //a function to allow control of the contact time, tc
  void setContactTime(double contactTime)
  {
    tc = contactTime;
  }
  //a pair of functions to allow the frequency and amplitude 
  //with which the base vibrates to be set by the user...
  void setBaseAmplitude(double amp) 
  {
    baseAmp = amp;
  }
  void setBaseFrequency(double freq) 
  {
    baseFreq = freq;
  }
  //...and a series of similar functions for all other walls!
  void setTopAmplitude(double amp) 
  {
    topAmp = amp;
  }
  void setTopFrequency(double freq) 
  {
    topFreq = freq;
  }
  void setLeftAmplitude(double amp) 
  {
    leftAmp = amp;
  }
  void setLeftFrequency(double freq) 
  {
    leftFreq = freq;
  }
void setRightAmplitude(double amp) 
  {
    rightAmp = amp;
  }
  void setRightFrequency(double freq) 
  {
    rightFreq = freq;
  }
  void setFrontAmplitude(double amp) 
  {
    frontAmp = amp;
  }
  void setFrontFrequency(double freq) 
  {
    frontFreq = freq;
  }
  void setBackAmplitude(double amp) 
  {
    backAmp = amp;
  }
  void setBackFrequency(double freq) 
  {
    backFreq = freq;
  }
  void setRightWallCOR(double cor)
  {
      rightWallCOR=cor;
  }
  void setLeftWallCOR(double cor)
  {
      leftWallCOR=cor;
  }
  void setRestWallCOR(double cor)
  {
      restWallCOR=cor;
  }
    
    


  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////////GET FUNCTIONS/////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  //a function to 'get' tc such that we can have an adaptive
  //time step
  double getContactTime() 
  {
    return tc;
  }
  
  ////////////////////////////////////////////////////////////////////////////
  ///////////////////////ADDITIONAL CLASS FUNCTIONS///////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  
  //Defining a function to superimpose a s secondary, long-time sine wave on 
  //the motion of each wall in order to allow the central, solid-like 'ball'
  //of particles formed under certain conditions to be moved and manipulated!
  //SHOULD HAVE SOME KIND OF SET MINIMUM IF I DON'T WANT IT HITTING
  //THE WALL!!

  //ask ant why I am setting this to private!
private:
  //declaring doubles to store the radii of large (l) and small (s) particles 
  //these are stored in the 'private'section of the code so that the may not 
  //accessed directly from the main code.
  double radiusS;
  double radiusL;
  //declaring **unsigned** ints to store the user-requested number
  //of each particle species (Small/Large).
  //the use of an unsigned int not only reduces memory requirements
  //but also prevents the use of negative numbers causing madness!
  unsigned int particleNumberSmall;
  unsigned int particleNumberLarge;
  //declaring values to record the restitution coefficients of
  //differing species
  double smallCOR;
  double largeCOR;
  double wallCOR;
  //declaring values of material density (rho) for the two
  //particle species used
  double rhoLarge;
  double rhoSmall;
  //The contact time, tc, for the current problem
  double tc;
  //Declaring a variable to store the frequency and amplitude 
  //with which the base vibrates, thus alllowing the details
  //and overall strength of shaking to be determined by 
  //the user
  double baseFreq;
  double baseAmp;
  //a similar series of variables for all the other system walls
  double topFreq;
  double topAmp;
  double leftFreq;
  double leftAmp;
  double rightFreq;
  double rightAmp;
  double frontFreq;
  double frontAmp;
  double backFreq;
  double backAmp;
    
  double leftWallCOR;
  double rightWallCOR;
  double restWallCOR;
    
  //declaring a base wall that may be made moveable
  InfiniteWall* baseWall;
  InfiniteWall* topWall;
  InfiniteWall* rightWall;
  InfiniteWall* leftWall;
  InfiniteWall* frontWall;
  InfiniteWall* backWall;
};



int main()
{
  //declaring an ***object(?)*** of the "GSHCube" class
  GSHCube problem;
  //giving the output file an appropriate name
  //(note that the '.' means that we are accessing a property 
  //belonging to an item of the above class and can then
  //assign it a value)
  problem.setName("GSHCube");
  //choosing the duration over which each
  //simulation runs
  problem.setTimeMax(51.0);
  
  ////////////////////////////////////////////////////////////////
  ///////////////////SETTING DRIVING STRENGTH/////////////////////
  ////////////////////////////////////////////////////////////////
  //setting the frequency and amplitude with which the system's 
  //base vibrates, allowing the user to determine both the 
  //overall strength of vibration as well as the specific details
  //of the driving
  //Can set amplitude to zero for a stationary wall.
  problem.setBaseAmplitude(0.0);
  //Note, it is generally advisable to keep f high and A low
  //such that the base vibrations can be simply treated as 
  //an energy source!
  problem.setBaseFrequency(0.0);
  //setting amplitude and frequency for system's top wall...
  problem.setTopAmplitude(0.0);
  problem.setTopFrequency(0.0);
  //...right-hand (XMax) wall...
  problem.setRightAmplitude(0.0);
  problem.setRightFrequency(0.0);
  //...left-hand (XMin) wall...
  problem.setLeftAmplitude(0.04);
  problem.setLeftFrequency(100.0);
  //...front wall...
  problem.setFrontAmplitude(0.0);
  problem.setFrontFrequency(0.0);
  //...and back wall.
  problem.setBackAmplitude(0.0);
  problem.setBackFrequency(0.0);
    
 
  problem.setLeftWallCOR(1.0);
  problem.setRightWallCOR(1.0);
  problem.setRestWallCOR(1.0);

  ////////////////////////////////////////////////////////////////
  ///////////////////SETTING PARTICLE DETAILS/////////////////////
  ////////////////////////////////////////////////////////////////
  //setting the **radii** of large and small particles for the
  //current 'problem'!
  //the values should be in SI units - m
  problem.setSmallParticleRadius(0.05);
  problem.setLargeParticleRadius(0.1);
  //setting the **number** of small and large particles to 
  //insert into the system
  problem.setSmallParticleNumber(3500);
  problem.setLargeParticleNumber(700);
  //setting the **coefficients of restitution** for all
  //species of particle and the walls
  //A value 1 is perfectly elastic, 
  //A value 0 is totally inelastic
  problem.setCoefficientOfRestitutionLarge(0.83);
  problem.setCoefficientOfRestitutionSmall(0.83);
  problem.setCoefficientOfRestitutionWall(0.83);
  //setting the **densities** of small and large particles
  //(do not need to worry about walls, as these have infinite
  //mass!)
  //The value should be in the SI units - kg/m^3
  problem.setDensityLarge(2000);
  problem.setDensitySmall(2000);
  
  ////////////////////////////////////////////////////////////////
  //////////////////DEFINING SYSTEM DIMENSIONS////////////////////
  ////////////////////////////////////////////////////////////////
  //defining the extent of the box in the x- y- and
  //z-dimensions (i.e. defining the maxima and minima
  //in these directions)
  problem.setXMax(2.5);
  problem.setYMax(2.5);
  problem.setZMax(2.5);
  problem.setXMin(-2.5);
  problem.setYMin(-2.5);
  problem.setZMin(-2.5);

  ////////////////////////////////////////////////////////////////
  ////////////////TYPICALLY CONSTANT PARAMETERS///////////////////
  ////////////////////////////////////////////////////////////////
  //setting tc, the contact time, which will strongly
  //influence the stiffness of the particles
  //a smaller tc means stiffer particles, but also
  //slower code!
  problem.setContactTime(1.0/200);
  //defining the strength of gravity in, respectively, 
  //the (x,y,z) directions (i.e. (0,0,-A) will 
  //give a gravitational force of strength 'A'
  //in the downward (negative z) direction
  //Note that, since I am using real values, g =9.81 (not 1!)
  problem.setGravity(Vec3D(0.0,0.0,0.0));
  //***ASK ANT FOR CLARIFICATION***
  problem.setSaveCount(100);
  //setting the duration of a time step
  problem.setTimeStep( ( problem.getContactTime() / 50.0 ) );
    
    
  problem.autoNumber();
    
  problem.solve();
    
  //Also write the paraview files in case I want to view it this way as well.
  problem.setParticlesWriteVTK(true);
  problem.setWallsWriteVTK(FileType::ONE_FILE);

  return 0;
}
