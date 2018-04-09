//NEED TO MAKE A SMALL BALLS SPECIES!!!
#include <iostream>
#include <Species/LinearViscoelasticSpecies.h>
#include "Species/LinearViscoelasticFrictionSpecies.h"
#include <Walls/InfiniteWall.h>
#include <Mercury3D.h>
#include <cmath>
#include <Boundaries/PeriodicBoundary.h>
//Defining a class for a vibrofluidised bed which is energised in all spatial
//dimensions, as opposed to conventional vertical vibration
class Shaker3d: public Mercury3D {
//ask Ant what the point of making this public is! 
public:
  
 
  //declaring 3 integers to control the positioning
  //of particles within the system
  int xUp = 0, yUp = 0, zUp = 0;
  void setupInitialConditions () {

    ////////////////////////////////////////////////////////////////////////////
    //////////////////////DEFINING PARTICLE PROPERTIES//////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    //************************DEFINING INDIVIDUAL PARTICLES*********************
    //Computing the masses of particles such that they can then b
    //to calculate relevant elastic and frictional properties
    double massLarge = (4 / 3) * constants::pi * pow(radiusL, 3.0) * rhoLarge;
    double massSmall = (4 / 3) * constants::pi * pow(radiusS, 3.0) * rhoSmall;
    //setting up an initial species of large particles and naming 
    //them appropriately!
    auto largeBalls = new LinearViscoelasticFrictionSpecies();
    largeBalls->setDensity(rhoLarge);
    largeBalls->setCollisionTimeAndRestitutionCoefficient(tc, largeCOR, massLarge);
    //DO LATER!!
    //largeBalls->setSlidingDissipation(speciesDrumWall->getDissipation());
    //speciesHandler.copyAndAddObject(largeBalls);
    speciesHandler.addObject(largeBalls);
    //setting up a second species of smaller particle
    auto smallBalls = new LinearViscoelasticFrictionSpecies();
    smallBalls->setDensity(rhoSmall);
    smallBalls->setCollisionTimeAndRestitutionCoefficient(tc, smallCOR, massSmall);
    //speciesHandler.copyAndAddObject(smallBalls);
    speciesHandler.addObject(smallBalls);
    //DO LATER!!
    //largeBalls->setSlidingDissipation(speciesDrumWall->getDissipation());

    //************************DEFINING MIXED INTERACTIONS***********************
    //FIND OUT WHAT THIS DOES FROM ANT OR THOMAS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //defining the relevant parameters for collisions between unlike species
     LargeAndSmallBalls = speciesHandler.getMixedObject(largeBalls, smallBalls);
     LargeAndSmallBalls->setCollisionTimeAndRestitutionCoefficient(tc, (0.5*(largeCOR+smallCOR)), massLarge, massSmall);
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////DEFINING WALL PROPERTIES////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
 
    //Defining a moving wall bottom (base) wall
    //Declaring the wall as a 'new InfiniteWall'
    baseWall = new InfiniteWall;
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
    rightWall->set(Vec3D(1.0,0.0,0.0),Vec3D(getXMax(),0,0));
    rightWall->setPrescribedPosition([this] (double time)
				    {
				      return Vec3D( (getXMax() + (rightAmp * std::sin(time * 2 * constants::pi * rightFreq ) ) ) ,0.0,0.0);
				    });
    wallHandler.addObject(rightWall);
    //...and left wall...
    leftWall = new InfiniteWall;
    leftWall->set(Vec3D(-1.0,0.0,0.0),Vec3D(getXMin(),0,0));
    leftWall->setPrescribedPosition([this] (double time)
    				    {
    				      return Vec3D( (getXMin() - (leftAmp * std::sin(time * 2 * constants::pi * leftFreq ) ) ) ,0.0,0.0);
    				    });
    wallHandler.addObject(leftWall);
    //...and the top wall, allowing a quasi-2D system to be 
    //symmetrically shaken.
    topWall = new InfiniteWall;
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
    backWall->set(Vec3D(0.0,-1.0,0.0),Vec3D(0.0,getYMin(),0.0));
    backWall->setPrescribedPosition([this] (double time)
    				    {
    				      return Vec3D(0.0,getYMin()+(backAmp * std::sin(time * 2 * constants::pi * backFreq ) ),0.0 );
    				    });
    wallHandler.addObject(backWall);
    //front wall
    frontWall = new InfiniteWall;
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
    for (int i=0; i < particleNumberLarge;i++)	    {
      //checking to ensure that incrementing the x-position of
      //insertion will not take particles outside of the 
      //experimental volume and cause cock-ups
      if ( (getXMin() + (2.1 * radiusL * xUp) ) < (getXMax() - (1.1 * radiusL) ) ) {
	position.X = getXMin() + (1.1 * radiusL) + (2.1 * radiusL * xUp);
	position.Y = getYMin() + (1.1 * radiusL) + (2.1 * radiusL * yUp);
	position.Z = getZMin() + (1.1 * radiusL) + (2.1 * radiusL * zUp);
	//incrementing xUp means that the next particle will be placed slightly more
	//than a particle diameter away in the x direction
	xUp++;
      }
      //if the next increase in 'xUp' would place the next particle outside the 
      //experimental volume, xUp is reset to 0 and yUp is incremented by 1
      //(i.e. next particle will be placed at the start of a new row in the 
      //y (depth) direction) 
      else if ( (getYMin() + (2.1 * radiusL * (yUp + 1 ) ) ) < (getYMax() - (1.1 * radiusL) ) ) 
	{
	  xUp = 0;
	  //incrementing yUp starts a new row in the y (depth) direction
	  yUp++;
	  position.X = getXMin() + (1.1 * radiusL) + (2.1 * radiusL * xUp);
	  position.Y = getYMin() + (1.1 * radiusL) + (2.1 * radiusL * yUp);
	  position.Z = getZMin() + (1.1 * radiusL) + (2.1 * radiusL * zUp);
	  xUp++;
	}
      //finally, if an increase in yUp woulf place the next particle in an
      //impossible position, starts a new vertical row by incrementing zUp
      else {
	xUp = 0;
	yUp = 0;
	zUp++;
	position.X = getXMin() + (1.1 * radiusL) + (2.1 * radiusL * xUp);
	position.Y = getYMin() + (1.1 * radiusL) + (2.1 * radiusL * yUp);
	position.Z = getZMin() + (1.1 * radiusL) + (2.1 * radiusL * zUp);
	xUp++;
      }
      p0.setRadius(radiusL);
      p0.setPosition(position);
      //giving a small initial velocity with random orientation 
      p0.setVelocity(Vec3D(random.getRandomNumber(-0.01,0.01)
			   ,random.getRandomNumber(-0.01,0.01)
			   ,random.getRandomNumber(-0.01,0.01)
			   ));
      p0.setSpecies(largeBalls);
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
    for (int i=0; i < particleNumberSmall;i++)	    {
      //checking to ensure that incrementing the x-position of
      //insertion will not take particles outside of the 
      //experimental volume and cause cock-ups
      if ( (getXMin() + (2.1 * radiusS * xUp) ) < (getXMax() - (1.1 * radiusS) ) ) {
	position.X = getXMin() + (1.1 * radiusS) + (2.1 * radiusS * xUp);
	position.Y = getYMin() + (1.1 * radiusS) + (2.1 * radiusS * yUp);
	position.Z = getZMin() + (1.1 * radiusS) + (2.1 * radiusS * zUp);
	//incrementing xUp means that the next particle will be placed slightly more
	//than a particle diameter away in the x direction
	xUp++;
      }
      //if the next increase in 'xUp' would place the next particle outside the 
      //experimental volume, xUp is reset to 0 and yUp is incremented by 1
      //(i.e. next particle will be placed at the start of a new row in the 
      //y (depth) direction) 
      else if ( (getYMin() + (2.1 * radiusS * (yUp + 1 ) ) ) < (getYMax() - (1.1 * radiusS) ) ) {
	xUp = 0;
	//incrementing yUp starts a new row in the y (depth) direction
	yUp++;
	position.X = getXMin() + (1.1 * radiusS) + (2.1 * radiusS * xUp);
	position.Y = getYMin() + (1.1 * radiusS) + (2.1 * radiusS * yUp);
	position.Z = getZMin() + (1.1 * radiusS) + (2.1 * radiusS * zUp);
	xUp++;
      }
      //finally, if an increase in yUp would place the next particle in an
      //impossible position, starts a new vertical row by incrementing zUp
      else {
	xUp = 0;
	yUp = 0;
	zUp++;
	position.X = getXMin() + (1.1 * radiusS) + (2.1 * radiusS * xUp);
	position.Y = getYMin() + (1.1 * radiusS) + (2.1 * radiusS * yUp);
	position.Z = getZMin() + (1.1 * radiusS) + (2.1 * radiusS * zUp);
	xUp++;
      }   
      p0.setRadius(radiusS);
      p0.setPosition(position);
      //giving a small initial velocity with random orientation 
      p0.setVelocity(Vec3D(random.getRandomNumber(-0.01,0.01)
			   ,random.getRandomNumber(-0.01,0.01)
			   ,random.getRandomNumber(-0.01,0.01)
			   ));
      p0.setSpecies(smallBalls);
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
  //Providing the necessary set functions for the various parameters 
  //relating to the long-term sinusoidal variation in the driving
  //strength of each moving wall.
  void setLongPeriodBase(double period)
  {
    longPeriodBase = period;
  }
  void setLongPeriodTop(double period)
  {
    longPeriodTop = period;
  }
  void setLongPeriodRight(double period)
  {
    longPeriodRight = period;
  }
  void setLongPeriodLeft(double period)
  {
    longPeriodLeft = period;
  }
  void setLongPeriodFront(double period)
  {
    longPeriodFront = period;
  }
  void setLongPeriodBack(double period)
  {
    longPeriodBack = period;
  }
  void setMinAmpFracBase(double amp)
  {
    minAmpFracBase = amp;
  }
  void setMinAmpFracTop(double amp)
  {
    minAmpFracTop = amp;
  }
  void setMinAmpFracRight(double amp)
  {
    minAmpFracRight = amp;
  }
  void setMinAmpFracLeft(double amp)
  {
    minAmpFracLeft = amp;
  }
  void setMinAmpFracFront(double amp)
  {
    minAmpFracFront = amp;
  }
  void setMinAmpFracBack(double amp)
  {
    minAmpFracBack = amp;
  }
  void setOffsetBase(double off) 
  {
    offsetBase = off;
  }
  void setOffsetTop(double off) 
  {
    offsetTop = off;
  }
  void setOffsetRight(double off) 
  {
    offsetRight = off;
  }
  void setOffsetLeft(double off) 
  {
    offsetLeft = off;
  }
  void setOffsetFront(double off) 
  {
    offsetFront = off;
  }
  void setOffsetBack(double off) 
  {
    offsetBack = off;
  }



  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////////GET FUNCTIONS/////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  //a function to 'get' tc such that we can have an adaptive
  //timestep
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
  //declaring as a double type as the function will simply return, at each
  //time step, a scalar value to re-adjust the current driving strength
  double superSine(double xOffset, double yOffset, double period, double currentTime) {
    //making sure the offset has a reasonable value - should not be greater
    //than 1 as the maximum value output by this function should be 1 --
    //it modifies the basal amplitude, not re-writes it!
    if ( (yOffset > 1) or (yOffset < 0) ) {
      std::cout << "Error: Invalid value in argument of function 'superSine'" << std::endl;
      std::cout << "minAmp must be between 0 and 1" << std::endl;
      exit (EXIT_FAILURE);
    }
    //letting the user select a constant amplitude by inputting a value of 
    //zero for the function's period
    if (period == 0) {
      return 1.0;
    }
    else {
      //Setting the amplitude of the sine wave such that its maximum will be 
      //1 no matter what the offset -- i.e. using this function will not result
      //in an amplitude higher than the basal amplitude set in the main program
      double maxAmp = 1 - yOffset; 
      //defining the frequency of the superimposed wave based on the period input
      double f = (1 / period);
      //defining the output for a given time step!
      //note we are actually using cos (not a real sine wave!) as this gives us 
      //the 'normal' (basal) value of peak amplitude at t = 0;
      double longSine = yOffset + 
	( maxAmp * std::cos( (2 * constants::pi * f * currentTime) - (2 * constants::pi * f * xOffset) ) );
      return longSine;
    }
  }

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
  //declaring a base wall that may be made moveable
  InfiniteWall* baseWall;
  InfiniteWall* topWall;
  InfiniteWall* rightWall;
  InfiniteWall* leftWall;
  InfiniteWall* frontWall;
  InfiniteWall* backWall;
  //declaring the mixed species in order that collisions between
  //unlike species may be correctly modelled
  LinearViscoelasticFrictionMixedSpecies* LargeAndSmallBalls;
  //declaring a series of parameters to allow the strength
  //of vibration for each wall to be individually controlled
  //in a manner chosen by the user
  //For each wall, we can set:
  //the period of the sinusoidal variation in driving strength...
  double longPeriodBase;
  double longPeriodTop;
  double longPeriodRight;
  double longPeriodLeft;
  double longPeriodFront;
  double longPeriodBack;
  //...the minimum fraction of the basal amplitude to be reached
  //during the variation (i.e. allowing wall collisions to be
  //prevented, if desirable, or to limit motion of central 
  //particle agglomerate)...
  double minAmpFracBase;
  double minAmpFracTop;
  double minAmpFracRight;
  double minAmpFracLeft;
  double minAmpFracFront;
  double minAmpFracBack;
  //...and the offset of the sinusoidal variation, such that walls
  //can be moved out of phase etc. etc.
  double offsetBase;
  double offsetTop;
  double offsetRight;
  double offsetLeft;
  double offsetFront;
  double offsetBack;  
};



int main()
{
  //declaring an ***object(?)*** of the "Shaker3d" class
  Shaker3d problem;
  //giving the output file an appropriate name
  //(note that the '.' means that we are accessing a property 
  //belonging to an item of the above class and can then
  //assign it a value)
  problem.setName("npd1-A0.03-N2000");
  //choosing the duration over which each
  //simulation runs
  problem.setTimeMax(101.0);  
  ////////////////////////////////////////////////////////////////
  ///////////////////SETTING DRIVING STRENGTH/////////////////////
  ////////////////////////////////////////////////////////////////
  //setting the frequency and amplitude with which the system's 
  //base vibrates, allowing the user to determine both the 
  //overall strength of vibration as well as the specific details
  //of the driving
  //Can set amplitude to zero for a stationary wall.
  problem.setBaseAmplitude(0.03);
  //Note, it is generally advisable to keep f high and A low
  //such that the base vibrations can be simply treated as 
  //an energy source!
  problem.setBaseFrequency(200.0);
  //setting amplitude and frequency for system's top wall...
  problem.setTopAmplitude(0.03);
  problem.setTopFrequency(200.0);
  //...right-hand (XMax) wall...
  problem.setRightAmplitude(0.03);
  problem.setRightFrequency(200.0);
  //...left-hand (XMin) wall...
  problem.setLeftAmplitude(0.03);
  problem.setLeftFrequency(200.0);
  //...front wall...
  problem.setFrontAmplitude(0.00);
  problem.setFrontFrequency(200.0);
  //...and back wall.
  problem.setBackAmplitude(0.00);
  problem.setBackFrequency(200.0);
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
  problem.setSmallParticleNumber(2000);
  problem.setLargeParticleNumber(0);
  //setting the **coefficients of restitution** for all
  //species of particle and the walls
  //A value 1 is perfectly elastic, 
  //A value 0 is totally inelastic
  problem.setCoefficientOfRestitutionLarge(0.91);
  problem.setCoefficientOfRestitutionSmall(0.91);
  //CURRENTLY DOES NOTHING!
  problem.setCoefficientOfRestitutionWall(0.91);
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
  problem.setXMax(4.0);
  problem.setYMax(0.5);
  problem.setZMax(4.0);
  problem.setXMin(0.0);
  problem.setYMin(0.0);
  problem.setZMin(0.0);
  ////////////////////////////////////////////////////////////////
  ////////////////TYPICALLY CONSTANT PARAMETERS///////////////////
  ////////////////////////////////////////////////////////////////
  //setting tc, the contact time, which will strongly
  //influence the stiffness of the particles
  //a smaller tc means stiffer particles, but also
  //slower code!
  problem.setContactTime(1.0/800);
  //defining the strength of gravity in, respectively, 
  //the (x,y,z) directions (i.e. (0,0,-A) will 
  //give a gravitational force of strength 'A'
  //in the downward (negative z) direction
  //Note that, since I am using real values, g =9.81 (not 1!)
  problem.setGravity(Vec3D(0.0,0.0,-1.0));
  //***ASK ANT FOR CLARIFICATION***
  problem.setSaveCount(1000);
  //setting the duration of a timestep
  problem.setTimeStep( ( problem.getContactTime() / 50.0 ) );
  problem.solve();

  std::cout << "Hello Driver" << std::endl;
  return 0;
}


