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
  
  //  Shaker3d() : radius_s(0.5)
  //  {
    
  //  }
  //declaring 3 integers to control the positioning
  //of particles within the system
  int xUp = 0, yUp = 0, zUp = 0;

  void setupInitialConditions (){
	  
    //Computing the masses of particles such that they can then be used
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
    speciesHandler.copyAndAddObject(largeBalls);
    //setting up a second species of smaller particle
    auto smallBalls = new LinearViscoelasticFrictionSpecies();
    largeBalls->setDensity(rhoSmall);
    largeBalls->setCollisionTimeAndRestitutionCoefficient(tc, smallCOR, massSmall);
    //DO LATER!!
    //largeBalls->setSlidingDissipation(speciesDrumWall->getDissipation());


    InfiniteWall w1;
    w1.set(Vec3D(1.0,0.0,0.0),Vec3D(getXMax(),0,0));
    wallHandler.copyAndAddObject(w1);
    w1.set(Vec3D(-1.0,0.0,0.0),Vec3D(getXMin(),0,0));
    wallHandler.copyAndAddObject(w1);
    w1.set(Vec3D(0.0,1.0,0.0),Vec3D(0,getYMax(),0));
    wallHandler.copyAndAddObject(w1);
    w1.set(Vec3D(0.0,-1.0,0.0),Vec3D(0,getYMin(),0));
    wallHandler.copyAndAddObject(w1);
    w1.set(Vec3D(0.0,0.0,1.0),Vec3D(0,0,getZMax()));
    wallHandler.copyAndAddObject(w1);
    //w1.set(Vec3D(0.0,0.0,-1.0),Vec3D(0,0,getZMin()));
    //wallHandler.copyAndAddObject(w1);
    baseWall = new InfiniteWall;
    baseWall->set(Vec3D(0.0,0.0,-1.0),Vec3D(0.0,0.0,getZMin()));
    baseWall->setPrescribedPosition([this] (double time)
				    {
				      return Vec3D(0.0,0.0,getZMin()+(0.05 * std::sin(time * 2 * constants::pi * 25.0 ) ) );
				    });
    wallHandler.addObject(baseWall);
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
      else if ( (getYMin() + (2.1 * radiusL * (yUp + 1 ) ) ) < (getYMax() - (1.1 * radiusL) ) ) {
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
      p0.setVelocity(Vec3D(0.0,0.0,0.0));
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
      p0.setVelocity(Vec3D(0.0,0.0,0.0));
      p0.setSpecies(smallBalls);
      particleHandler.copyAndAddObject(p0);
    }
  }
  
  
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
  //a function to 'get' tc such that we can have an adaptive
  //timestep
  double getContactTime() 
  {
    return tc;
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
  //declaring a base wall that may be made moveable
  InfiniteWall* baseWall;

};



int main()
{
  //declaring an ***object(?)*** of the "Shaker3d" class
  Shaker3d problem;
  //giving the output file an appropriate name
  //(note that the '.' means that we are accessing a property 
  //belonging to an item of the above class and can then
  //assign it a value)
  problem.setName("balls");
  //setting the **radii** of large and small particles for the
  //current 'problem'!
  //the values should be in SI units - m
  problem.setSmallParticleRadius(0.05);
  problem.setLargeParticleRadius(0.1);
  //setting the **number** of small and large particles to 
  //insert into the system
  problem.setSmallParticleNumber(100);
  problem.setLargeParticleNumber(100);
  //setting the **coefficients of restitution** for all
  //species of particle and the walls
  //A value 1 is perfectly elastic, 
  //A value 0 is totally inelastic
  problem.setCoefficientOfRestitutionLarge(0.5);
  problem.setCoefficientOfRestitutionSmall(0.5);
  problem.setCoefficientOfRestitutionWall(0.5);
  //setting the **densities** of small and large particles
  //(do not need to worry about walls, as these have infinite
  //mass!)
  //The value should be in the SI units - kg/m^3
  problem.setDensityLarge(2000);
  problem.setDensitySmall(2000);
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
  problem.setGravity(Vec3D(0.0,0.0,-1.0));
  ///////////////////////////////////////////////////////
  /////////////////LESS IMPORTANT PARAMETERS/////////////
  ///////////////////////////////////////////////////////
  //defining the extent of the box in the x- y- and
  //z-dimensions (i.e. defining the maxima and minima
  //in these directions)
  problem.setXMax(1.0);
  problem.setYMax(1.0);
  problem.setZMax(1.0);
  problem.setXMin(-1.0);
  problem.setYMin(-1.0);
  problem.setZMin(-1.0);
  //problem.setTimeMax(50.0);
  //***ASK ANT FOR CLARIFICATION***
  problem.setSaveCount(100);
  //setting the duration of a timestep
  problem.setTimeStep( ( problem.getContactTime() / 50.0 ) );
  //choosing the duration over which each
  //simulation runs
  problem.setTimeMax(10.0);
  problem.solve();

  std::cout << "Hello Driver" << std::endl;
  return 0;
}
