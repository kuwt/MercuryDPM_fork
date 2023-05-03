#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include <stdlib.h>
#include<string>
#include<sstream>
using namespace std;

///#define DEBUG_OUTPUT

class my_problem : public Mercury3D{

public:

  void setupInitialConditions() override
  {

    double particle_mass1=4.0/3.0*particle_density1*constants::pi*pow(particle_radius1,3);

    double particle_mass2=4.0/3.0*particle_density2*constants::pi*pow(particle_radius2,3);


    auto walls = new LinearViscoelasticSlidingFrictionSpecies;
    speciesHandler.addObject(walls);

    auto light = new LinearViscoelasticSlidingFrictionSpecies;
    speciesHandler.addObject(light);

    auto heavy = new LinearViscoelasticSlidingFrictionSpecies;
    speciesHandler.addObject(heavy);

    auto baseWallSpecies = new LinearViscoelasticSlidingFrictionSpecies;
    speciesHandler.addObject(baseWallSpecies);

    walls->setDensity(particle_density1);
    baseWallSpecies->setDensity(particle_density1);
    light->setDensity(particle_density1);


    //Setup properties of the light particles.


    light->setCollisionTimeAndRestitutionCoefficient(tc, r_particle11, particle_mass1);

    light->setSlidingDissipation(light->getDissipation()*2./7.);
    light->setSlidingStiffness(light->getStiffness()*2./7.);
    light->setSlidingFrictionCoefficient(mu_p);

    speciesHandler.getMixedObject(walls, light)->setCollisionTimeAndRestitutionCoefficient(tc, r_wall, 0.5, 2*particle_mass1);

    speciesHandler.getMixedObject(walls, light)->setSlidingFrictionCoefficient(mu_wall);

    //Setup properties of the heavy particles.
    heavy->setDensity(particle_density2);

    heavy->setCollisionTimeAndRestitutionCoefficient(tc, r_particle22, particle_mass2);

    heavy->setSlidingDissipation(heavy->getDissipation()*2./7.);
    heavy->setSlidingStiffness(heavy->getStiffness()*2./7.);
    heavy->setSlidingFrictionCoefficient(mu_p);

    speciesHandler.getMixedObject(heavy, walls)->setCollisionTimeAndRestitutionCoefficient(tc,r_wall,0.5, 2*particle_mass2); //note: here, tangential cor is set to 1 (i.e. no tangential friction). Can alternatively set to r_wall, for instance.

    speciesHandler.getMixedObject(heavy, walls)->setSlidingFrictionCoefficient(mu_wall);

    speciesHandler.getMixedObject(heavy, light)->setCollisionTimeAndRestitutionCoefficient(tc,r_particle12, 2*particle_mass1*particle_mass2 / (particle_mass1 + particle_mass2));

    speciesHandler.getMixedObject(heavy, light)->setSlidingFrictionCoefficient(mu_p);
    speciesHandler.getMixedObject(heavy, light)->setSlidingDissipation(light->getDissipation()*2./7.);
    speciesHandler.getMixedObject(heavy, light)->setSlidingStiffness(light->getStiffness()*2./7.);

    //Setup properties of the baseWall particles.

    speciesHandler.getMixedObject(baseWallSpecies,light)->setCollisionTimeAndRestitutionCoefficient(tc,r_base,2*particle_mass1);
    speciesHandler.getMixedObject(baseWallSpecies,heavy)->setCollisionTimeAndRestitutionCoefficient(tc,r_base,2*particle_mass2);
    speciesHandler.getMixedObject(baseWallSpecies,light)->setSlidingFrictionCoefficient(mu_wall);
    speciesHandler.getMixedObject(baseWallSpecies,heavy)->setSlidingFrictionCoefficient(mu_wall);


    //Make it a 3D problem
    setDimension(3);


    //Five solid walls - open top


    leftWall = new InfiniteWall;
    rightWall = new InfiniteWall;
    frontWall = new InfiniteWall;
    backWall = new InfiniteWall;
    baseWall = new InfiniteWall;

    leftWall->set(Vec3D(-1.0, 0.0, 0.0),getMin());
    leftWall->setSpecies(walls);
    leftWall->setPrescribedPosition([this] (double time)
				    {
				      double t = time-1.0;
				      if (t > 0.0)
					{
					  return Vec3D(getXMin(),0.0,shaker_amp * std::sin(t * 2.0 * shaker_freq * constants::pi));
					}
				      else
					{
					  return Vec3D(getXMin(),0.0,0.0);
					}
				    });
    wallHandler.addObject(leftWall);


    rightWall->set(Vec3D(+1.0, 0.0, 0.0),getMax());
    rightWall->setSpecies(walls);
    rightWall->setPrescribedPosition([this] (double time)
				     {
				       double t = time-1.0;
				       if (t>0.0)
					 return Vec3D(+getXMax(),0.0,shaker_amp * std::sin(t * 2.0 * shaker_freq * constants::pi));
				       else
					 {
					   return Vec3D(getXMax(),0.0,0.0);
					 }
				     });
    wallHandler.addObject(rightWall);

    frontWall->set(Vec3D( 0.0,-1.0, 0.0),getMin());
    frontWall->setSpecies(walls);
    frontWall->setPrescribedPosition([this] (double time)
				     {
				       double t = time-1.0;
				       if (t>0.0)
					 return Vec3D(0.0,getYMin(),shaker_amp * std::sin(t * 2.0 * shaker_freq * constants::pi));
				       else
					 {
					   return Vec3D(0.0,getYMin(),0.0);
					 }
				     });
    wallHandler.addObject(frontWall);

    backWall->set(Vec3D( 0.0,+1.0, 0.0),getMax());
    backWall->setSpecies(walls);
    backWall->setPrescribedPosition([this] (double time)
				    {
				      double t = time - 1.0;
				      if (t>0.0)
					{
					  return Vec3D(0.0,getYMax(),shaker_amp * std::sin(t * 2.0 * shaker_freq * constants::pi));
					}
				      else
					{
					  return Vec3D(0.0,getYMax(),0.0);
					}
				    });
    wallHandler.addObject(backWall);

    baseWall->set(Vec3D( 0.0, 0.0,-1.0),getMin());
    baseWall->setSpecies(baseWallSpecies);
    //TAP - SETTING THE VIBRATION TO GIVE A SINGLE, HALF-SINE TAP!
    baseWall->setPrescribedPosition([this] (double time)
				    {
				      //delaying the tap to allow system to initially settle
				      double t = time - 3.0;
				      //applying a single tap of duration 0.5 / f, i.e. corresponding
				      //to one half sine wave at the desired frequency!
				      if (t>0.0 && t <= (0.5 / shaker_freq) )
					{
					  return Vec3D(0.0,0.0,shaker_amp * std::sin(t * 2.0 * shaker_freq * constants::pi));
					}
				      else
					{
					  return Vec3D(0.0,0.0,getZMin());
					}

				    });
    wallHandler.addObject(baseWall);


    //Put the particles on a grid with small random velocities
    //TAP initialising a random number generator to give random initial positions
    random.randomise();
    SphericalParticle p0;
    double max_radius = max(particle_radius1,particle_radius2);
    unsigned int N=numberOfParticles;
    double x=max_radius*1.0000000;
    double y=max_radius*1.0000000;
    double z=getZMin()+max_radius*1.0000000;
    for (int i=0;i<N;i++)
      {
	//setting particles at slightly offset positions
	//to give a random initial packing
	cout << "adding particle" << endl;
	x += max_radius*2.0000000*random.getRandomNumber(1,1.05);
	if (x >= (getXMax()-max_radius)) {
	  cout << "adding row" << endl;
	  x = max_radius*1.0000000;
	  y += max_radius*2.0000000;
	}
	//TAP SETTING A SPACE BETWEEN LAYERS
	//SO THAT THE INITIAL "DROP" WILL RANDOMISE
	//TO AN EXTENT
	if (y > getYMax()-max_radius) {
	  x = max_radius*1.0000000;
	  y = max_radius*1.0000000;
	  z += max_radius*3.0000000;
	}

	p0.setPosition(Vec3D(x,y,z));
	//TAP - SETTING RANDOM INITIAL VELOCITIES TO CREATE
	//A DEGREE OF DISORDER!
	random.randomise(); //TURN THIS OFF IF I WANT REPRODUCIBLE DATA!!
	p0.setVelocity(Vec3D(random.getRandomNumber(-0.01,0.01),random.getRandomNumber(-0.01,0.01),random.getRandomNumber(-0.01,0.01)));
	//TAP - ZERO INITIAL VELOCITY
	//p0.setVelocity(Vec3D(0,0,0));
	//species one for even numbers, species 2 for odd --> should give me an initially well-mixed system!
	if (i%2 == 0) {
	  p0.setRadius(particle_radius1);
	  p0.setSpecies(light);
	} else {
	  p0.setRadius(particle_radius2);
	  p0.setSpecies(heavy);
	}
	particleHandler.copyAndAddObject(p0);
      }
    setGravity(Vec3D(0.0,0.0,-9.81));




  }


  void set_ParticleRadius(double pr1, double pr2){particle_radius1=pr1;particle_radius2=pr2;}

  void set_WallCOR(double cor){r_wall=cor;}

  void set_WallFriction(double mu){mu_wall=mu;}

  void set_ParticleFriction(double mu){mu_p=mu;}

  void set_ParticleCOR(double c11,double c22,double c12){r_particle11=c11; r_particle22=c22; r_particle12=c12;}

  void set_ParticleCOR(double cor){r_particle11=cor; r_particle12=cor; r_particle22=cor;}

  void set_ParticleDensity(double density1, double density2){particle_density1=density1; particle_density2=density2;}

  void set_CollisionTime(double tc_in){tc=tc_in;}

  void set_NumberOfParticles(int np){numberOfParticles=np;}

  void set_Frequency(double f){shaker_freq=f;}

  void set_Amplitude(double a){shaker_amp=a;}

  void set_BaseCOR(double cor){r_base=cor;}

  void set_switch_plate_amplitude(double time_, double amp_){switch_time=time_; switch_amp=amp_;}

protected:

  void actionsBeforeTimeStep() override
  {
    //After t=1.0 start to move the bottom wall
    double t=getTime();

    if (t>switch_time) set_Amplitude(switch_amp);

  }


private:

  double particle_radius1;
  double particle_radius2;
  double particle_density1;
  double particle_density2;
  double tc;

  double r_wall;
  double r_base;
  double mu_wall;
  double mu_p;
  double r_particle12;
  double r_particle22;
  double r_particle11;


  unsigned int numberOfParticles;

  double shaker_amp;
  double shaker_freq;

  double switch_time;
  double switch_amp;

  InfiniteWall* leftWall;
  InfiniteWall* rightWall;
  InfiniteWall* backWall;
  InfiniteWall* frontWall;
  InfiniteWall* baseWall;


};


int main(int argc UNUSED, char *argv[] UNUSED)
{


  //Set problem up
  my_problem problem;
  //defining the particle number as a local variable so it can also be output to file
  int NP = 400*10;
  problem.set_NumberOfParticles(NP);
  //defining the CoR for ALL PARTICLES
  double corAll = 0.99;
  //defining the CoR for particle-wall interactions
  double corWall = 1.00;
  //defining the friction coefficient for particle-particle interactions
  double pFric = 1.0;
  //setting the vibration amplitude in a manner that can also be output as a filename
  double Amplitude = 20./1000;
  problem.set_ParticleCOR(corAll); //(11, 22, 12) OR just leave as a single value for equal COR.
  //creating a name that automatically gives key information regarding file
  //setting up the system's x and y dimensions
  double xSize = 100.0001;
  double ySize = 100.0001;
  //implicitly converting to ints for simple output in filename
  int xInt = xSize;
  int yInt = ySize;
  stringstream nameStream;
  std::string base = "fricinTap";
  nameStream << base << "-A" << Amplitude << "-e" << corAll << "-eW" << corWall
	     << "-mu_p" << pFric << "-size" << xInt << "x" << yInt << "-n" << NP ;
  problem.setName(nameStream.str());
  problem.setTimeMax(15.00);


  //Set Container Geometry
  problem.setXMax(xSize/1000.);
  problem.setYMax(ySize/1000.);
  problem.setZMax(250.0/1000.);


  //creating an output filestream through which the particle number can be retrieved.
  ofstream nOut;
  //re-initialising stringstream
  nameStream.str(std::string());
  //creating a relevant file name to output particle number
  nameStream << "particleNumber" <<  "-e" << corAll << "-size" << xInt << "x" << yInt
	     << "-n" << NP << ".txt";
  nOut.open(nameStream.str());
  //simply outputting the particle number to file so it can be read in by analysis programme
  nOut << NP;

  //Set Particle and Wall properties
  problem.set_ParticleRadius(2.5/1000 , 2.5/1000);
  problem.set_ParticleDensity(7850,7850); //2500 & 7850 --> glass & steel.
  double tc=1e-4;
  problem.set_CollisionTime(tc);
  problem.set_WallCOR(corWall); //normally 0.7
  problem.set_WallFriction(pFric); //
  problem.set_ParticleFriction(pFric); //

  problem.set_BaseCOR(corAll);

  problem.set_switch_plate_amplitude(0.0,Amplitude); //(time, new amplitude)


  problem.set_Frequency(30.0);
  problem.set_Amplitude(Amplitude); //normally 0.862
  //problem.set_Amplitude(2.5/1000);


  //Now run the code and solve  - time is set here because I am using the autodetect
  problem.setTimeStep(tc/50);
  problem.setSaveCount(1000);



  problem.solve();
}


