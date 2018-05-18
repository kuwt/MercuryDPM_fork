//NEED TO FIND A SUITABLE VALUE OF KE BELOW WHICH WE RE-START THE EXCITATION OF PARTICLES!!
#include "Mercury3D.h"
#include "Particles/BaseParticle.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
#include<fstream>
using namespace std;

//**KIT** Attempting to set up an output file recording 
//only the maximal density points
//Currently declaring it **globally**, but may be able to
//find a better way to do this!
ofstream maxOut("maxDensitiesA-1.txt");
///#define DEBUG_OUTPUT

class my_problem : public Mercury3D{

public:
	
	void setupInitialConditions()
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
   
	//setting particle densities
    walls->setDensity(particle_density1);
    baseWallSpecies->setDensity(particle_density1);
    light->setDensity(particle_density1);
    
    //**KIT** hardwiring a friction coefficient for **BOTH SPECIES**
    //equal values, as we only typically use one species
    //in this code!
    light->setSlidingFrictionCoefficient(0.1);
    heavy->setSlidingFrictionCoefficient(0.1);
	//Setup properties of the light particles.
	

	light->setCollisionTimeAndRestitutionCoefficient(tc, r_particle11, particle_mass1);

	speciesHandler.getMixedObject(walls, light)->setCollisionTimeAndRestitutionCoefficient(tc, r_wall, 0.5, 2*particle_mass1);

	speciesHandler.getMixedObject(walls, light)->setSlidingFrictionCoefficient(mu_wall);
	
	//**KIT**setting up friction for interactions between differing particle species
	//speciesHandler.getMixedObject(light, heavy)->setSlidingFrictionCoefficient(0.1);

	//Setup properties of the heavy particles.
	heavy->setDensity(particle_density2);

	heavy->setCollisionTimeAndRestitutionCoefficient(tc, r_particle22, particle_mass2);

	speciesHandler.getMixedObject(heavy, walls)->setCollisionTimeAndRestitutionCoefficient(tc,r_wall,0.5, 2*particle_mass2); //note: here, tangential cor is set to 1 (i.e. no tangential friction). Can alternatively set to r_wall, for instance.

  speciesHandler.getMixedObject(heavy, walls)->setSlidingFrictionCoefficient\
  (mu_wall);

	speciesHandler.getMixedObject(heavy, light)->setCollisionTimeAndRestitutionCoefficient(tc,r_particle12, 2*particle_mass1*particle_mass2 / (particle_mass1 + particle_mass2));

	speciesHandler.getMixedObject(heavy, light)->setSlidingFrictionCoefficient(0.1);
	

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
        //**KIT** triggering a tap every time tapTime is reset to the current 
	//time, i.e. when the system's energy
	//NOTE, we can also **delay the taps**
	//by adding an additional value to (t - tapTime)
	if ( (t - tapTime) <= (1 / (4.0 * shaker_freq * constants::pi) ) ) {
	//**KIT** alternatively, if a fixed-period tap is required:
	//if (0 <= fmod(t,5) <= (1 / (4.0 * shaker_freq * constants::pi) ))        {
	  //*******************REMEMBER TO CHANGE T TO T - TAPTIME!!!!!!!**************
          return Vec3D(getXMin(),0.0,shaker_amp * 
		       std::sin( (t - tapTime) * 2.0 * shaker_freq * constants::pi));
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
        //initial test - tapping every 5 seconds
	if (0 <= fmod(t,5) <= (1 / (4.0 * shaker_freq * constants::pi) ))
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
        //initial test - tapping every 5 seconds
	if (0 <= fmod(t,5) <= (1 / (4.0 * shaker_freq * constants::pi) ))
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
        //initial test - tapping every 5 seconds
	if (0 <= fmod(5.0,t) <= (1 / (4.0 * shaker_freq * constants::pi) ))
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
    baseWall->setPrescribedPosition([this] (double time)
    {
        double t = time - 1.0;
	//**KIT** triggering a tap every time tapTime is reset to the current 
	//time, i.e. when the system's energy
	//NOTE, we can also **delay the taps**
	//by adding an additional value to (t - tapTime)
	if ( (t - tapTime) <= (1 / (4.0 * shaker_freq * constants::pi) ) ) {
        //initial test - tapping every 5 seconds
	//if (fmod(t,4.0) < ( 1 / (4.0 * shaker_freq * constants::pi) ) ){
	  return Vec3D(0.0,0.0,shaker_amp * std::sin( (t - tapTime) * 2.0 * shaker_freq * constants::pi));
        }
        else
        {
            return Vec3D(0.0,0.0,getZMin());
	 }

        
    });
    wallHandler.addObject(baseWall);


//Put the partilces on a grid with small random velocities
	BaseParticle p0;
	double max_radius = max(particle_radius1,particle_radius2);
		unsigned int N=numberOfParticles;
		double x=max_radius*1.01;
		double y=max_radius*1.01;
		double z=getZMin()+max_radius*1.01;
		for (int i=0;i<N;i++)
		{
			x += max_radius*2.01;
			if (x > getXMax()-max_radius) {
				x = max_radius*1.01;
				y += max_radius*2.01;
			}
			if (y > getYMax()-max_radius) {
				x = max_radius*1.01;
				y = max_radius*1.01;
				z += max_radius*2.01;
			}
			
			p0.setPosition(Vec3D(x,y,z));
			random.randomise(); //TURN THIS OFF IF I WANT REPRODUCIBLE DATA!!
			p0.setVelocity(Vec3D(random.getRandomNumber(-0.01,0.01),random.getRandomNumber(-0.01,0.01),random.getRandomNumber(-0.01,0.01)));
		
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

  void set_ParticleCOR(double c11,double c22,double c12){r_particle11=c11; r_particle22=c22; r_particle12=c12;}

  void set_ParticleCOR(double cor){r_particle11=cor; r_particle12=cor; r_particle22=cor;}

  void set_ParticleDensity(double density1, double density2){particle_density1=density1; particle_density2=density2;}

  void set_CollisionTime(double tc_in){tc=tc_in;}

  void set_NumberOfParticles(int np){numberOfParticles=np;}
    
  void set_Frequency(double f){shaker_freq=f;}

  void set_Amplitude(double a){shaker_amp=a;}

  void set_BaseCOR(double cor){r_base=cor;}

  void set_switch_plate_amplitude(double time_, double amp_){switch_time=time_; switch_amp=amp_;}

  //**KIT** defining a function to allow access to a simple variable which
  //can act as a trigger to initiate a tap
  void set_tapTime(double tTime){tapTime = tTime;}
  //**defining a simple boolean operator to uise as a "flag"
  //to ensure that taps are only triggered once per relaxation cycle!
  bool zeroFlag = false;
protected:

	void actionsBeforeTimeStep()
	{
		//After t=1.0 start to move the bottom wall
		double t=getTime();
		//**KIT* setting a "flag"to make sure that the onset
		//of ~zero energy only triggers a tap **ONCE**
		//bool zeroFlag = false; 
		//**KIT** Hardwiring a value of KE below which the
		//system is considered at rest, i.e. a new tap can 
		//be applied
		//a value **per particle** is chosen so as to, 
		//hopefully, provide a value applicable to all
		//systems!
		if ( (getKineticEnergy() / (particleHandler.getNumberOfObjects() ) < (pow(10,-14)) )
		     and (zeroFlag == false) ){
		  //if (zeroFlag == false) {
		    //**KIT**sets the 'zero flag' to true, to distinguish the FIRST point at which 
		    //energy drops to zero from all others!
		    zeroFlag = true;
		    //**KIT** Triggering a tap by setting tapTime to the current time
		    //i.e. zeroing time in the eyes of the vibration function.
		    set_tapTime(t);
		    //**KIT** Checking:
		    cout << "tapTime = " << tapTime << endl; 
		    //records the time corresponding to the instant before a 
		    //tap is applied - i.e. the local maximum in packing density
		    maxOut << t << "\t" << getGravitationalEnergy() << endl;
		  }
		//**KIT**checks!
		//cout << "t = " << t << endl;
		//cout << "KE = " 
		//<< getKineticEnergy() / particleHandler.getNumberOfObjects() 
		//<< endl;
     
		  //**KIT** resetting the zero flag to allow the next zero-energy
		  //point to be detected
		  if ((zeroFlag == true) and 
		      (getKineticEnergy() / particleHandler.getNumberOfObjects() 
		       > (pow(10,-10)) ) ) {
		    zeroFlag = false;
		  }
		  if (t>switch_time) {
		    set_Amplitude(switch_amp);
		  }
		  //cout << fmod(t,5) << endl;
		  //cout << "KE = " 
		  //<< getKineticEnergy() / particleHandler.getNumberOfObjects() 
		  //   << endl;
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
	double r_particle12;
	double r_particle22;
	double r_particle11;


	unsigned int numberOfParticles;
    
  double shaker_amp;
  double shaker_freq;
  //**KIT** defining a variable "tapTime" which acts to 
  //trigger a tap.  
  double tapTime = 0;
    
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
 	problem.setName("Atap-G2p75-f15-H20-W10-1");
 	problem.setTimeMax(30000);

	//Set Container Geometry
 	problem.setXMax(10./1000.);
 	//vector<int> study_num=problem.get_numbers(1,12);
	//problem.set_xmax(40.0+5.0*study_num[1]);
 	
 	
 	problem.setYMax(6./1000.);
 	problem.setZMax(300.0/1000.);


	//Set Partilce and Wall properties
    //
    //214 5mm glass particles
	problem.set_NumberOfParticles(20);
	problem.set_ParticleRadius(2.50001/1000 , 2.5/1000);
	problem.set_ParticleDensity(1200,1200); //2500 & 7850 --> glass & steel.
    	double tc=1e-4;
	problem.set_CollisionTime(tc);
	problem.set_WallCOR(0.9); //normally 0.7
	problem.set_WallFriction(0.1); //
	problem.set_ParticleCOR(0.9,0.9,0.9); //(11, 22, 12) OR just leave as a single value for equal COR.
	problem.set_BaseCOR(0.9);
	
	problem.set_switch_plate_amplitude(1.0,3.04/1000); //(time, new amplitude)
	
	
	//if (study_num[0] > 0)
	//		{
	//			std::cout << "Whole study has started" << std::endl;
	//			exit(0);
	//		}
	//	else
		//If the study is not complete save the data to disk and move on
	//		{
				//std::cout << "Going to launch a second/third/... code" <<std::endl;
				//problem.launch_new("Param");
	//		}
	
	
    
    
    //Shaker prop
//This is fake lowering the frequence and raising amplitude by a factor of 5.

    problem.set_Frequency(15);
    problem.set_Amplitude(3.04/1000); //normally 0.862
    //problem.set_Amplitude(2.5/1000);
    
    
	//Now run the code and solve  - time is set here because I am using the autodetect
    problem.setTimeStep(tc/50);
    problem.setSaveCount(1000*21);

    problem.solve();
}
