#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
#include <iostream>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>
using namespace std;

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
   
    walls->setDensity(particle_density1);
    baseWallSpecies->setDensity(particle_density1);
    light->setDensity(particle_density1);
        

	//Setup properties of the light particles.
	

	light->setCollisionTimeAndRestitutionCoefficient(tc, r_particle11, particle_mass1);

	speciesHandler.getMixedObject(walls, light)->setCollisionTimeAndRestitutionCoefficient(tc, r_wall, 0.5, 2*particle_mass1);

	speciesHandler.getMixedObject(walls, light)->setSlidingFrictionCoefficient(mu_wall);

	//Setup properties of the heavy particles.
	heavy->setDensity(particle_density2);

	heavy->setCollisionTimeAndRestitutionCoefficient(tc, r_particle22, particle_mass2);

	speciesHandler.getMixedObject(heavy, walls)->setCollisionTimeAndRestitutionCoefficient(tc,r_wall,0.5, 2*particle_mass2); //note: here, tangential cor is set to 1 (i.e. no tangential friction). Can alternatively set to r_wall, for instance.

  speciesHandler.getMixedObject(heavy, walls)->setSlidingFrictionCoefficient\
  (mu_wall);

	speciesHandler.getMixedObject(heavy, light)->setCollisionTimeAndRestitutionCoefficient(tc,r_particle12, 2*particle_mass1*particle_mass2 / (particle_mass1 + particle_mass2));

	speciesHandler.getMixedObject(heavy, light)->setSlidingFrictionCoefficient(0.0);
	

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
    baseWall->setPrescribedPosition([this] (double time)
    {
        double t = time - 1.0;
        if (t>0.0)
        {
            return Vec3D(0.0,0.0,shaker_amp * std::sin(t * 2.0 * shaker_freq * constants::pi));
        }
        else
        {
            return Vec3D(0.0,0.0,getZMin());
        }
        
    });
    wallHandler.addObject(baseWall);


//Put the partilces on a grid with small random velocities
	SphericalParticle p0;
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

protected:

	void actionsBeforeTimeStep()
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
 	problem.setName("ThomasPics1");
 	problem.setTimeMax(100);


	//Set Container Geometry
 	problem.setXMax(100./1000.);
 	//vector<int> study_num=problem.get_numbers(1,12);
	//problem.set_xmax(40.0+5.0*study_num[1]);
 	
 	
 	problem.setYMax(100./1000.);
 	problem.setZMax(250.0/1000.);


	//Set Partilce and Wall properties
    //
    //214 5mm glass particles
	problem.set_NumberOfParticles(400);
	problem.set_ParticleRadius(2.50001/1000 , 2.5/1000);
	problem.set_ParticleDensity(7850,7850); //2500 & 7850 --> glass & steel.
    	double tc=1e-4;
	problem.set_CollisionTime(tc);
	problem.set_WallCOR(0.9); //normally 0.7
	problem.set_WallFriction(0.0); //
	problem.set_ParticleCOR(0.9,0.9,0.9); //(11, 22, 12) OR just leave as a single value for equal COR.
	problem.set_BaseCOR(0.6);
	
	problem.set_switch_plate_amplitude(1.0,1.55/1000); //(time, new amplitude)
	
	
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

    problem.set_Frequency(80.0);
    problem.set_Amplitude(1.55/1000); //normally 0.862
    //problem.set_Amplitude(2.5/1000);
    
    
	//Now run the code and solve  - time is set here because I am using the autodetect
    problem.setTimeStep(tc/50);
    problem.setSaveCount(1000*21);

    problem.solve();
}
