#include<iostream>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
//#define DEBUG_OUTPUT

class my_problem : public Mercury3D{

public:
	
	void setupInitialConditions() override
	{

	double particle_mass=4.0/3.0*constants::pi*pow(particle_radius,3);	

	//Six solid walls on all sides	

    auto species0 = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
    auto species1 = speciesHandler.copyAndAddObject(species0);
    auto species10 = speciesHandler.getMixedObject(species1, species0);

    InfiniteWall w0;
    w0.set(Vec3D(-1.0, 0.0, 0.0),getMin());
    wallHandler.copyAndAddObject(w0);

    w0.set(Vec3D(+1.0, 0.0, 0.0),getMax());
    wallHandler.copyAndAddObject(w0);

    w0.set(Vec3D( 0.0,-1.0, 0.0),getMin());
    wallHandler.copyAndAddObject(w0);

    w0.set(Vec3D( 0.0,+1.0, 0.0),getMax());
    wallHandler.copyAndAddObject(w0);

    w0.set(Vec3D( 0.0, 0.0,-1.0),getMin());
    w0.setPrescribedPosition([this] (double time)
        {
            double t = getTime()-1.0;
            if (t > 0.0)
            {
                return Vec3D(0.0,0.0,getZMin() + shaker_amp * std::sin(t * 2.0 * shaker_freq * constants::pi));
            }
            else
            {
                return Vec3D(0.0,0.0,getZMin());
            }
        });
     wallHandler.copyAndAddObject(w0);


//Put the partilces on a grid with small random velocities
		SphericalParticle p0;
		
		unsigned int N=numberOfParticles;
		double x=particle_radius*1.01;
		double y=particle_radius*1.01;
		double z=getZMin()+particle_radius*1.01;
		for (int i=0;i<N;i++)
		{
			x += particle_radius*2.01;
			if (x > getXMax()-particle_radius) {
				x = particle_radius*1.01;
				y += particle_radius*2.01;
			}
			if (y > getYMax()-particle_radius) {
				x = particle_radius*1.01;
				y = particle_radius*1.01;
				z += particle_radius*2.01;
			}
			
			p0.setPosition(Vec3D(x,y,z));
			p0.setVelocity(Vec3D(random.getRandomNumber(0.0,0.1),random.getRandomNumber(0.0,0.1),random.getRandomNumber(0.0,0.1)));
			p0.setRadius(particle_radius);
			p0.setSpecies(speciesHandler.getObject(1));
			speciesHandler.getObject(1)->setDensity(particle_density);
			p0.setMassForP3Statistics(0.163/1000);
			particleHandler.copyAndAddObject(p0);
		}
		setGravity(Vec3D(0.0,0.0,-9.81));


		species1->setCollisionTimeAndRestitutionCoefficient(tc,r_particle, particle_mass);
		species10->setCollisionTimeAndRestitutionCoefficient(tc,r_wall, 2*particle_mass);
	}


void set_ParticleRadius(double pr){particle_radius=pr;}

void set_WallCOR(double cor){r_wall=cor;}

void set_ParticleCOR(double cor){r_particle=cor;}

void set_ParticleDensity(double density){particle_density=density;}

void set_CollisionTime(double tc_in){tc=tc_in;}

void set_NumberOfParticles(int np){numberOfParticles=np;}
    
void set_Frequency(double f){shaker_freq=f;}

void set_Amplitude(double a){shaker_amp=a;}

protected:

	void actionsBeforeTimeStep() override
	{
    }


private:

	double particle_radius;
	double particle_density;
	double tc;

	double r_wall;
	double r_particle;

	unsigned int numberOfParticles;
    
    double shaker_amp;
    double shaker_freq;

};


int main(int argc UNUSED, char *argv[] UNUSED)
{

	
	//Set problem up
 	my_problem problem;
 	problem.setName("leidenfrost");
 	problem.setTimeMax(21);


	//Set Container Geometry
 	problem.setXMax(100.0/1000);
 	problem.setYMax(25.0/1000);
 	problem.setZMax(200.0/1000);


	//Set Partilce and Wall properties
    //
    //214 5mm glass particles
	problem.set_NumberOfParticles(214);
	problem.set_ParticleRadius(2.5/1000);
	problem.set_ParticleDensity(2500);
    	double tc=1e-5;
	problem.set_CollisionTime(tc);
	problem.set_WallCOR(.1);
	problem.set_ParticleCOR(0.91);
    
    
    //Shaker prop
    //This is fake lowering the frequence and raising amplitude by a factor of 5.

    problem.set_Frequency(70);
    problem.set_Amplitude(0.862/1000);
    //problem.set_Amplitude(2.5/1000);
    
    
	//Now run the code and solve  - time is set here because I am using the autodetect
    problem.autoNumber();
    problem.setTimeStep(tc/50);
    problem.setSaveCount(1000*21);

    problem.solve();
}
