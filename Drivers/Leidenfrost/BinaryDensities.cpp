//Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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

#include<iostream>
#include <Species/LinearViscoelasticSlidingFrictionSpecies.h>

#include "Mercury3D.h"
#include "Walls/InfiniteWall.h"
//#define DEBUG_OUTPUT

class my_problem : public Mercury3D{

public:
	
	void setupInitialConditions()
	{


	

	double particle_mass=4.0/3.0*constants::pi*pow(particle_radius,3);	

	
	//Six solid walls on all sides	
    auto species0 = speciesHandler.copyAndAddObject(LinearViscoelasticSlidingFrictionSpecies());
    auto species1 = speciesHandler.copyAndAddObject(species0);
    auto species2 = speciesHandler.copyAndAddObject(species0);
    auto species10 = speciesHandler.getMixedObject(species1, species0);
    auto species20 = speciesHandler.getMixedObject(species2, species0);
    auto species21 = speciesHandler.getMixedObject(species2, species1);
    
    species1->setCollisionTimeAndRestitutionCoefficient(tc,r_particle, particle_mass);
    species10->setCollisionTimeAndRestitutionCoefficient(tc,r_wall, 2*particle_mass);
    species10->setSlidingFrictionCoefficient(mu_wall);
    species1->setDensity(particle_density1);
    
    species2->setCollisionTimeAndRestitutionCoefficient(tc,r_particle, particle_mass);
    species20->setCollisionTimeAndRestitutionCoefficient(tc,r_wall, 2*particle_mass);
    species21->setCollisionTimeAndRestitutionCoefficient(tc,r_wall, 2*particle_mass);
    species20->setSlidingFrictionCoefficient(mu_wall);
    species2->setDensity(particle_density2);
    
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
			p0.setVelocity(Vec3D(drand48()*0.01,drand48()*0.01,drand48()*0.01));
			if (i<N/2) {
				p0.setRadius(particle_radius);
				p0.setSpecies(speciesHandler.getObject(1));
			} else {
				p0.setRadius(particle_radius);
				p0.setSpecies(speciesHandler.getObject(2));
			}
			///\todo check whether setMass is needed here
			//p0.setMass(0.163/1000);
			particleHandler.copyAndAddObject(p0);
		}
		setGravity(Vec3D(0.0,0.0,-9.81));


		
	}


void set_ParticleRadius(double pr){particle_radius=pr;}

void set_WallCOR(double cor){r_wall=cor;}

void set_WallFriction(double mu){mu_wall=mu;}

void set_ParticleCOR(double cor){r_particle=cor;}

void set_ParticleDensity(double density1, double density2){particle_density1=density1; particle_density2=density2;}

void set_CollisionTime(double tc_in){tc=tc_in;}

void set_NumberOfParticles(int np){numberOfParticles=np;}
    
void set_Frequency(double f){shaker_freq=f;}

void set_Amplitude(double a){shaker_amp=a;}

///todo{DK: This is the old moving wall implementation, the new one has to be tested}
/*
protected:

	void actionsBeforeTimeStep()
	{
		//After t=1.0 start to move the bottom wall
		double t=getTime()-1.0;
		if (t>0.0)
		{
			wallHandler.getObject(4)->move(getZMin()-shaker_amp*sin(t*2.0*shaker_freq*constants::pi));
		}
	}
*/

private:

	double particle_radius;
	double particle_density1;
	double particle_density2;
	double tc;

	double r_wall;
	double mu_wall;
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
 	problem.setYMax(100.0/1000);
 	problem.setZMax(200.0/1000);


	//Set Partilce and Wall properties
    //
    //214 5mm glass particles
	problem.set_NumberOfParticles(1000);
	problem.set_ParticleRadius(2.5/1000);
	problem.set_ParticleDensity(2500,3000);
    	double tc=1e-5;
	problem.set_CollisionTime(tc);
	problem.set_WallCOR(.6);
	problem.set_WallFriction(4.0);
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
