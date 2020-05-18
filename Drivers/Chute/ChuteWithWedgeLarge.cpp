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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include "scr/Chute.h"
#include "scr/Time.h"
//#include "scr/fstatistics.h"
using namespace std;

class ChutePeriodic : public Chute
{
public:

	void actionsBeforeTimeStep()
	{
	}
	
	void setupInitialConditions()
	{
		int N=get_N();
		add_particles();
		cout << "added " << get_N()-N << " particles at max height " << InflowHeight << endl;
		set_HGRID_num_buckets_to_power();
		write(cout,false);

		cout 
			<< "tc=" << getCollisionTime() 
			<< ", eps="	<< getRestitutionCoefficient()
			<< ", vmax=" << getMaximumVelocity()
			<< ", InflowHeight/getZMax()=" << getInflowHeight()/getZMax()
			<< endl << endl;
	}

	void add_particles() 
	{
		double Nmax = get_Nmax();
		//added after dinant update
		HGridActionsBeforeTimeLoop();
		HGridActionsBeforeTimeStep();
		InflowHeight = getZMax()/1.2;
		
		//try max_failed times to find new insertable particle
		while (Particles.size()<Nmax) {
			//cout << "Particle added" << Particles.size();
			create_inflow_particle();
			//added after dinant update
			if (!IsInsertable(P0)) InflowHeight += .00001* MaxInflowParticleRadius;
		};
	}

		
	void create_inflow_particle()
	{
		P0.Radius = MaxInflowParticleRadius;
		P0.computeMass(Species);
		P0.Position.X = random(getXMin()+P0.Radius,getChuteLength()-ContractionLength-P0.Radius);
		P0.Position.Y = random(getYMin()+2.0*P0.Radius,getYMax());
		P0.Position.Z = random(getZMin()+2.0*P0.Radius,getInflowHeight());
		P0.Velocity = Vec3D(0.0,0.0,0.0);
	}

	///sets parameters of particles and time stepping to the L3 type used in Silbert's papers
	void set_silbert_parameters() {
		//time stepping
		setTimeStep(1e-4*10);
		setTimeMax(2e7/10*getTimeStep());
		//particle properties
		setInflowParticleRadius(.5);
		setFixedParticleRadius(.5);
		setDensity(6/pi);
		setStiffness(2e5/100);
		setSlidingStiffness(2.0/7.0*getStiffness());
		setDissipation(25.0);
		setSlidingDissipation(0);
		setSlidingFrictionCoefficient(0.5);
	}

	///sets parameters of the chute
	void set_chute_parameters() {
		unsigned int NBottomX = 25; //*20
		unsigned int NBottomY = 10; //*10

		if (readDataFile("../ini_files/bottom0.5.ini")) {
			//chute properties
			setChuteLength(20*NBottomX);
			setChuteWidth(10*NBottomY);
			setChuteAngle(28.0, 1.0);
			set_H(13);
			
			unsigned int Nold = get_N();
			set_Nmax(NBottomY*NBottomX*Nold+getChuteLength()*getChuteWidth()*getZMax()/1.2);
			set_N(NBottomY*NBottomX*Nold);
			for (unsigned int k=0; k<NBottomY; k++)
			for (unsigned int j=0; j<NBottomX; j++)
			for (unsigned int i=0; i<Nold; i++) {
				unsigned int ind = (k*NBottomX+j)*Nold+i;
				Particles[ind]=Particles[i];
				Particles[ind].Position.X += j*20;
				Particles[ind].Position.Y += k*20;
				Particles[ind].fixParticle();
			}
			cout << "fixed particles: " << get_N()<< endl;
		} else {
			cerr << "Input data not found exiting " << endl;
			exit(-1);
		}

		cout << "!" << endl;
		//create a wall below the chute
		set_NWall(2);
		Walls[0].set(Vec3D( 0.0, 0.0, -1.0), 3.4*MaxInflowParticleRadius);
		
		// wedge at an angle (90 degrees blocks the chute)
		double relativeContractionWidth = 0.4;
		ContractionLength=150;
		double angle = atan(relativeContractionWidth*getChuteWidth()/2/ContractionLength);
		cout << "ContractionAngle: " << angle*180/pi << endl;
		Vec3D normalIntoWall = Vec3D(-1,0,0);
		Vec3D point = Vec3D(getChuteLength(),0,0);
		Walls[1].addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));
		normalIntoWall = Vec3D(sin(angle),cos(angle),0);
		point = Vec3D(getChuteLength()-ContractionLength,getChuteWidth()/2,0);
		Walls[1].addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));
		normalIntoWall = Vec3D(sin(angle),-cos(angle),0);
		Walls[1].addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));

		//create periodic walls
		set_NWallPeriodic(2);
		WallsPeriodic[0].set(Vec3D( 1.0, 0.0, 0.0), getXMin(), getXMax());
		WallsPeriodic[1].set(Vec3D( 0.0, 1.0, 0.0), getYMin(), getYMax());
		
	}
		
	double ContractionLength;

	//set approximate height of flow
	void set_H(double H) {setZMax(H*1.2);}

	void printTime() const {
		static Time2Finish timer;
		cout << "t=" << setprecision(3) << left << setw(6) << getTime() 
			<< ", tmax=" << setprecision(3) << left << setw(6) << getTimeMax()
			<< ", finish by " << setprecision(3) << left << setw(6) << timer.getFinishTime(getTime())
			<< endl;
	}

	
};

int main(int argc, char *argv[])
{
	ChutePeriodic problem;
	problem.setName("ChuteWithWedge");
	problem.set_silbert_parameters();
	problem.set_chute_parameters();
	problem.setSaveCount(500);
	//problem.setTimeMax(500*problem.getTimeStep());

	problem.readArguments(argc, argv);
	problem.solve();
}
