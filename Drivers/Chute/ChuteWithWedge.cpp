//Copyright (c) 2013-2018, The MercuryDPM Developers Team. All rights reserved.
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
	}

	void add_particles() 
	{
		double Nmax = get_N()+getXMax()*getYMax()*getZMax();
		set_Nmax(Nmax);
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
		P0.Position.X = random(getXMin()+2.0*P0.Radius,getXMax());
		P0.Position.Y = random(getYMin()+2.0*P0.Radius,getYMax());
		P0.Position.Z = random(getZMin()+2.0*P0.Radius,getInflowHeight());
		P0.Velocity = Vec3D(0.0,0.0,0.0);
	}

	///sets parameters of particles and time stepping to the L3 type used in Silbert's papers
	void set_silbert_parameters() {
		//time stepping
		setTimeStep(1e-4);
		setTimeMax(2e7*getTimeStep());
		//particle properties
		setInflowParticleRadius(.5);
		setFixedParticleRadius(.5);
		setDensity(6/pi);
		setStiffness(2e5);
		setSlidingStiffness(2.0/7.0*getStiffness());
		setDissipation(25.0);
		setSlidingDissipation(0);
		setSlidingFrictionCoefficient(0.5);
	}

	///sets parameters of the chute
	void set_chute_parameters() {
		unsigned int NFirstBottom = 2;
		vector<CParticle> FirstBottom;
		if (readDataFile("../ini_files/bottom0.5.ini")) {
			FirstBottom.resize(NFirstBottom*Particles.size());
			for (unsigned int j=0; j<NFirstBottom; j++)
			for (unsigned int i=0; i<Particles.size(); i++) {
				FirstBottom[j*Particles.size()+i]=Particles[i];
				FirstBottom[j*Particles.size()+i].Position.X += j*20;
			}
			setChuteLength(20*NFirstBottom);
		} else {
			cerr << "Input data not found exiting " << endl;
			exit(-1);
		}
		unsigned int NSecondBottom = 2;
		vector<CParticle> SecondBottom;
		if (readDataFile("../ini_files/bottom0.5.ini")) {
			SecondBottom.resize(NSecondBottom*Particles.size());
			for (unsigned int j=0; j<NSecondBottom; j++)
			for (unsigned int i=0; i<Particles.size(); i++) {
				SecondBottom[j*Particles.size()+i]=Particles[i];
				SecondBottom[j*Particles.size()+i].Position.X += (NFirstBottom+j)*20;
			}
			setChuteLength(20*(NFirstBottom+NSecondBottom));
		} else {
			cerr << "Input data not found exiting " << endl;
			exit(-1);
		}
		set_Nmax(FirstBottom.size()+SecondBottom.size()+getChuteLength()*getChuteWidth()*getZMax());
		Particles.resize(FirstBottom.size()+SecondBottom.size());
		for (unsigned int i=0; i<FirstBottom.size(); i++) {
			Particles[i]=FirstBottom[i];
			Particles[i].fixParticle();
		}
		for (unsigned int i=0; i<SecondBottom.size(); i++) {
			Particles[FirstBottom.size()+i]=SecondBottom[i];
			Particles[FirstBottom.size()+i].fixParticle();
		}

		//chute properties
		setChuteAngle(27.0, 1.0);
		setChuteWidth(10);
		set_H(10);
		//create a wall below the chute
		set_NWall(2);
		Walls[0].set(Vec3D( 0.0, 0.0, -1.0), 3.4*MaxInflowParticleRadius);
		
		// square 
		//~ Vec3D normalIntoWall = Vec3D(1,0,0);
		//~ Vec3D point = Vec3D(60,0,0);
		//~ Walls[1].addObject(normalIntoWall, Vec3D::Dot(normalIntoWall,point));
		//~ normalIntoWall = Vec3D(-1,0,0);
		//~ point = Vec3D(70,0,0);
		//~ Walls[1].addObject(normalIntoWall, Vec3D::Dot(normalIntoWall,point));
		//~ normalIntoWall = Vec3D(0,1,0);
		//~ point = Vec3D(0,2.5,0);
		//~ Walls[1].addObject(normalIntoWall, Vec3D::Dot(normalIntoWall,point));
		//~ normalIntoWall = Vec3D(0,-1,0);
		//~ point = Vec3D(0,7.5,0);
		//~ Walls[1].addObject(normalIntoWall, Vec3D::Dot(normalIntoWall,point));

		// wedge at an angle (90 degrees blocks the chute)
		double angle = 15*pi/180;
		Vec3D normalIntoWall = Vec3D(-1,0,0);
		Vec3D point = Vec3D(70,0,0);
		Walls[1].addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));
		normalIntoWall = Vec3D(sin(angle),cos(angle),0);
		point = Vec3D(60,5,0);
		Walls[1].addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));
		normalIntoWall = Vec3D(sin(angle),-cos(angle),0);
		point = Vec3D(60,5,0);
		Walls[1].addObject(normalIntoWall, Vec3D::dot(normalIntoWall,point));

		//create periodic walls
		set_NWallPeriodic(2);
		WallsPeriodic[0].set(Vec3D( 1.0, 0.0, 0.0), getXMin(), getXMax());
		WallsPeriodic[1].set(Vec3D( 0.0, 1.0, 0.0), getYMin(), getYMax());
		
	}

	//set approximate height of flow
	void set_H(double H) {setZMax(H*1.2);}

	//~ void printTime() const {
		//~ static Time2Finish timer;
		//~ cout << "t=" << setprecision(3) << left << setw(6) << getTime() 
			//~ << ", tmax=" << setprecision(3) << left << setw(6) << tmax
			//~ << ", N=" << setprecision(3) << left << setw(6) << Particles.size()
			//~ //<< ", time left=" << setprecision(3) << left << setw(6) << timer.getTime2Finish(t)
			//~ << ", finish by " << setprecision(3) << left << setw(6) << timer.getFinishTime(t)
			//~ << endl;
		//~ cout.flush();
	//~ }

	
};

int main(int argc, char *argv[])
{
	ChutePeriodic problem;
	problem.setName("ChuteWithWedge");
	problem.set_silbert_parameters();
	problem.set_chute_parameters();
	problem.setSaveCount(2000);
	//problem.setTimeMax(500*problem.getTimeStep());

	problem.readArguments(argc, argv);
	problem.solve();
}
