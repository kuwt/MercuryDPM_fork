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
#include "scr/StatisticsVector.h"
using namespace std;

class ChutePeriodic : public Chute
{
public:

	void actionsBeforeTimeStep(){};
		
	void setupInitialConditions()
	{
		set_NWall(0);
		set_NWallPeriodic(2);
		WallsPeriodic[0].set(Vec3D( 1.0, 0.0, 0.0), getXMin(), xmax);
		WallsPeriodic[1].set(Vec3D( 0.0, 1.0, 0.0), getYMin(), getYMax());
		
		add_particles();
		write(cout,true);
	}

	void add_particles() 
	{
		set_Nmax(4);
		set_HGRID_num_buckets_to_power(get_Nmax());
		HGridActionsBeforeTimeLoop();
		HGridActionsBeforeTimeStep();

		Vec3D Gravity = getGravity();
		Vec3D Tangential = Vec3D(-Gravity.Z,0,Gravity.X);

		switch (4) {
		case 1: //normal forces
			P0.Radius = MaxInflowParticleRadius;
			P0.computeMass(Species);
			P0.Position = Vec3D(0,5,0);
			P0.fixParticle();
			add_particle(P0);
			P0.unfix(Species);

			P0.Radius = MaxInflowParticleRadius;
			P0.computeMass(Species);
			P0.Position += -(1.-2./getStiffness())*getGravity();
			P0.Velocity = Vec3D(0.0,0.0,0.0);
			cout << P0 << endl;
			add_particle(P0);

			P0.Position += -(1.-1./getStiffness())*getGravity();
			P0.Velocity = Vec3D(0.0,0.0,0.0);
			add_particle(P0);
			break;
		case 2: //normal and tangential forces
			P0.Radius = MaxInflowParticleRadius;
			P0.computeMass(Species);
			P0.Position = Vec3D(10,5,1);
			P0.Velocity = Vec3D(0.0,0.0,-1.0);
			cout << P0 << endl;
			add_particle(P0);

			P0.Radius = MaxInflowParticleRadius;
			P0.computeMass(Species);
			P0.Position -= (1.-2./getStiffness())*Tangential;
			P0.fixParticle();
			add_particle(P0);
			break;
		case 3: ///steady case
			P0.Radius = MaxInflowParticleRadius;
			P0.computeMass(Species);
			P0.Position = Vec3D(0,5,0);
			P0.fixParticle();
			add_particle(P0);
			P0.unfix(Species);

			P0.Radius = MaxInflowParticleRadius;
			P0.computeMass(Species);
			P0.Position += -(1.5-1./getStiffness())*getGravity();
			P0.Velocity = Vec3D(0.0,0.0,0.0);
			cout << P0 << endl;
			add_particle(P0);

			setStiffness(2e5);
			setDissipation(0.0);
			cout << getDissipation() << endl;
			setTimeStep(1e-4);
			setTimeMax(9.93);
			break;
		case 4: //normal forces
			setChuteAngle(0.0, 1.0);
			
			P0.Radius = MaxInflowParticleRadius;
			P0.computeMass(Species);
			P0.Position = Vec3D(0,5,0);
			P0.fixParticle();
			add_particle(P0);
			P0.unfix(Species);

			P0.Radius = MaxInflowParticleRadius;
			P0.computeMass(Species);
			P0.Position += -(1.-1./getStiffness())*getGravity();
			P0.Velocity = Vec3D(0.0,0.0,0.0);
			cout << P0 << endl;
			add_particle(P0);

			P0.Position = Vec3D(0,5,3);
			P0.Velocity = Vec3D(0.0,0.0,0.0);
			add_particle(P0);
			P0.Position += -(1.-1./getStiffness())*getGravity();
			P0.Velocity = Vec3D(0.0,0.0,0.0);
			add_particle(P0);
			break;
		}			
	}
		
	void setup() {
		// Problem parameters
		setName("statisticsTest");
		//auto_number();
		//fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
		//dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
		
		//time stepping
		setTimeStep(1e-44);
		setTimeMax(0);
	 
		//particle radii
		setInflowParticleRadius(.5);
		setFixedParticleRadius(getInflowParticleRadius());
		setRoughBottomType(MULTILAYER);

		//particle properties
		setDensity(6/pi);
		setStiffness(100);
		setSlidingStiffness(2.0/7.0*getStiffness());
		setDissipation(25.0);
		setSlidingDissipation(0.0);
		setSlidingFrictionCoefficient(0.5);
		
		//chute properties
		setChuteAngle(45.0, 1.0);
		setChuteLength(20);
		setChuteWidth(10);
		set_H(20);
			
		//output parameters
		setSaveCount(25);
		//setXBallsColourMode(7);
		//setXBallsVectorScale(1);
		
	}

	//set approximate height of flow
	void set_H(double H) {setZMax(H);}

	void printTime() const {
		//~ cout << "t=" << setprecision(3) << left << setw(6) << t 
			//~ << ", tmax=" << setprecision(3) << left << setw(6) << tmax
			//~ << endl;
	}

};

int main(int argc, char *argv[])
{
	ChutePeriodic problem;
	problem.setup();
	problem.solve();
	
	cout << "Now do statistics" << endl;
	StatisticsVector<Z> stats("statisticsTest");
	stats.setN(200);
	stats.setCGWidth(0.0997355701);
	stats.doTimeAverage(true);
	stats.setZMinStat(-1);
	stats.setZMaxStat(3);
	stats.setZMaxStat(5);
	stats.setDoPeriodicWalls(false);
	//stats.verbose();
	stats.set_infiniteStressForFixedParticles(true);
	stats.readStatArguments(argc+1, argv-1);
	stats.statistics_from_fstat_and_data();
}
