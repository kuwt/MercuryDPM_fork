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

	void actionsBeforeTimeStep(){
		static bool done = false;
		if (t>.8*tmax && !done) {cout<< "decreased save_count" << endl; dataFile.getSaveCount() /= 4; done = true;}
	};
		
	void setupInitialConditions()
	{

		Chute::setupInitialConditions();
		set_NWall(1);
		Walls[0].set(Vec3D( 0.0, 0.0, -1.0), 3.4*MaxInflowParticleRadius);
		set_NWallPeriodic(2);
		WallsPeriodic[0].set(Vec3D( 1.0, 0.0, 0.0), getXMin(), xmax);
		WallsPeriodic[1].set(Vec3D( 0.0, 1.0, 0.0), getYMin(), getYMax());
		
		particleHandler.set_StorageCapacity(particleHandler.getNumberOfObjects()+getChuteLength()*getChuteWidth()*getZMax());
		add_particles();

		cout << endl << "Status before solve:" << endl;
		write(std::cout,false);
		cout 
			<< "tc=" << getCollisionTime() 
			<< ", eps="	<< getRestitutionCoefficient()
			<< ", vmax=" << getMaximumVelocity()
			<< ", InflowHeight/getZMax()=" << getInflowHeight()/getZMax()
			<< endl << endl;
		timer.set(t,tmax);
	}

	void add_particles() 
	{
		set_HGRID_num_buckets_to_power(particleHandler.getStorageCapacity());
		HGridActionsBeforeTimeLoop();
		HGridActionsBeforeTimeStep();
		InflowHeight = getZMax();
		setZMax(1.2*getZMax());
		
		writeRestartFile();
		//try to find new insertable particles
		while (particleHandler.getNumberOfObjects()<Particles.capacity()){
			create_inflow_particle();
			if (IsInsertable(P0)) {
				num_created++;
			} else InflowHeight += .0001* MaxInflowParticleRadius;
		}
		set_HGRID_num_buckets_to_power();
	}
		
	void create_inflow_particle()
	{
		P0.Radius = MaxInflowParticleRadius;
		P0.computeMass(Species);
		
		P0.Position.X = random.get_RN(getXMin()+2.0*P0.Radius,getXMax());
		P0.Position.Y = random.get_RN(getYMin()+2.0*P0.Radius,getYMax());
		P0.Position.Z = random.get_RN(getZMin()+2.0*P0.Radius,getInflowHeight());
		P0.Velocity = Vec3D(0.0,0.0,0.0);
	}

	void setup() {
		// Problem parameters
		setName("silbert");
		//auto_number();
		fStatFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
		dataFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
		restartFile.setFileType(FileType::MULTIPLE_FILES_PADDED);
		
		//time stepconstants::ping
		setTimeStep(1e-4);
		setTimeMax(2e7*getTimeStep());
		setSaveCount(1e6); //save every 100 time units

		//particle radii
		setInflowParticleRadius(.5);
		setFixedParticleRadius(getInflowParticleRadius());
		setRoughBottomType(MULTILAYER);

		//particle properties
		setDensity(6/constants::pi);
		setStiffness(2e5);
		setSlidingStiffness(2.0/7.0*getStiffness());
		setDissipation(25.0);
		setSlidingDissipation(0);
		setSlidingFrictionCoefficient(0.5);
		
		//chute properties
		setChuteAngle(24.0, 1.0);
		setChuteLength(20);
		setChuteWidth(10);
		set_H(20);
			
		//output parameters
		//setXBallsColourMode(7);
		//setXBallsVectorScale(1);
	}

	//set approximate height of flow
	void set_H(Mdouble H) {setZMax(H);}

	void printTime() const {
		cout << "t=" << setprecision(3) << left << setw(6) << t 
			<< ", tmax=" << setprecision(3) << left << setw(6) << tmax
			<< ", N=" << setprecision(3) << left << setw(6) << particleHandler.getNumberOfObjects()
			//<< ", time left=" << setprecision(3) << left << setw(6) << timer.getTime2Finish(t)
			<< ", finish by " << setprecision(3) << left << setw(6) << timer.getFinishTime(t)
			<< endl;
		cout.flush();
	}

	Time2Finish timer;
};

int main(int argc, char *argv[])
{
	ChutePeriodic problem;
	problem.setup();
	problem.readArguments(argc, argv);
	//~ problem.set_counter(5);
	//~ problem.setXMax(40);
	//~ problem.setYMax(20);
	//problem.writeRestartFile();
	problem.solve();
	//problem.writeRestartFile();
}
